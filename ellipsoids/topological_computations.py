import json
import os
import time
import typing
from typing import Optional
from enum import Enum
import logging

import numpy as np
import gudhi as gd
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree
from scipy import spatial
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar, brentq
from scipy.spatial import Delaunay, Voronoi
from multiprocessing import Pool
from multiprocessing import cpu_count

from ellipsoids.common import (
    AxesMode,
    Dataset,
    EllipsoidBuildSpec,
    LengthsMode,
    LengthsNormalization,
    Parameters,
    EllipsoidParameters,
    ComplexSubtype,
    Results,
    EllipsoidResults,
    Ellipsoid,
    ConversionType,
    convert,
    spherisize_axes,
    convert,
    ConversionType
)

logger = logging.getLogger(__name__)


############################
## NOT EVERYTHING WAS PORTED TO THE FASTER ALGO:
## - ALPHA is currently broken
## - there are functions here that are not used
############################


def normalize_lengths(lengths: np.ndarray,
                      normalization: LengthsNormalization = LengthsNormalization.MAX_TO_1,
                      min_length: float = 1e-4) -> np.ndarray:
    '''
    normalize axes lengths for consistency
    '''
    lengths = np.asarray(lengths, dtype=float)
    if normalization is LengthsNormalization.NONE:
        result = lengths
    elif normalization is LengthsNormalization.MAX_TO_1:
        m = float(np.max(lengths))
        if m <= 0:
            raise ValueError("Cannot normalize lengths: max(lengths) must be > 0.")
        result = lengths / m
    else:
        raise ValueError(f"Unknown LengthsNormalization: {normalization}")

    return np.maximum(result, min_length)


def normalize_axes(axes_per_point: np.ndarray) -> np.ndarray:
    """
    normalize vectors to unit lengths for consistency
    """
    norms = np.linalg.norm(axes_per_point, axis=-1, keepdims=True)
    return axes_per_point / norms

    

def _compute_axes(
        points: np.ndarray,
        spec: EllipsoidBuildSpec,
) -> tuple[np.ndarray, Optional[np.ndarray]]:
    
    match spec.axes_mode:
        
        case AxesMode.PCA:
            assert spec.nbhd_size is not None

            axes = np.zeros(shape=[len(points), len(points[0]), len(points[0])])
            pca_singulars_per_point = np.zeros(shape=[len(points), len(points[0])])

            kdTree = spatial.KDTree(points)
            _, nbhd_idx = kdTree.query(points, spec.nbhd_size)
            nbhds = points[nbhd_idx]

            for i in range(points):
                pca_axes_point, pca_lengths_pointwise = pca_per_point(points[i], nbhds[i])
                pca_singulars_per_point[i] = np.asarray(pca_lengths_pointwise)
                axes[i] = np.asarray(pca_axes_point)
            return axes, pca_singulars_per_point
        
        case AxesMode.GIVEN:
            assert spec.axes is not None
            return spec.axes, None
        
        case _:
            raise RuntimeError(f"Unknown axes_mode: {spec.axes_mode}")
    

def _compute_axes_lengths(
        points: np.ndarray,
        spec: EllipsoidBuildSpec,
        pca_singulars_per_point: Optional[np.ndarray] = None
        ) -> np.ndarray:

    match spec.lengths_mode:
        case LengthsMode.RATIO:
            assert spec.axes_ratios is not None
            return np.tile(spec.axes_ratios, (len(points),1))
        case LengthsMode.GIVEN:
            assert spec.axes_lengths is not None
            return np.asarray(spec.axes_lengths)
        case LengthsMode.PCA_SINGULAR_RATIO:
            if pca_singulars_per_point is None:
                raise ValueError("pca_singulars_per_point not calculated.")
            return pca_singulars_per_point
        case _:
            raise RuntimeError("Invalid mode combination reached.")


        
def build_ellipsoids(
        points: np.ndarray,
        axes: np.ndarray,
        axes_lengths,
        normalization,
        spec: EllipsoidBuildSpec
) -> list[Ellipsoid]:

    points = np.asarray(points, dtype=float)
    if points.ndim != 2:
        raise ValueError("points must have shape (N, d)")
    n_points, ambient_dim = points.shape

    spec.validate(n_points=n_points, dim=ambient_dim)

    axes, pca_singulars_per_point = _compute_axes(points, spec)
    axes_lengths = _compute_axes_lengths( points, spec, pca_singulars_per_point )
        
    axes_lengths_normalized = np.asarray(
        [normalize_lengths(axes_lengths[i], spec.lengths_normalization) for i in range(n_points)],
        dtype = float)
    axes_normalized = np.asarray(
        [normalize_axes(axes[i]) for i in range(n_points)],
        dtype= float)

    return [
        Ellipsoid(points[i], axes_normalized[i], axes_lengths_normalized[i])
        for i in range(n_points)
        ]



def pca_per_point(
        center: np.ndarray,
        nbhd_pts: np.ndarray
) -> list[np.ndarray]:

    pca = PCA(n_components=len(center))
    pca.fit(nbhd_pts)
    pca_vectors = pca.components_

    pca_singular_values = pca.singular_values_

    if np.any(np.isnan(pca_singular_values)):
        raise ValueError("NaNs in PCA singular values during ellipsoid fitting. \
            There might be some collinear points in the dataset.")

    return [pca_vectors, pca_singular_values]


def pca_axes_and_lengths(
        points: np.ndarray,
        nbhd_size: int
) -> tuple[np.ndarray, np.ndarray]:

    n_points, dim = points.shape 
    axes = np.zeros((n_points, dim, dim))
    axes_lengths = np.zeros((n_points, dim))
    
    kdTree = spatial.KDTree(points)
    _, nbhd_idx = kdTree.query(points, nbhd_size)
    nbhds = points[nbhd_idx]

    for i in range(n_points):
        axes[i], axes_lengths[i] = pca_per_point(points[i], nbhds[i])
    return axes, axes_lengths
    


# legacy
def axes_lengths_from_axes_ratios(axes_ratios):
    axes_ratios = np.asarray(axes_ratios)
    if np.any(axes_ratios <= 0):
        raise ValueError("Invalid axes_ratios: axes_ratios must contain positive values.")
    return axes_ratios / axes_ratios[0]


def fit_ellipsoid(center: np.ndarray,
                  nbhd_pts: np.ndarray,
                  axes_ratios: Optional[np.ndarray] = None,
                  min_axes_length: float = 1e-4) -> Ellipsoid:
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=len(center))
    pca.fit(nbhd_pts)
    axes = pca.components_

    if axes_ratios is None:
        singular_values = pca.singular_values_

        if np.any(np.isnan(singular_values)):
            raise ValueError("NaNs in PCA singular values during ellipsoid fitting. \
             There might be some collinear points in the dataset.")
            # logger.warning("NaNs in the PCA axes lengths - setting all axes_lengths to 1.")
            # axes_lengths = np.ones(len(axes))

        axes_lengths = singular_values / singular_values[0]
        axes_lengths = np.maximum(axes_lengths, min_axes_length)

    else:
        axes_lengths = axes_lengths_from_axes_ratios(axes_ratios)
        # axes_ratios = np.asarray(axes_ratios)
        # if np.any(axes_ratios <= 0):
        #     raise ValueError("Invalid axes_ratios: axes_ratios may not contain a zero.")
        # axes_lengths = axes_ratios / axes_ratios[0] # r determines the long axis (normalising the long axis to 1)

    return Ellipsoid(center, axes, axes_lengths)



def fit_ellipsoids(
        points: np.ndarray,
        nbhd_size: int,
        axes_ratios: Optional[np.ndarray] = None
        ) -> list[Ellipsoid]:
    logger.info('Creating KD tree... ')
    kdTree = spatial.KDTree(points)
    logger.info('KD tree created.')
    logger.info('Fitting ellipsoids... ')

    if len(points) < nbhd_size:
        raise ValueError(f"The neighbourhood size {nbhd_size=} \
         is larger than the total number of points {len(points)=}.")

    _,nbhd_idx = kdTree.query(points, nbhd_size)
    nbhds = points[nbhd_idx]
    ellipsoid_list \
        = [fit_ellipsoid(point, nbhd, axes_ratios)
           for point, nbhd in zip(points, nbhds)]
    logger.info('Fitting ellipsoids completed.')

    return ellipsoid_list



# def spherisize(axes_lengths: np.ndarray, s:float=0):
#     '''
#     Returns axes lengths of a "spherisized" ellipsoid.
#     For s=0, the function will output axes_lengths and for s=1 all elements of axes_lengths
#     will be equal to the longest one.

#     For example, for axes_ratios=np.array([3,2,1]), we get:
#     - s=0: np.array([3,2,1])
#     - s=0.5: np.array([3,2.5,2])
#     - s=1: np.array([3,3,3])
#     '''

#     major_axis = np.max(axes_lengths)
#     spherisized_axes_lengths = np.zeros(np.size(axes_lengths))
#     for idx, axis in enumerate(axes_lengths):
#         spherisized_axes_lengths[idx] = s*major_axis + (1-s)*axis
#     return spherisized_axes_lengths



# def scale_to_01(x: float, min: float = 0, max: float = 1):
#     '''
#     Scales the input x in the interval [min, max] to the interval [0,1]
#     If x is bigger than max, returns 1.
#     If x is smaller than min, returns 0.
#     '''
#     if x > max:
#         return 1
#     elif x < min:
#         return 0
#     else:
#         return (x-min)/(max-min)



# def spherisize_axes(axes_lengths: np.ndarray, r: float, r_spherisize: float = 10):
#     '''
#     At r=0 should have axes_ratio from the start.
#     At r=r_spherisize, should get balls.
#     '''

#     return spherisize(axes_lengths, scale_to_01(r,max=r_spherisize))



def K(s, lambdas, v_squared, r):
    ''' Auxiliary function needed in ellipsoidIntersection
    '''
    return 1.-(1./r**2)*np.sum(v_squared*((s*(1.-s))/(1.+s*(lambdas-1.))))



def _ellipsoid_intersection(center_1: np.ndarray, axes_lengths_1: np.ndarray, axes_1: np.ndarray,
                            center_2: np.ndarray, axes_lengths_2: np.ndarray, axes_2: np.ndarray,
                            r: float):
    ''' Checks whether two ellipsoids at filtration level r intersect.
    The method is from https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect
    i.e. this paper: https://tisl.cs.toronto.edu/publication/201207-fusion-kalman_filter_fault_detection/fusion12-kalman_filter_fault_detection.pdf
    :return: true or false
    '''
    Sigma_A = axes_1.T @ np.diag(axes_lengths_1**2) @ axes_1
    Sigma_B = axes_2.T @ np.diag(axes_lengths_2**2) @ axes_2
    mu_A = center_1
    mu_B = center_2

    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    res = minimize_scalar(K,
                          bracket=[0.0, 0.5, 1.0],
                          args=(lambdas, v_squared, r))
    return res.fun >= 0



def ellipsoid_intersection(ellipsoid_1: Ellipsoid,
                           ellipsoid_2: Ellipsoid,
                           r: float,
                           r_spherisize: float = np.inf):
    ''' Checks whether ellipsoid_1 and ellipsoid_2 at the filtration level r intersect.
    If r_spherisize is not infinite, the axes_lengths will be linearly adapted so that
    for filtrations above spherisize_filtration, the ellipsoids become spheres
    '''

    if ellipsoid_1 == ellipsoid_2:
        return True

    axes_lengths_1 = spherisize_axes(ellipsoid_1.axes_lengths, r, r_spherisize)
    axes_lengths_2 = spherisize_axes(ellipsoid_2.axes_lengths, r, r_spherisize)

    return _ellipsoid_intersection(ellipsoid_1.center, axes_lengths_1, ellipsoid_1.axes,
                                   ellipsoid_2.center, axes_lengths_2, ellipsoid_2.axes,
                                   r)



def ellipsoid_intersection_cached(
        ellipsoid_1: Ellipsoid,
        ellipsoid_2: Ellipsoid,
        r: float,
        r_spherisize: float = np.inf
        ):
    ''' Checks whether ellipsoid_1 and ellipsoid_2 at the filtration level r intersect.
    If r_spherisize is not infinite, the axes_lengths will be linearly adapted so that
    for filtrations above spherisize_filtration, the ellipsoids become spheres
    '''

    if ellipsoid_1 == ellipsoid_2:
        return True

    Sigma_A = ellipsoid_1.get_sigma(r, r_spherisize)
    Sigma_B = ellipsoid_2.get_sigma(r, r_spherisize)
    mu_A = ellipsoid_1.center
    mu_B = ellipsoid_2.center

    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    res = minimize_scalar(K,
                          bracket=[0.0, 0.5, 1.0],
                          args=(lambdas, v_squared, r))
    return res.fun >= 0



def get_max_axes_ratio(ellipsoid: Ellipsoid):
    max_axis_length = max(ellipsoid.axes_lengths)
    min_axis_length = min(ellipsoid.axes_lengths)

    return max_axis_length / min_axis_length



# def find_intersection_radius(
#         ellipsoid_1: Ellipsoid,
#         ellipsoid_2: Ellipsoid,
#         threshold = 0.001,
#         epsilon = 0.001,
#         r_spherisize: float = np.inf,
#         use_cache: bool = True
#         ):

#     dist = np.linalg.norm(ellipsoid_1.center - ellipsoid_2.center)
#     max_axes_ratio = max(get_max_axes_ratio(ellipsoid_1), get_max_axes_ratio(ellipsoid_2))
#     lower_bound_r = (dist / 2) - epsilon              # maximum filtration at which ellipsoids can not intersect
#     upper_bound_r = dist/2 * max_axes_ratio + epsilon # minimum filtration at which ellipsoids can intersect
#     r = (upper_bound_r - lower_bound_r)/2

#     intersection_fn = (
#         ellipsoid_intersection_cached if use_cache
#         else ellipsoid_intersection
#         )

#     while True:
#         if intersection_fn(ellipsoid_1, ellipsoid_2, r, r_spherisize):
#             upper_bound_r = r
#         else: lower_bound_r = r

#         if (upper_bound_r - lower_bound_r) < threshold:
#             return r
#         else: r = (upper_bound_r + lower_bound_r)/2


def _g_max(Sigma_A: np.ndarray, Sigma_B: np.ndarray, delta_mu: np.ndarray) -> float:
    """max over s in [0,1] of g(s) = sum_i v_i^2 s(1-s)/(1+s(lambda_i-1)).

    The ellipsoids (with these fixed shape matrices) scaled by r intersect
    iff r^2 >= g_max, so sqrt(g_max) is their exact intersection radius.
    """
    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, delta_mu) ** 2
    res = minimize_scalar(
        lambda s: -np.sum(v_squared * s * (1. - s) / (1. + s * (lambdas - 1.))),
        bounds = (0.0, 1.0),
        method = 'bounded',
    )

    return max(-res.fun, 0.0)



    

def find_intersection_radius_closed_form(
        ellipsoid_1: Ellipsoid,
        ellipsoid_2: Ellipsoid,
        ) -> float:
    """Exact intersection radius for r_spherisize == np.inf. No bisection."""
    if ellipsoid_1 == ellipsoid_2:
        return 0.0
    return float(np.sqrt(_g_max(
        ellipsoid_1.get_sigma(0.0, np.inf),
        ellipsoid_2.get_sigma(0.0, np.inf),
        ellipsoid_1.center - ellipsoid_2.center)))


def find_intersection_radius_spherisized(
        ellipsoid_1: Ellipsoid,
        ellipsoid_2: Ellipsoid,
        r_spherisize: float,
        xtol: float = 1e-9,
) -> float:
    """Intersection radius for finite r_spherisize.

    The shape matrices depend on r, so we root-find h(r) = r^2 - g_max(r)
    with Brent's method. Bracket: ellipsoids sit inside balls of radius
    r * max_axis_length, so no intersection below dist/(2*max_len); and the
    unspherisized closed-form radius is an upper bound because spherisizing
    only enlarges each ellipsoid. Since ellipsoids grow monotonically with r,
    h has a single sign change and brentq converges superlinearly.
    """
    if ellipsoid_1 == ellipsoid_2:
        return 0.0
    delta_mu = ellipsoid_1.center - ellipsoid_2.center
    dist = np.linalg.norm(delta_mu)
    if dist == 0.0:
        return 0.0

    def h(r):
        return r * r - _g_max(ellipsoid_1.get_sigma(r, r_spherisize),
                              ellipsoid_2.get_sigma(r, r_spherisize),
                              delta_mu)

    amax = max(np.max(ellipsoid_1.axes_lengths), np.max(ellipsoid_2.axes_lengths))
    lower = dist / (2.0 * amax)
    upper = max(find_intersection_radius_closed_form(ellipsoid_1, ellipsoid_2),
                lower) * (1.0 + 1e-9) + 1e-12

    if h(lower) >= 0:  # touching exactly at the theoretical minimum
        return float(lower)
    while h(upper) < 0:  # numerical safety net; essentially never loops
        upper *= 1.5

    return float(brentq(h, lower, upper, xtol=xtol))




def find_intersection_radius(
        ellipsoid_1: Ellipsoid,
        ellipsoid_2: Ellipsoid,
        r_spherisize: float = np.inf,
        ) -> float:
    """Dispatcher depending on r_spherisize: closed form (inf) or Brent root-find (finite)."""
    if r_spherisize == np.inf:
        return find_intersection_radius_closed_form(ellipsoid_1, ellipsoid_2)
    return find_intersection_radius_spherisized(ellipsoid_1, ellipsoid_2, r_spherisize)



def wrapper_find_intersection_radius(
        ellipsoid_1,
        ellipsoid_2,
        r_spherisize: float,
        # use_cache: bool,
        ):
    '''
    Wrapper that allows for an additional (non-keyword) argument.
    Necessary to divide tasks into multiple cores for multiprocessing.
    '''
    return find_intersection_radius(
        ellipsoid_1,
        ellipsoid_2,
        r_spherisize=r_spherisize)
        # use_cache=use_cache)


def generate_ellipsoid_list(
        points: np.ndarray,
        axes: np.ndarray,
        axes_ratios: Optional[np.ndarray] = None,
        ):

    ellipsoid_list = []

    normalized_axes = np.array([
        [
            axis/np.linalg.norm(axis) for axis in axis_tuple
        ]
        for axis_tuple in axes
    ])

    if axes_ratios is None:
        axes_ratios = np.asarray([[2,1]] * len(points))

    if len(axes_ratios) != len(points):
        axes_ratios = np.asarray([axes_ratios] * len(points))


    normalized_axes_ratios = np.array([
        abs(axis_ratio) / max(abs(axis_ratio)) for axis_ratio in axes_ratios 
    ])

    for idx, point in enumerate(points):
        ellipsoid_list.append(Ellipsoid(point, normalized_axes[idx], normalized_axes_ratios[idx]))
   
    return ellipsoid_list
    
def _batch_insert_vertices(simplex_tree, n_points: int):
    """Insert vertices 0..n-1 at filtration 0 in one call when possible."""
    if hasattr(simplex_tree, "insert_batch"):
        simplex_tree.insert_batch(
            np.arange(n_points, dtype=np.int64).reshape(1, -1),
            np.zeros(n_points))
    else:
        [simplex_tree.insert([i], 0) for i in np.arange(n_points)]

def _batch_insert_edges(simplex_tree, pairs, filtrations):
    """Insert edges in one call via gudhi's insert_batch when available.

    insert_batch expects vertex_array of shape (2, n) for edges and a length-n
    filtration array; repeated simplices keep the smallest filtration and
    missing faces are inserted automatically (harmless here since all vertices
    are pre-inserted at filtration 0). Falls back to per-edge insert on old
    gudhi versions.
    """
    pairs = np.asarray(pairs, dtype=np.int64)
    filtrations = np.asarray(filtrations, dtype=np.float64)
    if len(pairs) == 0:
        return
    if hasattr(simplex_tree, "insert_batch"):
        simplex_tree.insert_batch(np.ascontiguousarray(pairs.T), filtrations)
    else:
        for pair, f in zip(pairs, filtrations):
            simplex_tree.insert(list(pair), float(f))



def batch_intersection_radii(
        Sigma_A: np.ndarray,  # (N, d, d)
        Sigma_B: np.ndarray,  # (N, d, d)
        delta_mu: np.ndarray,  # (N, d)
        n_golden_iters: int = 60,
        base_chunk_size=400_000,
) -> np.ndarray:
    """Vectorized closed form over N pairs of points (r_spherisize == inf only).

    Cholesky whitening turns the generalized eigenproblem into a batched
    np.linalg.eigh; g(s) is unimodal on [0,1], so a vectorized golden-section
    search maximizes all pairs simultaneously.
    """

    _inv_phi = (np.sqrt(5.0) - 1.0) / 2.0

    n_total = len(Sigma_A)
    dim = Sigma_A.shape[1]
    chunk_size = max(1, int(base_chunk_size // (dim * dim)))
    if n_total > chunk_size:
        out = np.empty(n_total)
        for s in range(0, n_total, chunk_size):
            e = s + chunk_size
            out[s:e] = batch_intersection_radii(
                Sigma_A[s:e], Sigma_B[s:e], delta_mu[s:e],
                n_golden_iters=n_golden_iters, base_chunk_size=base_chunk_size)
        return out

    L = np.linalg.cholesky(Sigma_B)
    L_inv = np.linalg.inv(L)
    M = L_inv @ Sigma_A @ np.transpose(L_inv, (0, 2, 1))
    lambdas, Q = np.linalg.eigh(M)
    w = np.einsum('nij,nj->ni', L_inv, delta_mu)
    v2 = np.einsum('nji,nj->ni', Q, w) ** 2

    def g(s):
        s_ = s[:, None]
        return np.sum(v2 * s_ * (1. - s_) / (1. + s_ * (lambdas - 1.)), axis=1)

    n = len(Sigma_A)
    a = np.zeros(n)
    b = np.ones(n)
    c = b - _inv_phi * (b - a)
    d = a + _inv_phi * (b - a)
    gc, gd = g(c), g(d)
    for _ in range(n_golden_iters):
        keep_left = gc > gd
        b = np.where(keep_left, d, b)
        a = np.where(keep_left, a, c)
        c = b - _inv_phi * (b - a)
        d = a + _inv_phi * (b - a)
        gc, gd = g(c), g(d)

    return np.sqrt(np.maximum(g((a + b) / 2.0), 0.0))



def generate_ELLIPSOID_RIPS_simplex_tree(
        ellipsoid_list: list[Ellipsoid],
        r_spherisize: float = np.inf,
        use_cache: bool = True,
        radius_conversion_fn = lambda f: f,
        ):

    logger.info('Calculating ellipsoid simplex tree... ')

    new = True

    if new:
        points = np.asarray([e.center for e in ellipsoid_list])
        simplexTree = gd.SimplexTree()
        _batch_insert_vertices(simplexTree, len(points))

        iu = np.triu_indices(len(points), k=1)
        pairs = np.stack(iu, axis=1)

        if len(pairs) == 0:
            return simplexTree

        cpuCores = int(os.environ.get("SLURM_NTASKS", cpu_count()))

        if r_spherisize == np.inf:

            sigmas = np.stack([e.get_sigma(0.0, np.inf) for e in ellipsoid_list])
            centers = np.stack([e.center for e in ellipsoid_list])
            i_idx, j_idx = pairs[:, 0], pairs[:, 1]
            radii = batch_intersection_radii(
                sigmas[i_idx], sigmas[j_idx], centers[i_idx] - centers[j_idx])
            task_cpu_times = []
        else:
            tasks = [(ellipsoid_list[i], ellipsoid_list[j], r_spherisize)
                    for i, j in pairs]
            with Pool(cpuCores) as p:
                radii = list(p.starmap(wrapper_find_intersection_radius, tasks))

        filtrations = np.asarray([radius_conversion_fn(r) for r in radii],
                                dtype=np.float64)
        _batch_insert_edges(simplexTree, pairs, filtrations)

        return simplexTree, ellipsoid_list

    else:
        points = np.asarray([e.center for e in ellipsoid_list])
        simplexTree = gd.SimplexTree()
        [simplexTree.insert([i],0) for i in np.arange(len(points))]

        pairs = np.array([[i,j] for i in np.arange(len(points)) for j in np.arange(i+1,len(points))])
        tasks = zip([ellipsoid_list[i] for i in pairs[:,0]], \
                    [ellipsoid_list[j] for j in pairs[:,1]])
        tasks = [(*x, r_spherisize, use_cache) for x in tasks]

        cpuCores = int(os.environ.get("SLURM_NTASKS", cpu_count()))
        with Pool(cpuCores) as p:
            radii = list(p.starmap(wrapper_find_intersection_radius, tasks))

        # filtrations = [convert(r, radius_conversion_type) for r in radii]
        filtrations = [radius_conversion_fn(r) for r in radii]
        [simplexTree.insert(pair,f) for pair, f in zip(pairs,filtrations)]

        logger.info('Ellipsoid simplex tree calculated.')
        return simplexTree



def adjacent_delaunay_vertices(triangulation, vertex_index):
    adjacent_vertices = set()

    for simplex in triangulation.simplices:
        if vertex_index in simplex:
            adjacent_vertices.update(simplex[simplex != vertex_index])

    return sorted(adjacent_vertices)



def generate_ELLIPSOID_ALPHA_simplex_tree(
        ellipsoid_list: list[Ellipsoid],
        r_spherisize: float,
        use_cache: bool,
        radius_conversion_fn = lambda v: convert(v, ConversionType.NONE), # default: r, not as in GUDHI alpha
        ):

    logger.info("Calculating ELLIPSOID ALPHA simplex tree... ")

    points = np.asarray([e.center for e in ellipsoid_list])
    simplex_tree = gd.SimplexTree()
    [simplex_tree.insert([i],0) for i in np.arange(len(points))]

    delaunay = Delaunay(points)
    for vertex_index, _ in enumerate(delaunay.points):
        adjacent_vertices = adjacent_delaunay_vertices(delaunay, vertex_index)

        for adjacent_vertex in adjacent_vertices:
            vertex_i = vertex_index
            vertex_j = adjacent_vertex

            if not simplex_tree.find([vertex_i, vertex_j]):
                intersection_radius = find_intersection_radius(
                    ellipsoid_list[vertex_i],
                    ellipsoid_list[vertex_j],
                    r_spherisize=r_spherisize,
                    use_cache=use_cache,
                    )
                # filtration = convert(intersection_radius, radius_conversion_type)
                filtration = radius_conversion_fn(intersection_radius)
                simplex_tree.insert([vertex_i,vertex_j], filtration)

    logger.info("Simplex tree calculated.")
    return simplex_tree



def expand_simplex_tree(simplex_tree, expansion_dim=2):
    logger.info("Expanding the simplex tree...")
    simplex_tree.expansion(expansion_dim)
    logger.info("Simplex tree expanded.")
    return simplex_tree



def collapse_edges(simplex_tree):
    logger.info("Collapsing edges...")
    simplex_tree.collapse_edges()
    logger.info("Edges collapsed.")
    return simplex_tree



def generate_BALL_RIPS_simplex_tree(
        points,
        max_dimension=1,
        filtration_conversion_fn = lambda v: convert(v, ConversionType.RIPS_TO_RADIUS),
        ):

    logger.info("Creating the Rips complex... ")
    rips_complex = gd.RipsComplex(points=points)
    logger.info("Rips complex created.")

    logger.info("Creating the Rips simplex tree...")
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dimension)
    logger.info("Simplex tree created.")

    map_filtration(simplex_tree, filtration_conversion_fn)

    return simplex_tree



def double_filtrations(simplex_tree):
    for simplex in simplex_tree.get_simplices():
        simplex_set = simplex[0]  # Get the simplex (set of vertices)
        filtration_value = simplex[1]  # Get the current filtration value
        simplex_tree.assign_filtration(simplex_set, 2*filtration_value)
    return simplex_tree



def map_filtration(st: gd.SimplexTree, f):
    for simplex, val in st.get_filtration():
        st.assign_filtration(simplex, f(val))
    st.make_filtration_non_decreasing()


def generate_BALL_ALPHA_simplex_tree(
        points,
        filtration_conversion_fn = lambda f: np.sqrt(f)
        ):
    logger.info("Creating the alpha complex...")
    alpha_complex = gd.AlphaComplex(points=points)
    logger.info("Alpha complex created.")

    logger.info("Creating the alpha simplex tree...")
    simplex_tree = alpha_complex.create_simplex_tree()
    logger.info("Alpha simplex tree created.")

    map_filtration(simplex_tree, filtration_conversion_fn)

    return simplex_tree



def calculate_barcode(simplex_tree):
    logger.info("Calculating the barcode of the expanded tree...")
    barcode = simplex_tree.persistence()
    logger.info("Barcode calculated.")
    return barcode



def maxFiltration(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    return max(splx[1] for splx in simplexList)



def set_max_bar_end(bar, max_bar_end):
    bar_end = bar[1][1]
    if bar_end != float('inf') and bar_end > max_bar_end:
        max_bar_end = bar_end
    return max_bar_end



def set_axes_ratios_to_tangent(axes_ratios: np.ndarray, dim: int, manifold_dim: int):
    ratio = axes_ratios[0]
    return pad_axes_ratios(np.repeat(ratio,manifold_dim), dim)



def pad_axes_ratios(axesRatios: np.ndarray, dim: int):
    ''' For high dimensional ellipsoids, it is enough for the user to specify 
    the first few axes. This function will set the remaining axes to 1.'''
    if dim > len(axesRatios):
        return np.pad(axesRatios, (0,dim - len(axesRatios)), constant_values=1) # creates an array of length dim
    else:
        return axesRatios[0:dim]



def calculate_BALL_Results(
        dataset: Dataset,
        parameters: Parameters
        ) -> Results:

    t0_simplex_tree = time.time()
    if parameters.complex_subtype == ComplexSubtype.RIPS:
        simplex_tree = generate_BALL_RIPS_simplex_tree(dataset.points)
    elif parameters.complex_subtype == ComplexSubtype.ALPHA:
        simplex_tree = generate_BALL_ALPHA_simplex_tree(dataset.points)
    else:
        exit(f"{parameters.complex_subtype} is an invalid complex type.")

    if parameters.collapse_edges:
        collapse_edges(simplex_tree)
    if parameters.expansion_dim > 1:
        expand_simplex_tree(simplex_tree, expansion_dim=parameters.expansion_dim)
    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = calculate_barcode(simplex_tree)
    t1_barcode = time.time()
    execution_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    results = Results()

    if parameters.save_simplex_tree: results.simplex_tree = simplex_tree
    results.barcode = barcode
    results.execution_time = execution_time

    return results



def calculate_ELLIPSOID_Results(
        dataset: Dataset,
        ellipsoid_parameters: EllipsoidParameters
        ) -> EllipsoidResults:

    ambient_dim = dataset.ambient_dim
    points = dataset.points

    if ellipsoid_parameters.axes_ratios is None:
        axes_ratios = None
    else:
        axes_ratios = set_axes_ratios_to_tangent(
            ellipsoid_parameters.axes_ratios,
            ambient_dim,
            ambient_dim-1)
        
    nbhd_size = ellipsoid_parameters.nbhd_size
    axes = ellipsoid_parameters.axes
    if axes is None:
        if nbhd_size is None:
            raise TypeError("If axes are not provided, nhbd_size must be.")
        ellipsoid_list = fit_ellipsoids(points, nbhd_size, axes_ratios)
    else:
        ellipsoid_list = generate_ellipsoid_list(points, axes, axes_ratios) 

    r_spherisize = ellipsoid_parameters.r_spherisize
    use_cache = ellipsoid_parameters.use_cache
    axes = ellipsoid_parameters.axes

    t0_simplex_tree = time.time()
    if (ellipsoid_parameters.complex_subtype == ComplexSubtype.ALPHA):
        simplex_tree = generate_ELLIPSOID_ALPHA_simplex_tree(
            ellipsoid_list = ellipsoid_list,
            r_spherisize=r_spherisize,
            use_cache=use_cache,
            )
    elif ellipsoid_parameters.complex_subtype == ComplexSubtype.RIPS:
        [simplex_tree, ellipsoid_list] = generate_ELLIPSOID_RIPS_simplex_tree(
            ellipsoid_list=ellipsoid_list,
            r_spherisize=r_spherisize,
            use_cache=use_cache,
            )
    else:
        raise ValueError(f"{ellipsoid_parameters.complex_subtype} is an invalid complex subtype.")

    if ellipsoid_parameters.collapse_edges:
        collapse_edges(simplex_tree)
    if ellipsoid_parameters.expansion_dim > 1:
        expand_simplex_tree(simplex_tree, expansion_dim=ellipsoid_parameters.expansion_dim)
    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = calculate_barcode(simplex_tree)
    t1_barcode = time.time()

    execution_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    results = EllipsoidResults()

    if ellipsoid_parameters.save_simplex_tree: results.simplex_tree = simplex_tree
    if ellipsoid_parameters.save_ellipsoid_list: results.ellipsoid_list = ellipsoid_list
    results.barcode = barcode
    results.execution_time = execution_time

    return results
