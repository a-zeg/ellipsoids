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
from scipy.optimize import minimize_scalar
from scipy.spatial import Delaunay, Voronoi
from multiprocessing import Pool
from multiprocessing import cpu_count

from ellipsoids.common import Dataset
from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import Results
from ellipsoids.common import EllipsoidResults
from ellipsoids.common import Ellipsoid
from ellipsoids.common import ConversionType
from ellipsoids.common import convert
from ellipsoids.common import spherisize_axes

logger = logging.getLogger(__name__)



def fit_ellipsoid(center: np.ndarray,
                  nbhd_pts: np.ndarray,
                  axes_ratios: Optional[np.ndarray] = None,
                  min_axes_length: Optional[float] = 1e-4) -> Ellipsoid:
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=len(center))
    pca.fit(nbhd_pts)
    axes = pca.components_


    if axes_ratios is None:
        singular_values = pca.singular_values_

        if np.any(np.isnan(singular_values)):
            logger.warning("NaNs in the PCA axes lengths - setting all axes_lengths to 1.")
            axes_lengths = np.ones(len(axes))

        else:
            axes_lengths = singular_values / singular_values[0]
            axes_lengths = np.maximum(axes_lengths, min_axes_length)

    else:
        if axes_ratios.all() != 0:
            axes_lengths = axes_ratios / axes_ratios[0] # r determines the long axis (normalising the long axis to 1)
        else:
            raise ValueError("Invalid axes_ratios: axes_ratios may not contain a zero.")

    return Ellipsoid(center, axes, axes_lengths)



def fit_ellipsoids(points: np.ndarray, neighbourhood_size: int, axes_ratios: Optional[np.ndarray] = None) -> list[Ellipsoid]:
    logger.info('Creating KD tree... ')
    kdTree = spatial.KDTree(points)
    logger.info('KD tree created.')
    logger.info('Fitting ellipsoids... ')

    if len(points) < neighbourhood_size:
        logger.warning('WARNING: the chosen neighbourhood size is too small. \
              Setting the neighbhourhood size to the total number of points.')
        neighbourhood_size = len(points)

    _,neighbourhood_idx = kdTree.query(points, neighbourhood_size)
    neighbourhoods = points[neighbourhood_idx]
    ellipsoid_list \
        = [fit_ellipsoid(point, neighbourhood, axes_ratios) for point,neighbourhood in zip(points, neighbourhoods)]
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



def find_intersection_radius(
        ellipsoid_1: Ellipsoid,
        ellipsoid_2: Ellipsoid,
        threshold = 0.001,
        epsilon = 0.001,
        r_spherisize: float = np.inf,
        use_cache: bool = True
        ):

    dist = np.linalg.norm(ellipsoid_1.center - ellipsoid_2.center)
    max_axes_ratio = max(get_max_axes_ratio(ellipsoid_1), get_max_axes_ratio(ellipsoid_2))
    lower_bound_r = (dist / 2) - epsilon              # maximum filtration at which ellipsoids can not intersect
    upper_bound_r = dist/2 * max_axes_ratio + epsilon # minimum filtration at which ellipsoids can intersect
    r = (upper_bound_r + lower_bound_r)/2

    intersection_fn = (
        ellipsoid_intersection_cached if use_cache
        else ellipsoid_intersection
        )

    while True:
        if intersection_fn(ellipsoid_1, ellipsoid_2, r, r_spherisize):
            upper_bound_r = r
        else: lower_bound_r = r

        if (upper_bound_r - lower_bound_r) < threshold:
            return r
        else: r = (upper_bound_r + lower_bound_r)/2



def wrapper_find_intersection_radius(ellipsoid_1, ellipsoid_2, r_spherisize: float, use_cache: bool):
    '''
    Wrapper that allows for an additional (non-keyword) argument.
    Necessary to divide tasks into multiple cores for multiprocessing.
    '''
    return find_intersection_radius(ellipsoid_1, ellipsoid_2, r_spherisize=r_spherisize, use_cache=use_cache)



def generate_ELLIPSOID_RIPS_simplex_tree(
        points: np.ndarray,
        nbhd_size: int,
        axes_ratios: Optional[np.ndarray] = None,
        r_spherisize: float = np.inf,
        use_cache: bool = True,
        ):
    ''' multiprocessing '''
    ''' Creates a simplex tree from the ellipsoids by adding an edge between each two points whose 
    corresponding ellipsoids intersect.
    :kdTree: KD tree of the initial dataset
    :ellipsoidList: list of ellipsoids (output of ??)
    :queryRadius:
    :filtrationValues:

    :return: gudhi.SimplexTree
    '''

    ellipsoidList = fit_ellipsoids(points, nbhd_size, axes_ratios)

    logger.info('Calculating ellipsoid simplex tree... ')

    simplexTree = gd.SimplexTree()
    [simplexTree.insert([i],0) for i in np.arange(len(points))]
    
    pairs = np.array([[i,j] for i in np.arange(len(points)) for j in np.arange(i+1,len(points))])
    tasks = zip([ellipsoidList[i] for i in pairs[:,0]], \
                [ellipsoidList[j] for j in pairs[:,1]])
    tasks = [(*x, r_spherisize, use_cache) for x in tasks]

    cpuCores = int(os.environ.get("SLURM_NTASKS", cpu_count()))
    with Pool(cpuCores) as p:
        radii = list(p.starmap(wrapper_find_intersection_radius, tasks))

    filtrations = [convert(r, ConversionType.RADIUS_TO_RIPS) for r in radii]
    [simplexTree.insert(pair,f) for pair, f in zip(pairs,filtrations)]

    logger.info('Ellipsoid simplex tree calculated.')
    return [simplexTree, ellipsoidList]



def adjacent_delaunay_vertices(triangulation, vertex_index):
    adjacent_vertices = set()

    for simplex in triangulation.simplices:
        if vertex_index in simplex:
            adjacent_vertices.update(simplex[simplex != vertex_index])

    return sorted(adjacent_vertices)



def generate_ELLIPSOID_ALPHA_simplex_tree(
        points: np.ndarray,
        nbhd_size: int,
        axes_ratios: Optional[np.ndarray],
        r_spherisize: float,
        use_cache: bool,
        ):
    ellipsoid_list: list[Ellipsoid] = fit_ellipsoids(points, nbhd_size, axes_ratios)

    logger.info("Calculating ELLIPSOID ALPHA simplex tree... ")

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
                filtration = convert(intersection_radius, ConversionType.RADIUS_TO_ALPHA)
                simplex_tree.insert([vertex_i,vertex_j], filtration)

    logger.info("Simplex tree calculated.")
    return [simplex_tree, ellipsoid_list]



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



def generate_BALL_RIPS_simplex_tree(points, max_dimension=1):

    logger.info("Creating the Rips complex... ")
    rips_complex = gd.RipsComplex(points=points)
    logger.info("Rips complex created.")

    logger.info("Creating the Rips simplex tree...")
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dimension)
    logger.info("Simplex tree created.")

    return simplex_tree



def double_filtrations(simplex_tree):
    for simplex in simplex_tree.get_simplices():
        simplex_set = simplex[0]  # Get the simplex (set of vertices)
        filtration_value = simplex[1]  # Get the current filtration value
        simplex_tree.assign_filtration(simplex_set, 2*filtration_value)
    return simplex_tree



def generate_BALL_ALPHA_simplex_tree(points):
    logger.info("Creating the alpha complex...")
    alpha_complex = gd.AlphaComplex(points=points)
    logger.info("Alpha complex created.")

    logger.info("Creating the alpha simplex tree...")
    simplex_tree = alpha_complex.create_simplex_tree()
    logger.info("Alpha simplex tree created.")

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



# def reduce_barcode(barcode, nBarsDim0 = 10, nBarsDim1 = 10, nBarsDim2 = 10):
#     # return only the first nBarsDimk bars in each dimension k
#     reduced_barcode = []
#     max_bar_end = 0
#     for bar in barcode:
#         if bar[0] == 0 and nBarsDim0 > 0:
#             reduced_barcode.append(bar)
#             nBarsDim0 = nBarsDim0 - 1
#             max_bar_end = set_max_bar_end(bar, max_bar_end)
#         elif bar[0] == 1 and nBarsDim1 > 0:
#             reduced_barcode.append(bar)
#             nBarsDim1 = nBarsDim1 - 1
#             max_bar_end = set_max_bar_end(bar, max_bar_end)
#         elif bar[0] == 2 and nBarsDim2 > 0:
#             reduced_barcode.append(bar)
#             nBarsDim2 = nBarsDim2 - 1
#             max_bar_end = set_max_bar_end(bar, max_bar_end)

#     return reduced_barcode, max_bar_end
    


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
    if ellipsoid_parameters.axes_ratios is None:
        axes_ratios = None
    else:
        axes_ratios = set_axes_ratios_to_tangent(
            ellipsoid_parameters.axes_ratios,
            ambient_dim,
            ambient_dim-1)
    points = dataset.points
    nbhd_size = ellipsoid_parameters.nbhd_size
    r_spherisize = ellipsoid_parameters.r_spherisize
    use_cache = ellipsoid_parameters.use_cache

    t0_simplex_tree = time.time()
    if (ellipsoid_parameters.complex_subtype == ComplexSubtype.ALPHA):
        simplex_tree, ellipsoid_list = generate_ELLIPSOID_ALPHA_simplex_tree(
            points=points,
            nbhd_size=nbhd_size,
            axes_ratios=axes_ratios,
            r_spherisize=r_spherisize,
            use_cache=use_cache,
            )
    elif ellipsoid_parameters.complex_subtype == ComplexSubtype.RIPS:
        [simplex_tree, ellipsoid_list] = generate_ELLIPSOID_RIPS_simplex_tree(
            points=points,
            nbhd_size=nbhd_size,
            axes_ratios=axes_ratios,
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
