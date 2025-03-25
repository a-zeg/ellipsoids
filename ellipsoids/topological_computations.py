import json
import os
import time
import typing

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





def fit_ellipsoid(center: np.ndarray, nbhd_pts: np.ndarray, axes_ratios: np.ndarray) -> Ellipsoid:
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=len(center))
    pca.fit(nbhd_pts)
    axes = pca.components_
    axes_lengths = pca.singular_values_

    if axes_ratios.all() != 0:
        axes_lengths = axes_ratios / axes_ratios[0] # r determines the long axis (normalising the long axis to 1)
        # axesLengths = axesRatios / axesRatios[-1] # alt: r determines the short axis
    else: 
        exit("Error: axes ratios contain a zero.")
    return Ellipsoid(center, axes, axes_lengths)



def fit_ellipsoids(points, neighbourhood_size, axes_ratios) -> list[Ellipsoid]:
    print('Creating KD tree... ', end='', flush=True)
    kdTree = spatial.KDTree(points)
    print('Done.')
    print('Fitting ellipsoids... ', end='', flush=True)

    if len(points) < neighbourhood_size:
        print('WARNING: the chosen neighbourhood size is too small. \
              Setting the neighbhourhood size to the total number of points.')
        neighbourhood_size = len(points)

    _,neighbourhood_idx = kdTree.query(points, neighbourhood_size)
    neighbourhoods = points[neighbourhood_idx]
    ellipsoid_list \
        = [fit_ellipsoid(point, neighbourhood, axes_ratios) for point,neighbourhood in zip(points, neighbourhoods)]
    print('Done.')

    return ellipsoid_list


def spherisize(axes_lengths: np.ndarray, s:float=0):
    '''
    Returns axes lengths of a "spherisized" ellipsoid.
    For s=0, the function will output axes_lengths and for s=1 all elements of axes_lengths
    will be equal to the longest one.

    For example, for axes_ratios=np.array([3,2,1]), we get:
    - s=0: np.array([3,2,1])
    - s=0.5: np.array([3,2.5,2])
    - s=1: np.array([3,3,3])
    '''

    major_axis = np.max(axes_lengths)
    spherisized_axes_lengths = np.zeros(np.size(axes_lengths))
    for idx, axis in enumerate(axes_lengths):
        spherisized_axes_lengths[idx] = s*major_axis + (1-s)*axis
    return spherisized_axes_lengths



def scale_to_01(x: float, min: float = 0, max: float = 1):
    '''
    Scales the input x in the interval [min, max] to the interval [0,1]
    If x is bigger than max, returns 1.
    If x is smaller than min, returns 0.
    '''
    if x > max:
        return 1
    elif x < min:
        return 0
    else:
        return (x-min)/(max-min)



def spherisize_axes(axes_lengths: np.ndarray, r: float, r_spherisize: float = 10):
    '''
    At r=0 should have axes_ratio from the start.
    At r=scale, should get balls.
    '''

    return spherisize(axes_lengths, scale_to_01(r,max=r_spherisize))



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
    Sigma_A = np.linalg.multi_dot([\
        np.transpose(axes_1),\
        np.diag(axes_lengths_1**2),\
        axes_1])
    Sigma_B = np.linalg.multi_dot([\
        np.transpose(axes_2),\
        np.diag(axes_lengths_2**2),\
        axes_2])
    mu_A = center_1
    mu_B = center_2

    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    res = minimize_scalar(K,
                          bracket=[0.0, 0.5, 1.0],
                          args=(lambdas, v_squared, r))
    return (res.fun >= 0)



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



def get_max_axes_ratio(ellipsoid: Ellipsoid):
    max_axis_length = max(ellipsoid.axes_lengths)
    min_axis_length = min(ellipsoid.axes_lengths)

    return max_axis_length / min_axis_length



def find_intersection_radius(ellipsoid_1: Ellipsoid,
                           ellipsoid_2: Ellipsoid,
                           threshold = 0.001,
                           epsilon = 0.001,
                           r_spherisize: float = np.inf):

    dist = np.linalg.norm(ellipsoid_1.center - ellipsoid_2.center)
    
    max_axes_ratio = max(get_max_axes_ratio(ellipsoid_1),
                         get_max_axes_ratio(ellipsoid_2))

    maxNonIntersectionFiltration = (dist / 2) - epsilon
    minIntersectionFiltration = dist/2 * max_axes_ratio + epsilon

    r = (minIntersectionFiltration - maxNonIntersectionFiltration)/2

    while True:
        if ellipsoid_intersection(ellipsoid_1, ellipsoid_2, r, r_spherisize):
            minIntersectionFiltration = r
        else: maxNonIntersectionFiltration = r

        if (minIntersectionFiltration - maxNonIntersectionFiltration) < threshold:
            return 2*r # 2r so that it's comparable to Rips
        else: r = (minIntersectionFiltration + maxNonIntersectionFiltration)/2



def wrapper_find_intersection_radius(ellipsoid_1, ellipsoid_2, r_spherisize: float):
    '''
    Wrapper that allows for an additional (non-keyword) argument.
    Necessary to divide tasks into multiple cores for multiprocessing.
    '''
    return find_intersection_radius(ellipsoid_1, ellipsoid_2, r_spherisize = r_spherisize)



def generate_ellipsoid_simplex_tree(points: np.ndarray,
                                  nbhd_size: int,
                                  axes_ratios: np.ndarray,
                                  r_spherisize: float = np.inf):
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

    print('Calculating ellipsoid simplex tree... ', end='', flush=True)

    simplexTree = gd.SimplexTree()
    [simplexTree.insert([i],0) for i in np.arange(len(points))]
    
    pairs = np.array([[i,j] for i in np.arange(len(points)) for j in np.arange(i+1,len(points))])
    tasks = zip([ellipsoidList[i] for i in pairs[:,0]], \
                [ellipsoidList[j] for j in pairs[:,1]])
    tasks = [(*x, r_spherisize) for x in tasks]

    cpuCores = int(os.environ.get("SLURM_NTASKS", 4))
    with Pool(cpuCores) as p:
        radii = list(p.starmap(wrapper_find_intersection_radius, tasks))

    [simplexTree.insert(pair,r) for pair, r in zip(pairs,radii)]

    print('Done.\n')
    return [simplexTree, ellipsoidList]



def adjacent_voronoi_vertices(vertex, vor: Voronoi):
    adjacents = []
    for ridge_idx, ridge in enumerate(vor.ridge_vertices):
        if vertex in ridge:
            adjacents.extend(vor.ridge_points[ridge_idx])
    return list(set(adjacents))



def generate_alpha_ellipsoid_simplex_tree_old(dataset: Dataset,
                                          ellipsoid_parameters: EllipsoidParameters):
    points = dataset.points
    nbhd_size = ellipsoid_parameters.nbhd_size
    axes_ratios = ellipsoid_parameters.axes_ratios
    r_spherisize = ellipsoid_parameters.r_spherisize

    ellipsoid_list: list[Ellipsoid] = fit_ellipsoids(points, nbhd_size, axes_ratios)

    print('Calculating alpha ellipsoid simplex tree... ', end='', flush=True)

    simplex_tree = gd.SimplexTree()
    [simplex_tree.insert([i],0) for i in np.arange(len(points))]

    vor = Voronoi(points)
    for vertex, _ in enumerate(vor.vertices):
        adjacents = adjacent_voronoi_vertices(vertex, vor)

        for i, vertex_i in enumerate(adjacents):
            for _, vertex_j in enumerate(adjacents[i+1:], start=i+1):
                intersection_radius = find_intersection_radius(ellipsoid_list[vertex_i],
                                                               ellipsoid_list[vertex_j],
                                                               r_spherisize=r_spherisize)
                simplex_tree.insert([vertex_i,vertex_j], intersection_radius)

    print('Done.\n')
    return [simplex_tree, ellipsoid_list]


def adjacent_delaunay_vertices(tri, vertex_index):
    """
    Given a Delaunay triangulation and a vertex index, return a list of all adjacent vertices.

    Parameters:
    tri (Delaunay): A Delaunay object containing the triangulation.
    vertex_index (int): The index of the vertex whose adjacent vertices are to be found.

    Returns:
    List[int]: A list of indices of adjacent vertices.
    """
    adjacent_vertices = set()

    # Loop through all simplices (triangles)
    for i, simplex in enumerate(tri.simplices):
        if vertex_index in simplex:
            # Add the other two vertices (those adjacent to the given vertex)
            adjacent_vertices.update(simplex[simplex != vertex_index])

    # Convert set to sorted list for easy viewing
    return sorted(adjacent_vertices)



def generate_alpha_ellipsoid_simplex_tree(dataset: Dataset,
                                          ellipsoid_parameters: EllipsoidParameters):
    points = dataset.points
    nbhd_size = ellipsoid_parameters.nbhd_size
    axes_ratios = ellipsoid_parameters.axes_ratios
    r_spherisize = ellipsoid_parameters.r_spherisize

    ellipsoid_list: list[Ellipsoid] = fit_ellipsoids(points, nbhd_size, axes_ratios)

    print('Calculating alpha ellipsoid simplex tree... ', end='', flush=True)

    simplex_tree = gd.SimplexTree()
    [simplex_tree.insert([i],0) for i in np.arange(len(points))]

    delaunay = Delaunay(points)
    for vertex_index, _ in enumerate(delaunay.points):
        adjacent_vertices = adjacent_delaunay_vertices(delaunay, vertex_index)

        for adjacent_vertex in adjacent_vertices:
            vertex_i = vertex_index
            vertex_j = adjacent_vertex
            intersection_radius = find_intersection_radius(ellipsoid_list[vertex_i],
                                                            ellipsoid_list[vertex_j],
                                                            r_spherisize=r_spherisize)
            simplex_tree.insert([vertex_i,vertex_j], intersection_radius)

    print('Done.\n')
    return [simplex_tree, ellipsoid_list]


def expand_simplex_tree(simplex_tree, expansion_dim=2):
    print("Expanding the simplex tree...", end="", flush=True)
    simplex_tree.expansion(expansion_dim)
    print("Done.")
    return simplex_tree


def collapse_edges(simplex_tree):
    print('Collapsing edges...', end='', flush=True)
    simplex_tree.collapse_edges()
    print('Done.')
    return simplex_tree



def generate_ellipsoid_simplex_tree__single_process(points: np.ndarray, nbhd_size: int, axes_ratios: np.ndarray, r_spherisize: float = np.inf):

    ellipsoidList: list[Ellipsoid] = fit_ellipsoids(points, nbhd_size, axes_ratios)

    print('Calculating ellipsoid simplex tree... ', end='', flush=True)

    simplexTree = gd.SimplexTree()
    [simplexTree.insert([i],0) for i in np.arange(len(points))] #TODO check if the whole array can be assigned
    
    for i, ellipsoid1 in enumerate(ellipsoidList):
        for j, ellipsoid2 in enumerate(ellipsoidList[i:]):
            intersection_radius = find_intersection_radius(ellipsoid1, ellipsoid2)
            simplexTree.insert([i, j+i], intersection_radius)

    print('Done.')
    return [simplexTree, ellipsoidList]



def generate_rips_simplex_tree(points, expansion_dim=2):

    print('Creating the Rips complex... ', end='', flush=True)
    rips_complex = gd.RipsComplex(points=points)
    print('Done.')

    print('Creating the Rips simplex tree... ', end='', flush=True)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    print('Done.')


    return simplex_tree


def double_filtrations(simplex_tree):
    for simplex in simplex_tree.get_simplices():

        simplex_set = simplex[0]  # Get the simplex (set of vertices)
        filtration_value = simplex[1]  # Get the current filtration value
        simplex_tree.assign_filtration(simplex_set, 2*filtration_value)

    return simplex_tree



def generate_alpha_simplex_tree(points):

    print('Creating the alpha complex... ', end='', flush=True)
    alpha_complex = gd.AlphaComplex(points=points)
    print('Done.')

    print('Creating the alpha simplex tree... ', end='', flush=True)
    simplex_tree = alpha_complex.create_simplex_tree()
    print('Done.')

    simplex_tree = double_filtrations(simplex_tree)

    return simplex_tree



def calculate_barcode(simplex_tree, expansion_dim=2, collapse_edges=False):
    # simplexTreeExpanded = simplexTree.copy()
    # if collapse_edges:
    #     print('Collapsing edges...', end='', flush=True)
    #     simplexTreeExpanded.collapse_edges()
    #     print('Done.')

    # print('Expanding the simplex tree... ', end='', flush=True)
    # simplexTreeExpanded.expansion(expansion_dim) # expands the simplicial complex to include
    #                                             # dim-dimensional simplices whose 1-skeleton is in simplexTree
    # print('Done.')

    # print('Calculating the barcode of the expanded tree... ', end='', flush=True)
    # barcode = simplexTreeExpanded.persistence()
    # print('Done.\n')
    print('Calculating the barcode of the expanded tree... ', end='', flush=True)
    barcode = simplex_tree.persistence()
    print('Done.\n')

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



def reduce_barcode(barcode, nBarsDim0 = 10, nBarsDim1 = 10, nBarsDim2 = 10):
    # return only the first nBarsDimk bars in each dimension k
    reduced_barcode = []
    max_bar_end = 0
    for bar in barcode:
        if bar[0] == 0 and nBarsDim0 > 0:
            reduced_barcode.append(bar)
            nBarsDim0 = nBarsDim0 - 1
            max_bar_end = set_max_bar_end(bar, max_bar_end)
        elif bar[0] == 1 and nBarsDim1 > 0:
            reduced_barcode.append(bar)
            nBarsDim1 = nBarsDim1 - 1
            max_bar_end = set_max_bar_end(bar, max_bar_end)
        elif bar[0] == 2 and nBarsDim2 > 0:
            reduced_barcode.append(bar)
            nBarsDim2 = nBarsDim2 - 1
            max_bar_end = set_max_bar_end(bar, max_bar_end)
    
    return reduced_barcode, max_bar_end
    


# def calculateBottleeckDistance(barcode1, barcode2, dim):
#     npBarcode1 = np.array()
#     npBarcode2 = np.array()
#     for line in barcode1:
#         npBarcode1[line[0]].append(line[1])
#     for line in barcode2:
#         npBarcode2[line[0]].append(line[1])

#     bottleneckDistance = [gd.bottleneck_distance(i,j) for i,j in zip(npBarcode1, npBarcode2)]
#     return bottleneckDistance



def pad_axes_ratios(axesRatios: np.ndarray, dim: int):
    ''' For high dimensional ellipsoids, it is enough for the user to specify 
    the first few axes. This function will set the remaining axes to 1.'''
    if dim > len(axesRatios):
        return np.pad(axesRatios, (0,dim - len(axesRatios)), constant_values=1) # creates an array of length dim
    else: 
        return axesRatios[0:dim]



def calculate_ellipsoid_barcode(points: np.ndarray,
                                nbhd_size: int,
                                axes_ratios: np.ndarray,
                                expansion_dim: int = 2,
                                collapse_edges: bool = True,
                                r_spherisize: float = np.inf):
    '''
    Calculate ellipsoid barcode from given points

        Arguments:
            points (np.array): array of points
            nbhd_size (int): size of the neighbourhood for PCA when fitting ellipsoids
            axes_ratios: ratios of axes of ellipsoids (from largest to smallest - last one has to be 1).
                It is sufficient to specify only the non-1 entries, as the rest will be padded with 1's
                to ensure the correct dimensionality.
            expansion_dim (int): [gudhi parameter] dimension to which to expand simplex tree.
                If expansion_dim = n, then only n- and lower-dimensional features will be captured.
            collapse_edges (bool): [gudhi parameter] collapse edges that do not affect the result 
                for faster computation

        Returns:
            barcode (): 
            simplex_tree ():
            ellipsoid_list (list[Ellipsoids]): list of ellipsoids
            total_time (int): total execution time in seconds
    '''
    dim = len(points[0])
    axes_ratios = pad_axes_ratios(axes_ratios,dim)

    t0_simplex_tree = time.time()
    [simplex_tree, ellipsoid_list] = generate_ellipsoid_simplex_tree(points, nbhd_size, axes_ratios, r_spherisize)
    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = calculate_barcode(simplex_tree, expansion_dim, collapse_edges=collapse_edges)
    t1_barcode = time.time()

    total_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    return barcode, simplex_tree, ellipsoid_list, total_time



def calculate_rips_barcode(points: np.ndarray, expansion_dim=2, collapse_edges=True):

    t0_simplex_tree = time.time()
    simplex_tree = generate_rips_simplex_tree(points)
    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = calculate_barcode(simplex_tree, expansion_dim, collapse_edges=collapse_edges)
    t1_barcode = time.time()

    total_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    return barcode, simplex_tree, total_time



def calculate_rips(dataset: Dataset, parameters: Parameters) -> Results:

    t0_simplex_tree = time.time()
    if parameters.complex_subtype == ComplexSubtype.RIPS:
        simplex_tree = generate_rips_simplex_tree(dataset.points)
    elif parameters.complex_subtype == ComplexSubtype.ALPHA:
        simplex_tree = generate_alpha_simplex_tree(dataset.points)
    else:
        exit(f"{parameters.complex_subtype} is an invalid complex type.")
    t1_simplex_tree = time.time()

    if parameters.collapse_edges:
        collapse_edges(simplex_tree)

    if parameters.expansion_dim > 1:
        expand_simplex_tree(simplex_tree)

    t0_barcode = time.time()
    barcode = calculate_barcode(simplex_tree,
                                parameters.expansion_dim,
                                collapse_edges=parameters.collapse_edges)
    t1_barcode = time.time()
    execution_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    results = Results()

    if parameters.save_simplex_tree: results.simplex_tree = simplex_tree

    results.barcode = barcode
    results.execution_time = execution_time

    return results



def calculate_ellipsoids(dataset: Dataset, ellipsoid_parameters: EllipsoidParameters) -> EllipsoidResults:

    dim = dataset.ambient_dim()
    axes_ratios = pad_axes_ratios(ellipsoid_parameters.axes_ratios, dim)

    t0_simplex_tree = time.time()

    if (ellipsoid_parameters.complex_subtype == ComplexSubtype.ALPHA):
        simplex_tree, ellipsoid_list = generate_alpha_ellipsoid_simplex_tree(dataset, ellipsoid_parameters)
    else:
        points = dataset.points
        nbhd_size = ellipsoid_parameters.nbhd_size
        r_spherisize = ellipsoid_parameters.r_spherisize
        [simplex_tree, ellipsoid_list] = generate_ellipsoid_simplex_tree(points, nbhd_size, axes_ratios, r_spherisize)

    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = calculate_barcode(simplex_tree,
                                ellipsoid_parameters.expansion_dim,
                                collapse_edges=ellipsoid_parameters.collapse_edges)
    t1_barcode = time.time()

    execution_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    results = EllipsoidResults()


    if ellipsoid_parameters.save_simplex_tree: results.simplex_tree = simplex_tree
    if ellipsoid_parameters.save_ellipsoid_list: results.ellipsoid_list = ellipsoid_list

    results.barcode = barcode
    results.execution_time = execution_time

    return results
