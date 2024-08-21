import json
import os
import time

import numba
import numpy as np
import gudhi as gd
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree
from scipy import spatial
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from multiprocessing import Pool
from multiprocessing import cpu_count



class Ellipsoid:
    def __init__(self, center: np.array, axes: np.array, axesLengths: np.array):
        self.center = center
        self.axes = axes
        self.axesLengths = axesLengths

    def __eq__(self, other):
        if isinstance(other, Ellipsoid):
            return (np.array_equal(self.center, other.center) and
                    np.array_equal(self.axes, other.axes) and
                    np.array_equal(self.axesLengths, other.axesLengths))
        return False

    def toDict(self):
        obj_data = {
                "center": self.center,
                "axes": self.axes.tolist(),
                "axesLengths": self.axesLengths.tolist()
            }
        return obj_data
    
    def fromDict(self, obj_data: dict):
        self.var1 = obj_data['var1']
        self.var2 = obj_data['var2']
        self.np_array = np.array(obj_data['np_array'])

    def toJSON(self):
        return json.dumps(self.toDict(), 
            sort_keys=True, indent=4)


def fitEllipsoid(center: list, nbhd_pts: np.array, axesRatios: np.array) -> Ellipsoid:
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=len(center))
    pca.fit(nbhd_pts)
    axes = pca.components_
    axesLengths = pca.singular_values_

    if axesRatios.all() != 0:
        axesLengths = axesRatios / axesRatios[0] # r determines the long axis (normalising the long axis to 1)
        # axesLengths = axesRatios / axesRatios[-1] # alt: r determines the short axis
    else: 
        exit("Error: axes ratios contain a zero.")
    return Ellipsoid(center, axes, axesLengths)



def fitEllipsoids(points, neighbourhoodSize, axesRatios) -> list[Ellipsoid]:
    print('Creating KD tree... ', end='', flush=True)
    kdTree = spatial.KDTree(points)
    print('Done.')
    print('Fitting ellipsoids... ', end='', flush=True)

    if len(points) < neighbourhoodSize:
        print('WARNING: the chosen neighbourhood size is too small. \
              Setting the neighbhourhood size to the total number of points.')
        neighbourhoodSize = len(points)

    _,neighbourhoodIdx = kdTree.query(points, neighbourhoodSize)
    neighbourhoods = points[neighbourhoodIdx]
    ellipsoidList \
        = [fitEllipsoid(point, neighbourhood, axesRatios) for point,neighbourhood in zip(points, neighbourhoods)]
    print('Done.')

    return ellipsoidList



def ellipsoidIntersection(ellipsoid1: Ellipsoid, ellipsoid2: Ellipsoid, r):
    ''' Checks whether ellipsoid1 and ellipsoid2 at the filtration level r intersect.
    The method is from https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect
    i.e. this paper: https://tisl.cs.toronto.edu/publication/201207-fusion-kalman_filter_fault_detection/fusion12-kalman_filter_fault_detection.pdf
    :return: true or false
    '''

    if ellipsoid1 == ellipsoid2:
        return True

    Sigma_A = np.linalg.multi_dot([\
        np.transpose(ellipsoid1.axes),\
        np.diag(ellipsoid1.axesLengths**2),\
        ellipsoid1.axes])
    Sigma_B = np.linalg.multi_dot([\
        np.transpose(ellipsoid2.axes),\
        np.diag(ellipsoid2.axesLengths**2),\
        ellipsoid2.axes])
    mu_A = ellipsoid1.center
    mu_B = ellipsoid2.center

    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    res = minimize_scalar(K,
                          bracket=[0.0, 0.5, 1.0],
                          args=(lambdas, v_squared, r))
    return (res.fun >= 0)



def K(s, lambdas, v_squared, r):
    ''' Auxiliary function needed in ellipsoidIntersection
    '''
    return 1.-(1./r**2)*np.sum(v_squared*((s*(1.-s))/(1.+s*(lambdas-1.))))



def get_max_axes_ratio(ellipsoid: Ellipsoid):
    max_axis_length = max(ellipsoid.axesLengths)
    min_axis_length = min(ellipsoid.axesLengths)

    return max_axis_length / min_axis_length



def findIntersectionRadius(ellipsoid1: Ellipsoid, ellipsoid2: Ellipsoid, threshold=0.001, epsilon=0.001):
# def findIntersectionRadius(ellipsoid1: Ellipsoid, ellipsoid2: Ellipsoid, threshold=0.001, epsilon=0.001, *kwargs):

    dist = np.linalg.norm(ellipsoid1.center - ellipsoid2.center)
    
    max_axes_ratio = max(get_max_axes_ratio(ellipsoid1), 
                         get_max_axes_ratio(ellipsoid2))

    maxNonIntersectionFiltration = (dist / 2) - epsilon
    minIntersectionFiltration = dist/2 * max_axes_ratio + epsilon

    r = (minIntersectionFiltration - maxNonIntersectionFiltration)/2

    while True:
        if ellipsoidIntersection(ellipsoid1, ellipsoid2, r):
            minIntersectionFiltration = r
        else: maxNonIntersectionFiltration = r

        if (minIntersectionFiltration - maxNonIntersectionFiltration) < threshold:
            return 2*r # 2r so that it's comparable to Rips
        else: r = (minIntersectionFiltration + maxNonIntersectionFiltration)/2


    
def generateEllipsoidSimplexTree3(points, nbhdSize, axesRatios):
    ''' list comprehesion '''

    ellipsoidList = fitEllipsoids(points, nbhdSize, axesRatios)

    print('Calculating ellipsoid simplex tree... ', end='', flush=True)

    simplexTree = gd.SimplexTree()
    [simplexTree.insert([i],0) for i in np.arange(len(points))]
    
    pairs = [(i,j) for i in np.arange(len(points)) for j in np.arange(i+1,len(points))]
    [simplexTree.insert([i,j], findIntersectionRadius(ellipsoidList[i],ellipsoidList[j],axesRatios=axesRatios)) for (i,j) in pairs]

    print('Done.')

    return [simplexTree, ellipsoidList]



def generateEllipsoidSimplexTree4(points, nbhdSize, axesRatios):
    ''' multiprocessing '''
    ''' Creates a simplex tree from the ellipsoids by adding an edge between each two points whose 
    corresponding ellipsoids intersect.
    :kdTree: KD tree of the initial dataset
    :ellipsoidList: list of ellipsoids (output of ??)
    :queryRadius:
    :filtrationValues:

    :return: gudhi.SimplexTree
    '''

    ellipsoidList = fitEllipsoids(points, nbhdSize, axesRatios)

    print('Calculating ellipsoid simplex tree... ', end='', flush=True)

    simplexTree = gd.SimplexTree()
    [simplexTree.insert([i],0) for i in np.arange(len(points))]
    
    pairs = np.array([[i,j] for i in np.arange(len(points)) for j in np.arange(i+1,len(points))])

    tasks = zip([ellipsoidList[i] for i in pairs[:,0]], \
                [ellipsoidList[j] for j in pairs[:,1]])

    cpuCores = int(os.environ.get("SLURM_NTASKS", 4))

    with Pool(cpuCores) as p:
        radii = list(p.starmap(findIntersectionRadius, tasks))

    [simplexTree.insert(pair,r) for pair, r in zip(pairs,radii)]

    print('Done.')
    return [simplexTree, ellipsoidList]




def generateEllipsoidSimplexTree4_new(points, nbhdSize, axesRatios):
    ''' multiprocessing '''
    ''' Creates a simplex tree from the ellipsoids by adding an edge between each two points whose 
    corresponding ellipsoids intersect.
    :kdTree: KD tree of the initial dataset
    :ellipsoidList: list of ellipsoids (output of ??)
    :queryRadius:
    :filtrationValues:

    :return: gudhi.SimplexTree
    '''

    ellipsoidList: list[Ellipsoid] = fitEllipsoids(points, nbhdSize, axesRatios)

    print('Calculating ellipsoid simplex tree... ', end='', flush=True)

    simplexTree = gd.SimplexTree()
    [simplexTree.insert([i],0) for i in np.arange(len(points))] #TODO check if the whole array can be assigned
    
    for i, ellipsoid1 in enumerate(ellipsoidList):
        for j, ellipsoid2 in enumerate(ellipsoidList[i:]):
            intersection_radius = findIntersectionRadius(ellipsoid1, ellipsoid2)
            simplexTree.insert([i, j+i], intersection_radius)


    print('Done.')
    return [simplexTree, ellipsoidList]


def generateRipsSimplexTree(points):

    print('Creating the Rips complex... ', end='', flush=True)
    ripsComplex = gd.RipsComplex(points=points)
    print('Done.')

    print('Creating the Rips simplex tree... ', end='', flush=True)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=1)
    print('Done.')

    return simplexTreeRips


def expandTreeAndCalculateBarcode(simplexTree, expansionDim, collapseEdges=False):
    simplexTreeExpanded = simplexTree.copy()
    if collapseEdges:
        print('Collapsing edges...', end='', flush=True)
        simplexTreeExpanded.collapse_edges()
        print('Done.')

    print('Expanding the simplex tree... ', end='', flush=True)
    simplexTreeExpanded.expansion(expansionDim) # expands the simplicial complex to include 
                                                # dim-dimensional simplices whose 1-skeleton is in simplexTree
    print('Done.')

    print('Calculating the barcode of the expanded tree... ', end='', flush=True)
    barcode = simplexTreeExpanded.persistence()
    print('Done.')

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

def reduceBarcode(barcode, nBarsDim0 = 1, nBarsDim1 = 0, nBarsDim2 = 0):
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
    
# TODO those blocks below should be separate functions
# def reduceBarcode(barcode, nBarsDim0 = 1, nBarsDim1 = 0, nBarsDim2 = 0):
#     # return only nBarsDimk longest bars in each dimension k
#     reducedBarcode = []
#     maxBarEnd = 0
#     for bar in barcode:
#         if bar[0] == 0 and nBarsDim0 > 0:
#             reducedBarcode.append(bar)
#             nBarsDim0 = nBarsDim0 - 1
#             barEnd = bar[1][1]
#             if barEnd != float('inf') and barEnd > maxBarEnd:
#                 maxBarEnd = barEnd
#         elif bar[0] == 1 and nBarsDim1 > 0:
#             reducedBarcode.append(bar)
#             nBarsDim1 = nBarsDim1 - 1
#             barEnd = bar[1][1]
#             if barEnd != float('inf') and barEnd > maxBarEnd:
#                 maxBarEnd = barEnd
#         elif bar[0] == 2 and nBarsDim2 > 0:
#             reducedBarcode.append(bar)
#             nBarsDim2 = nBarsDim2 - 1
#             barEnd = bar[1][1]
#             if barEnd != float('inf') and barEnd > maxBarEnd:
#                 maxBarEnd = barEnd
    
#     return reducedBarcode, maxBarEnd

def calculateBottleeckDistance(barcode1, barcode2, dim):
    npBarcode1 = np.array()
    npBarcode2 = np.array()
    for line in barcode1:
        npBarcode1[line[0]].append(line[1])
    for line in barcode2:
        npBarcode2[line[0]].append(line[1])

    bottleneckDistance = [gd.bottleneck_distance(i,j) for i,j in zip(npBarcode1, npBarcode2)]
    return bottleneckDistance


def padAxesRatios(axesRatios: np.array, dim: int):
    ''' For high dimensional ellipsoids, it is enough for the user to specify 
    the first few axes. This function will set the remaining axes to 1.'''
    if dim > len(axesRatios):
        return np.pad(axesRatios, (0,dim - len(axesRatios)), constant_values=1) # creates an array of length dim
    else: 
        return axesRatios[0:dim]

    
def calculate_ellipsoid_barcode(points: np.array, nbhd_size: int, axes_ratios: np.array, expansion_dim=2, collapse_edges=True):
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
    axes_ratios = padAxesRatios(axes_ratios,dim)

    t0_simplex_tree = time.time()
    [simplex_tree, ellipsoid_list] = generateEllipsoidSimplexTree4(points, nbhd_size, axes_ratios)
    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = expandTreeAndCalculateBarcode(simplex_tree, expansion_dim, collapseEdges=collapse_edges)
    t1_barcode = time.time()

    total_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    return barcode, simplex_tree, ellipsoid_list, total_time

def calculate_rips_barcode(points: np.array, expansion_dim=2, collapse_edges=True):

    t0_simplex_tree = time.time()
    simplex_tree = generateRipsSimplexTree(points)
    t1_simplex_tree = time.time()

    t0_barcode = time.time()
    barcode = expandTreeAndCalculateBarcode(simplex_tree, expansion_dim, collapseEdges=collapse_edges)
    t1_barcode = time.time()

    total_time = t1_barcode - t0_barcode + t1_simplex_tree - t0_simplex_tree

    return barcode, simplex_tree, total_time
    

    
