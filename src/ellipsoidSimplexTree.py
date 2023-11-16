import numpy as np
import gudhi as gd
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree
from scipy import spatial
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import json

from multiprocessing import Pool
from multiprocessing import cpu_count
import itertools

class Ellipsoid:
    def __init__(self, center, axes, axesLengths):
        self.center = center
        self.axes = axes
        self.axesLengths = axesLengths

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

def fitEllipsoid(center, neighbourhood, axesRatios=0):
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=len(center))
    pca.fit(neighbourhood)
    axes = pca.components_
    axesLengths = pca.singular_values_
    if axesRatios.all() != 0:
        axesLengths = axesRatios / axesRatios[0]
        # alt20230927: if r should determine the short axis
        # axesLengths = axesRatios / axesRatios[-1] 
        # /alt
    return Ellipsoid(center, axes, axesLengths)

def fitEllipsoids(points, neighbourhoodSize, axesRatios = 0):
    print('Creating KD tree... ', end='', flush=True)
    kdTree = spatial.KDTree(points)
    print('Done.')
    print('Fitting ellipsoids... ', end='', flush=True)

    if len(points) < neighbourhoodSize:
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

def generateEllipsoidSimplexTree(kdTree, ellipsoidList, queryRadius, filtrationValues):
    ''' Creates a simplex tree from the ellipsoids by adding an edge between each two points whose 
    corresponding ellipsoids intersect.
    :kdTree: KD tree of the initial dataset
    :ellipsoidList: list of ellipsoids (output of ??)
    :queryRadius:
    :filtrationValues:

    :return: gudhi.SimplexTree
    '''
    points = kdTree.data
    nPoints = len(points)
    simplexTree = gd.SimplexTree()

    for r in filtrationValues:
        for i in range(nPoints):
            simplexTree.insert([i], 0)
            neighboursIdx = kdTree.query_ball_point(points[i], queryRadius, return_sorted=False)
            if len(neighboursIdx) > 1:
                for idx in neighboursIdx:
                    if (idx != i) and not(simplexTree.find([i,idx])):
                        if ellipsoidIntersection(ellipsoidList[i], ellipsoidList[idx],r):
                            simplexTree.insert([i,idx], r)

    return simplexTree

def generateEllipsoidSimplexTree2(points, nbhdSize, axesRatios):
    ''' Creates a simplex tree from the ellipsoids by adding an edge between each two points whose 
    corresponding ellipsoids intersect.
    :kdTree: KD tree of the initial dataset
    :ellipsoidList: list of ellipsoids (output of ??)
    :queryRadius:
    :filtrationValues:

    :return: gudhi.SimplexTree
    '''
    
    ellipsoidList = fitEllipsoids(points, nbhdSize, axesRatios)
    longestEllipsoidAxis = max(ellipsoid.axesLengths[0] for ellipsoid in ellipsoidList)

    simplexTree = gd.SimplexTree()
    distanceMatrix = spatial.distance.squareform(spatial.distance.pdist(points))
    threshold = 0.001
    epsilon = 0.001

    print('Calculating ellipsoid simplex tree... ', end='', flush=True)

    nPoints = len(points)
    for i in range(nPoints):
        simplexTree.insert([i], 0)
        for j in range(i+1,nPoints): # j goes from i+1 so as to check for intersections only once per pair of ellipsoids
            dist = distanceMatrix[i,j]
            if axesRatios.all() != 0:
                maxNonIntersectionFiltration = (dist / 2) - epsilon
                minIntersectionFiltration = dist/2 * max(axesRatios) + epsilon
            else:
                maxNonIntersectionFiltration = (dist / 2) - epsilon
                minIntersectionFiltration = dist/2 * longestEllipsoidAxis + epsilon
            # alt20230927: if r should determine the short axis
            # maxNonIntersectionFiltration = (dist / 2) / max(axesRatios) - epsilon
            # minIntersectionFiltration = dist/2 + epsilon
            # /alt

            r = (minIntersectionFiltration - maxNonIntersectionFiltration)/2

            while True:
                if ellipsoidIntersection(ellipsoidList[i], ellipsoidList[j], r):
                    minIntersectionFiltration = r
                else: maxNonIntersectionFiltration = r

                if (minIntersectionFiltration - maxNonIntersectionFiltration) < threshold:
                    # simplexTree.insert([i,j], r)
                    simplexTree.insert([i,j], 2*r) # alt20230927_2: 2r so that it's comparable to Rips
                    break
                else: r = (minIntersectionFiltration + maxNonIntersectionFiltration)/2

    print('Done.')
    return [simplexTree, ellipsoidList]


def findIntersectionRadius(ellipsoid1, ellipsoid2, axesRatios, **kwargs):
    if 'threshold' not in kwargs:
        threshold = 0.001
    if 'epsilon' not in kwargs:
        epsilon = 0.001
    if 'dist' not in kwargs:
        dist = np.linalg.norm(ellipsoid1.center - ellipsoid2.center)
    if 'axesRatios' in kwargs and kwargs['axesRatios'].all() != 0:
        maxNonIntersectionFiltration = (dist / 2) - epsilon
        minIntersectionFiltration = dist/2 * max(kwargs['axesRatios']) + epsilon
    elif 'longestEllipsoidAxes' in kwargs:
        maxNonIntersectionFiltration = (dist / 2) - epsilon
        minIntersectionFiltration = dist/2 * kwargs['longestEllipsoidAxis'] + epsilon
    else:
        maxAxesRatio = max(ellipsoid1.axesLengths / min(ellipsoid1.axesLengths) + \
                           ellipsoid2.axesLengths / min(ellipsoid2.axesLengths))
        maxNonIntersectionFiltration = (dist / 2) - epsilon
        minIntersectionFiltration = dist/2 * maxAxesRatio + epsilon

    # temp
    maxNonIntersectionFiltration = (dist / 2) - epsilon
    minIntersectionFiltration = dist/2 * max(axesRatios) + epsilon

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
                [ellipsoidList[j] for j in pairs[:,1]], \
                itertools.repeat(axesRatios, len(pairs)))

    cpus = cpu_count()

    with Pool(cpus) as p:
        radii = list(p.starmap(findIntersectionRadius, tasks))

    [simplexTree.insert(pair,r) for pair, r in zip(pairs,radii)]

    print('Done.')
    return [simplexTree, ellipsoidList]