import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy import spatial
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

import shapes
import figure_eight
import barcodePlotting

class Ellipsoid:
    def __init__(self, center, axes, axesLengths):
        self.center = center
        self.axes = axes
        self.axesLengths = axesLengths

def createData(nPoints, type, variation = 0.1, dim=2):
    ''' Generate point cloud data of type 'type' and consiting of 'nPoints' points.
    :nPoints: Number of points to generate
    :type: Type of data (e.g. 'circle', 'ellipse', 'Cassini_oval')
    :variation: Amount of variation from the :type:
    :dim: Not used yet
    :return: nPoints x dim numpy array
    '''
    np.random.seed(0)

    if type == 'circle':
        r = 1
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        x = r*np.cos(t) + variation * np.random.rand(nPoints)
        y = r*np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

    elif type == 'ellipse':
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        a = 2
        b = 1
        x = a * np.cos(t) + variation * np.random.rand(nPoints)
        y = b * np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

    elif type == 'Cassini_oval':
        t = np.linspace(-1, 1, int(nPoints/2))
        #t = np.sign(t)*np.abs(t)**(1/4)
        x = np.concatenate((t,t)) + variation * np.random.rand(nPoints)
        yh = (t**2 + 0.5) * np.sqrt(1 - t**2)
        y = np.concatenate((-yh, yh)) + variation * np.random.rand(nPoints)
        
        output = np.vstack((x,y)).transpose()

    else:
        raise Exception(type + " is an invalid type of data.")
    
    return output

def fitEllipsoid(dim, center, neighbourhood):
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=dim)
    pca.fit(neighbourhood)
    axes = pca.components_
    axesLengths = pca.singular_values_
    return Ellipsoid(center, axes, axesLengths)

def plotEllipse(ellipse: Ellipsoid, color='grey', r=1, axes=None):
    sampleRate = 100
    t = np.linspace(0, 2*np.pi, sampleRate)
    xTemp = r*ellipse.axesLengths[0]*np.cos(t)
    yTemp = r*ellipse.axesLengths[1]*np.sin(t)
    x = ellipse.center[0] + ellipse.axes[0,0]*xTemp + ellipse.axes[1,0]*yTemp
    y = ellipse.center[1] + ellipse.axes[0,1]*xTemp + ellipse.axes[1,1]*yTemp
    if axes is None:
        plt.plot(x,y,c=color)
    else:
        axes.plot(x,y,c=color)

def ellipsoidIntersection(ellipsoid1: Ellipsoid, ellipsoid2: Ellipsoid, r):
    ''' Checks whether ellipsoid1 and ellipsoid2 at the filtration level r intersect.
    The method is from https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect
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

def PCAtesting(neighbourhood, ellipse: Ellipsoid):
    # marking the neighbourhood in red:
    plt.scatter(neighbourhood[:,0], neighbourhood[:,1], c='red')

    # plotting the axes of the ellipse:
    plt.quiver(ellipse.center[0],ellipse.center[1], \
               ellipse.axes[0,0], ellipse.axes[0,1], \
               scale=5/(ellipse.axesLengths[0]))
    plt.quiver(ellipse.center[0],ellipse.center[1], \
               ellipse.axes[1,0], ellipse.axes[1,1], \
               scale=5/(ellipse.axesLengths[1]))

def generateEllipoidSimplexTree(kdTree, ellipsoidList, queryRadius, **kwargs):
    ''' Creates a Simplex Tree from the ellipsoids by adding an edge between each two points whose 
    corresponding ellipsoids intersect.
    As **kwargs, one must provide either:
    1)
    :filtrationStart: 
    :filtrationEnd:
    :filtrationStep:
    or
    2)
    :filtration:

    :return: gudhi.SimplexTree
    '''
    points = kdTree.data
    nPoints = len(points)
    simplexTree = gd.SimplexTree()

    if 'filtrationStart' in kwargs and 'filtrationEnd' in kwargs and 'filtrationStep' in kwargs: 
        filtrationStart = kwargs['filtrationStart']
        filtrationEnd = kwargs['filtrationEnd']
        filtrationStep = kwargs['filtrationStep']

        for r in np.arange(filtrationStart, filtrationEnd, filtrationStep):
            for i in range(nPoints):
                simplexTree.insert([i], 0)
                neighboursIdx = kdTree.query_ball_point(points[i], queryRadius, return_sorted=False)
                if len(neighboursIdx) > 1:
                    for idx in neighboursIdx:
                        if (idx != i) and not(simplexTree.find([i,idx])):
                            if ellipsoidIntersection(ellipsoidList[i], ellipsoidList[idx],r):
                                simplexTree.insert([i,idx], r)

    elif 'filtration' in kwargs:
        r = kwargs['filtration']
        for i in range(nPoints):
            simplexTree.insert([i], 0)
            neighboursIdx = kdTree.query_ball_point(points[i], 2, return_sorted=False)
            if len(neighboursIdx) > 1:
                for idx in neighboursIdx:
                    if (idx != i) and not(simplexTree.find([i,idx])):
                        if ellipsoidIntersection(ellipsoidList[i], ellipsoidList[idx],r):
                            simplexTree.insert([i,idx], r)

    else:
        raise Exception("Invalid input argument.")

    return simplexTree

def plotEllipses(ellipseList, r, axes=None):
    for ellipse in ellipseList:
        plotEllipse(ellipse, r=r, axes=axes)

def plotSimplexTree(points, simplexTree, r, axes):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    for splx in simplexList:
        if splx[1] <= r:
            vertices = splx[0]
            if len(vertices) == 1:
                if axes is None:
                    plt.scatter(points[vertices[0]][0], points[vertices[0]][1], c='k', zorder=100)
                else:
                    axes.scatter(points[vertices[0]][0], points[vertices[0]][1], c='k', zorder=100)
            elif len(splx[0]) == 2:
                if axes is None:
                    plt.plot(points[vertices,0], points[vertices,1], c='r')
                else:
                    axes.plot(points[vertices,0], points[vertices,1], c='r')
            elif len(splx[0]) == 3:
                if axes is None:
                    plt.fill(points[vertices,0], points[vertices,1], c='r', alpha=0.1)
                else:
                    axes.fill(points[vertices,0], points[vertices,1], c='r', alpha=0.1)
                    
def printListOfSimplices(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    for splx in simplexList:
        print(splx)

def plotDataPoints(points, axes=None):
    if axes is None:
        plt.scatter(points[:,0],points[:,1])
    else:
        axes.scatter(points[:,0],points[:,1])

def main():

    dim = 2                 # dimension of the ambient space
    r = 0                   # filtration parameter (if set to zero, 
                            # calculations will be performed for a range of filtrations)
    rStart = 0.1
    rEnd = 4
    rStep = 0.05
    rPlot = 0.7            # the value of r at which the simplex tree will be plotted
    if rPlot not in np.arange(rStart,rEnd,rStep):
        print('Warning: the simplex tree plot may be inaccurate since the calculations are ' \
              +'not performed for the chosen value of rPlot. To fix this, make sure that ' \
              +'rPlot is in np.arange(rStart,rEnd,rStep).')
    neighbourhoodSize = 3   # number of points for doing PCA

    nPoints = 10            # number of data points
    points = createData(nPoints,'circle')
    #points = shapes.sample_from_sphere(n=nPoints)
    #points = figure_eight.figure_eight(nPoints, 1, 0.2)

    kdTree = spatial.KDTree(points)

    ellipseList = []
    for point in points:
        [NULL,neighbourhoodIdx] = kdTree.query(point, min(nPoints,neighbourhoodSize))
        neighbourhood = points[neighbourhoodIdx]
        currentEllipsoid = fitEllipsoid(dim, point, neighbourhood)
        ellipseList.append(currentEllipsoid)

    longestEllipsoidAxis = max(ellipsoid.axesLengths[0] for ellipsoid in ellipseList)
    queryRadius = 2*longestEllipsoidAxis
        
    if r > 0: 
        simplexTree = generateEllipoidSimplexTree(kdTree, ellipseList, queryRadius, \
                                                  filtration=r)
    else:
        simplexTree = generateEllipoidSimplexTree(kdTree, ellipseList, queryRadius, \
                                                  filtrationStart=rStart, filtrationEnd=rEnd, \
                                                  filtrationStep=rStep)
    
    simplexTree.expansion(dim) # expands the simplicial complex to include 
                               # dim-dimensional simplices whose 1-skeleton is in simplexTree
    barcode = simplexTree.persistence()

    ripsComplex = gd.RipsComplex(points=points)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=dim)
    barcodeRips = simplexTreeRips.persistence()

    fig = plt.figure()
    gs = fig.add_gridspec(2,2)
    axData = fig.add_subplot(gs[:, 0])
    axBarE = fig.add_subplot(gs[0, 1])
    axBarR = fig.add_subplot(gs[1, 1])

    plotEllipses(ellipseList, rPlot, axes=axData)
    plotSimplexTree(points, simplexTree, rPlot, axes=axData)
    axData.set_title('Data and the ellipsoid simplex tree for r = %0.2f' %(rPlot), fontsize=12)

    barcodePlotting.plot_persistence_barcode(barcode, inf_delta=0.5, axes=axBarE, fontsize=12,\
                                             axis_start = -0.1, infinity = rEnd + 0.1)
    axBarE.set_title('Ellipsoid barcode', fontsize=12)
    barcodePlotting.plot_persistence_barcode(barcodeRips, inf_delta=0.5, axes=axBarR, fontsize=12,\
                                             axis_start = -0.1, infinity = rEnd + 0.1)
    axBarR.set_title('Rips barcode', fontsize=12)
    plt.show()

if __name__ == "__main__":
    main()

    


