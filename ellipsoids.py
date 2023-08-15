import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree
from scipy import spatial
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

import shapes
import figure_eight
import barcodePlotting

import json
from datetime import datetime
import time

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

    elif type == 'sphere':
        r = 1
        

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

def fitEllipsoids(dim, kdTree, neighbourhoodSize):
    points = kdTree.data
    nPoints = len(points)
    ellipseList = []
    for point in points:
        [NULL,neighbourhoodIdx] = kdTree.query(point, min(nPoints,neighbourhoodSize))
        neighbourhood = points[neighbourhoodIdx]
        currentEllipsoid = fitEllipsoid(dim, point, neighbourhood)
        ellipseList.append(currentEllipsoid)
    return ellipseList

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

def plotEllipsoid(ellipsoid: Ellipsoid, color='grey', r=1, axes=None):
    # NOT FINSIHED
    sampleRate = 100

    rx = ellipsoid.axesLengths[0]
    ry = ellipsoid.axesLengths[1]
    rz = ellipsoid.axesLengths[2]
    
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, sampleRate)
    v = np.linspace(0, np.pi, sampleRate)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    a1 = rx * np.outer(np.cos(u), np.sin(v))
    a2 = ry * np.outer(np.sin(u), np.sin(v))
    a3 = rz * np.outer(np.ones_like(u), np.cos(v))

    transformationMatrix = ellipsoid.axes()
    x = a1*transformationMatrix
    y = a2*transformationMatrix
    z = a3*transformationMatrix

    if axes is None:
        plt.plot_surface(x,y,z, rstride=4, cstride=4, c=color)
    else:
        axes.plot_surface(x,y, rstride=4, cstride=4, c=color)

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

def generateEllipoidSimplexTree(kdTree, ellipsoidList, queryRadius, filtrationValues):
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
                    #debug:print(f'{i=}' + ' ' f'{idx}')
                    if (idx != i) and not(simplexTree.find([i,idx])):
                        #debug:print(f'{ellipsoidList[i].axes=}' + '; ' + f'{ellipsoidList[idx].axes=}')
                        if ellipsoidIntersection(ellipsoidList[i], ellipsoidList[idx],r):
                            simplexTree.insert([i,idx], r)

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

def maxFiltration(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    return max(splx[1] for splx in simplexList)

def plotDataPoints(points, axes=None):
    if axes is None:
        plt.scatter(points[:,0],points[:,1])
    else:
        axes.scatter(points[:,0],points[:,1])

class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        print(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Ellipsoid):
            obj_data = {
                "center": obj.center,
                "axes": obj.axes.tolist(),
                "axesLengths": obj.axesLengths.tolist()
            }
            return obj_data
        elif isinstance(obj,gd.SimplexTree):
            return list(obj.get_filtration())
        return json.JSONEncoder.default(self, obj)
    
def saveVarsToFile(dim, rStart, rEnd, rStep, rValues, nbhdSize, nPts, points, ellipseList, simplexTreeEllipsoids, \
                simplexTreeRips, barcodeEllipsoids, barcodeRips, \
                filename=datetime.now().strftime("data/test.json")):
    print('Saving data to file...')
    dictOfVars = {
        'dim': dim,
        'rStart': rStart,
        'rEnd': rEnd,
        'rStep': rStep,
        'rValues': rValues,
        'nbhdSize': nbhdSize,
        'nPts': nPts,
        'points': points,
        'ellipseList': ellipseList,
        'simplexTreeEllipsoids': simplexTreeEllipsoids,
        'simplexTreeRips': simplexTreeRips,
        'barcodeEllipsoids': barcodeEllipsoids,
        'barcodeRips': barcodeRips
    }
    json_string = json.dumps(dictOfVars, cls=CustomEncoder, indent=4)
    with open(filename, 'w') as outfile:
        outfile.write(json_string)
    print("Data saved to file.")

def loadVarsFromFile(filename):
    with open(filename, "r") as f:
        jsonVars = json.load(f)
    points = np.asarray(jsonVars['points'])

    ellipseListRaw = jsonVars['ellipseList']
    ellipseList = []
    for ellipse in ellipseListRaw:
        ellipseList.append(Ellipsoid(ellipse['center'], np.asarray(ellipse['axes']), \
                                     np.asarray(ellipse['axesLengths'])))

    simplexTreeEllipsoidsRaw = jsonVars['simplexTreeEllipsoids']
    simplexTreeEllipsoids = gd.SimplexTree()
    for simplexTreeEntry in simplexTreeEllipsoidsRaw:
        simplexTreeEllipsoids.insert(simplexTreeEntry[0],simplexTreeEntry[1])
    
    simplexTreeRipsRaw = jsonVars['simplexTreeEllipsoids']
    simplexTreeRips = gd.SimplexTree()
    for simplexTreeEntry in simplexTreeRipsRaw:
        simplexTreeRips.insert(simplexTreeEntry[0],simplexTreeEntry[1])

    barcodeEllipsoids = jsonVars['barcodeEllipsoids']
    barcodeRips = jsonVars['barcodeRips']

    return [points, ellipseList, simplexTreeEllipsoids, \
            simplexTreeRips, barcodeEllipsoids, barcodeRips]

def visualisation(**kwargs):
    print('Generating plots...')

    points = kwargs['points']
    simplexTreeEllipsoids = kwargs['simplexTreeEllipsoids']
    simplexTreeRips = kwargs['simplexTreeRips']
    barcodeEllipsoids = kwargs['barcodeEllipsoids']
    barcodeRips = kwargs['barcodeRips']
    dim = len(points[0,:])

    if dim == 2 or dim == 3:
        fig = plt.figure(figsize=(14,7))
        gs = fig.add_gridspec(2,2)
        if dim == 2:
            axData = fig.add_subplot(gs[:, 0])
        else:
            axData = fig.add_subplot(gs[:,0], projection='3d')
        axBarE = fig.add_subplot(gs[0, 1])
        axBarR = fig.add_subplot(gs[1, 1])

        for point in points:
            axData.scatter(*point, c='k')
        axData.set_title('Data (%d points)' %len(points), fontsize=12)

    else:
        fig = plt.figure(figsize=(14,7))
        gs = fig.add_gridspec(2,1)
        axBarE = fig.add_subplot(gs[0, 1])
        axBarR = fig.add_subplot(gs[1, 1])

    if 'filename' in kwargs:
        filename = kwargs['filename']
    else: filename = 'data/plotTest.png'

    if ('ellipseList' in kwargs) and ('rPlot' in kwargs) and (kwargs['rPlot'] != 0):
        ellipseList = kwargs['ellipseList']
        rPlot = kwargs['rPlot']
        plotEllipses(ellipseList, rPlot, axes=axData)
        plotSimplexTree(points, simplexTreeEllipsoids, rPlot, axes=axData)
        axData.set_title('Data and the ellipsoid simplex tree for r = %0.2f' %(rPlot), fontsize=12)

    xAxisEnd = max(maxFiltration(simplexTreeEllipsoids), maxFiltration(simplexTreeRips)) + 0.1
    barcodePlotting.plot_persistence_barcode(barcodeEllipsoids, inf_delta=0.5, axes=axBarE, fontsize=12,\
                                             axis_start = -0.1, infinity = xAxisEnd)
    axBarE.set_title('Ellipsoid barcode', fontsize=12)
    barcodePlotting.plot_persistence_barcode(barcodeRips, inf_delta=0.5, axes=axBarR, fontsize=12,\
                                             axis_start = -0.1, infinity = xAxisEnd + 0.1)
    axBarR.set_title('Rips barcode', fontsize=12)

    if 'rValues' in kwargs:
        rValues = kwargs['rValues']
        for rValue in rValues:
            axBarE.axvline(x = rValue, color='gray', linewidth=0.5, linestyle='dashed')
    
    if 'showPlot' in kwargs and kwargs['showPlot'] is True:
        plt.show()
    
    if 'savePlot' in kwargs and kwargs['savePlot'] is True:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print('Plot saved to file.')

def visualisationFromFile(filename, rPlot=0.6):
    [points, ellipseList, simplexTreeEllipsoids, simplexTreeRips, \
     barcodeEllipsoids, barcodeRips] = loadVarsFromFile(filename)
    visualisation(points=points, \
                  simplexTreeEllipsoids=simplexTreeEllipsoids, \
                  simplexTreeRips=simplexTreeRips, \
                  barcodeEllipsoids=barcodeEllipsoids, \
                  barcodeRips=barcodeRips, \
                  showPlot = True,
                  rPlot = rPlot,
                  ellipseList = ellipseList
                  )

def readOFF(filename):
    file = open(filename, 'r')
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    nVerts, nFaces, nEdges = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')] for i_vert in range(nVerts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for i_face in range(nFaces)]
    return np.asarray(verts)

def COOsignature(persistanceBarcode):
    diagramPts = []
    for entry in persistanceBarcode:
        diagramPts.append(entry[1])

    diagramKDTree = KDTree(diagramPts, metric='chebyshev')   
    for diagramPt in diagramPts:
        distClosest, NULL = diagramKDTree.query(diagramPt,k=1)

    pass

def main():
    
    ###### User input ######
    boolSaveData = False
    boolShowPlot = True
    boolSavePlot = False
    rPlot = 0
    # -------------------- #
    dim = 2              # dimension of the ambient space
    rStart = 0.2
    rEnd = 6
    rStep = 0.5
    rValues = np.arange(rStart, rEnd, rStep)
    #rValues = [0.5, 1, 2, 5]
    nbhdSize = 3         # number of points for doing PCA
    nPts = 20            # number of data points
    #points = createData(nPts,'circle')
    #points = shapes.sample_from_sphere(n=nPts, d=(dim-1), seed=0)
    points = figure_eight.figure_eight(nPts, 1, 0.2)
    #points = readOFF('data/61.off')
    #points = points[1:100,:]
    ########################

    rStart = rValues[0]
    rEnd = rValues[-1]
    rStep = abs(rValues[1] - rValues[0])
    nPts = len(points)
    
    importantParameters = f'{nPts=}' + '_'+ f'{rStep=}' + '_' + f'{nbhdSize=}'
    dataPlotFilename='data/ellipsoids_'+importantParameters+datetime.now().strftime("_%Y%m%d_%H%M%S")
    
    if (boolShowPlot or boolSavePlot) and (rPlot not in rValues):
        print('\nWarning: the simplex tree plot may be inaccurate since the calculations are ' \
            +'not performed for the chosen value of rPlot. To fix this, make sure that ' \
            +'rPlot is in np.arange(rStart,rEnd,rStep).')

    tStart = time.time()    # for testing performace

    kdTree = spatial.KDTree(points)
    ellipseList = fitEllipsoids(dim, kdTree, nbhdSize)
    #debug:for ellipse in ellipseList:
    #debug:    print(f'{ellipse.axes=}')
    longestEllipsoidAxis = max(ellipsoid.axesLengths[0] for ellipsoid in ellipseList)
    queryRadius = 2*longestEllipsoidAxis
    try:
        simplexTreeEllipsoids = generateEllipoidSimplexTree(kdTree, ellipseList, queryRadius, \
                                                        filtrationValues = rValues)
    except np.linalg.LinAlgError:
        print("\nError: Attempting to add an edge in Ellipsoid Simplex from a vertex to itself. \n" + 
                "Reason: Two of the ellipsoids are the same because two points from " + \
                "the data set have the same neighbourhood. Please change the " + \
                "neighbourhood size (nbhdSize) and try again.")
        return 0
    simplexTreeEllipsoids.expansion(dim) # expands the simplicial complex to include 
                                        # dim-dimensional simplices whose 1-skeleton is in simplexTree
    barcodeEllipsoids = simplexTreeEllipsoids.persistence()

    ripsComplex = gd.RipsComplex(points=points)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=dim)
    barcodeRips = simplexTreeRips.persistence()

    tEnd = time.time()      # for testing performace
    print('The total execution time is ' + str(tEnd-tStart))
    
    if boolSaveData is True:
        saveVarsToFile(dim, rStart, rEnd, rStep, rValues, nbhdSize, nPts, \
                        points, ellipseList, simplexTreeEllipsoids, \
                        simplexTreeRips, barcodeEllipsoids, barcodeRips, \
                        filename = dataPlotFilename + '.json')

    if (boolShowPlot or boolSavePlot) is True:
        visualisation(points = points,\
                        ellipseList = ellipseList, rPlot = rPlot, \
                        simplexTreeEllipsoids = simplexTreeEllipsoids, \
                        simplexTreeRips = simplexTreeRips, \
                        barcodeEllipsoids = barcodeEllipsoids, \
                        barcodeRips = barcodeRips, \
                        showPlot = boolShowPlot, 
                        savePlot = boolSavePlot,
                        filename = dataPlotFilename + '.png', \
                        rValues = rValues)
        

if __name__ == "__main__":
    main()




    


