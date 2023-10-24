# problems:
# some close ones are touching but no edge (e.g. around 2am)
# ratios not always same 
# - maybe axes are not noramlised? they seem normalised

import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree
from scipy import spatial
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

import scipy.io

import shapes
import figure_eight
import barcodePlotting

import json
from datetime import datetime
import time

from sys import getsizeof

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

def createData(nPoints, type, variation = 0.1, dim=2, outliers=False):
    ''' Generate point cloud data of type 'type' and consiting of 'nPoints' points.
    :nPoints: Number of points to generate
    :type: Type of data (e.g. 'circle', 'ellipse', 'Cassini_oval')
    :variation: Amount of variation from the :type:
    :dim: Not used yet
    :return: nPoints x dim numpy array
    '''
    np.random.seed(0)

    if type == 'circle':
        if outliers is True: nPoints = nPoints - 1
        r = 1
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        x = r*np.cos(t) + variation * np.random.rand(nPoints)
        y = r*np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

        if outliers is True:
            output = np.append(output,[[0,0]],axis=0)

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

    elif type == 'torus':
        r = 1
        R = 2

        nPtsSampled = 0
        theta = np.zeros([nPoints])
        phi = np.zeros([nPoints])

        # rejection sampling
        while nPtsSampled < nPoints:
            thetaSample = 2 * np.pi * np.random.rand()
            phiSample = 2 * np.pi * np.random.rand()
            W = 2 * np.pi * np.random.rand()

            if W <= (R + r * np.cos(thetaSample))/(R+r):
                theta[nPtsSampled] = thetaSample
                phi[nPtsSampled] = phiSample
                nPtsSampled += 1


            x = (R + r * np.cos(theta)) * np.cos(phi)
            y = (R + r * np.cos(theta)) * np.sin(phi)
            z = r * np.sin(theta)

        output = np.vstack((x,y,z)).transpose()
        

    else:
        raise Exception(type + " is an invalid type of data.")
    
    return output

def fitEllipsoid(dim, center, neighbourhood, axesRatios=0):
    ''' Use PCA to fit an ellipsoid to the given neighbourhood
    :return: ellipsoid of dimension dim with axes obtained from PCA
    '''
    pca = PCA(n_components=dim)
    pca.fit(neighbourhood)
    axes = pca.components_
    axesLengths = pca.singular_values_
    if axesRatios.all() != 0:
        axesLengths = axesRatios / axesRatios[0]
        # alt20230927: if r should determine the short axis
        # axesLengths = axesRatios / axesRatios[-1] 
        # /alt
    return Ellipsoid(center, axes, axesLengths)

def fitEllipsoids(dim, kdTree, neighbourhoodSize, axesRatios = 0):
    points = kdTree.data
    nPoints = len(points)
    ellipseList = []
    for point in points:
        [NULL,neighbourhoodIdx] = kdTree.query(point, min(nPoints,neighbourhoodSize))
        neighbourhood = points[neighbourhoodIdx]
        currentEllipsoid = fitEllipsoid(dim, point, neighbourhood, axesRatios)
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
    # see https://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
    sampleRate = 100

    rx = r * ellipsoid.axesLengths[0]
    ry = r * ellipsoid.axesLengths[1]
    rz = r * ellipsoid.axesLengths[2]
    
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, sampleRate)
    v = np.linspace(0, np.pi, sampleRate)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    all = np.concatenate((np.reshape(x, [-1,1]), np.reshape(y, [-1,1]), np.reshape(z, [-1,1])), axis=1)
    allTransformed = ellipsoid.axes @ np.transpose(all)
    
    x = np.reshape(allTransformed[0,:],(100,100)) + ellipsoid.center[0]
    y = np.reshape(allTransformed[1,:],(100,100)) + ellipsoid.center[1]
    z = np.reshape(allTransformed[2,:],(100,100)) + ellipsoid.center[2]

    if axes is None:
        plt.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)
    else:
        axes.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)

def plotEllipses(ellipseList, r, axes=None):
    for ellipse in ellipseList:
        plotEllipse(ellipse, r=r, axes=axes)

def plotEllipsoids(ellipsoidList, r, axes=None):
    for ellipsoid in ellipsoidList:
        plotEllipsoid(ellipsoid, r = r, axes=axes)

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
                    if (idx != i) and not(simplexTree.find([i,idx])):
                        if ellipsoidIntersection(ellipsoidList[i], ellipsoidList[idx],r):
                            simplexTree.insert([i,idx], r)

    return simplexTree

def generateEllipoidSimplexTree2(kdTree, ellipsoidList, axesRatios):
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
    distanceMatrix = spatial.distance.squareform(spatial.distance.pdist(points))
    threshold = 0.001
    epsilon = 0.001

    for i in range(nPoints):
        simplexTree.insert([i], 0)
        for j in range(i+1,nPoints):
            dist = distanceMatrix[i,j]
            if axesRatios.all() != 0:
                maxNonIntersectionFiltration = (dist / 2) - epsilon
                minIntersectionFiltration = dist/2 * max(axesRatios) + epsilon
            else:
                longestEllipsoidAxis = max(ellipsoid.axesLengths[0] for ellipsoid in ellipsoidList)
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

    return simplexTree

def plotSimplexTree(points, simplexTree, r, axes):
    dim = len(points[0])
    if dim > 3:
        raise Exception('Error: Attempting to plot simplex tree in dimension higher than 3.')
    generator = simplexTree.get_filtration()
    simplexList = list(generator)

    if axes is None:
        for splx in simplexList:
            if splx[1] <= r:
                vertices = splx[0]
                match len(vertices):
                    case 1:
                        plt.scatter(*np.transpose(points[vertices]), c='k', zorder=100)
                        # points[vertices] gives us points forming the vertices of the simplex
                        # transposing them and taking the * operator returns x-, y-. and z-coords separately
                    case 2:
                        plt.plot(*np.transpose(points[vertices]), c='r')
                    case 3:
                        if dim == 2:
                            plt.fill(*np.transpose(points[vertices]), c='r', alpha=0.1)
    else:
        idx = 1
        for splx in simplexList:
            idx = idx + 1
            
            # if splx[1] <= r:
            if splx[1] <= 2*r: # alt20230927_2: 2r so that it's comparable to Rips
                vertices = splx[0]
                match len(vertices):
                    case 1:
                        axes.scatter(*np.transpose(points[vertices]), c='k', zorder=100)
                    case 2:
                        axes.plot(*np.transpose(points[vertices]), c='r')
                    case 3:
                        if dim == 2:
                            axes.fill(*np.transpose(points[vertices]), c='r', alpha=0.1)

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
    
def saveVarsToFile(dictOfVars,
                   filename=datetime.now().strftime("data/test.json")):
    print('Saving data to file...')
    json_string = json.dumps(dictOfVars, cls=CustomEncoder, indent=4)
    with open(filename, 'w') as outfile:
        outfile.write(json_string)
    print("Data saved to file " + filename + '.')

def loadVarsFromFile(filename):
    with open(filename, "r") as f:
        jsonVars = json.load(f)
    
    vars = {}

    if 'dim' in jsonVars:
        vars['dim'] = jsonVars['dim']

    if 'rStart' in jsonVars:
        vars['rStart'] = jsonVars['rStart']

    if 'rEnd' in jsonVars:
        vars['rEnd'] = jsonVars['rEnd']

    if 'rStep' in jsonVars:
        vars['rStep'] = jsonVars['rStep']

    if 'rValues' in jsonVars:
        vars['rValues'] = np.asarray(jsonVars['rValues'])

    if 'nbhdSize' in jsonVars:
        vars['nbhdSize'] = jsonVars['rStep']

    if 'nPts' in jsonVars:
        vars['nPts'] = jsonVars['nPts']

    if 'points' in jsonVars:
        vars['points'] = np.asarray(jsonVars['points'])

    if 'ellipseList' in jsonVars:
        ellipsoidListRaw = jsonVars['ellipseList']
        ellipsoidList = []
        for ellipsoid in ellipsoidListRaw:
            ellipsoidList.append(Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                     np.asarray(ellipsoid['axesLengths'])))
        vars['ellipsoidList'] = ellipsoidList

    if 'ellipsoidList' in jsonVars:
        ellipsoidListRaw = jsonVars['ellipsoidList']
        ellipsoidList = []
        for ellipsoid in ellipsoidListRaw:
            ellipsoidList.append(Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                     np.asarray(ellipsoid['axesLengths'])))
        vars['ellipsoidList'] = ellipsoidList

    if 'simplexTreeEllipsoids' in jsonVars:
        simplexTreeEllipsoidsRaw = jsonVars['simplexTreeEllipsoids']
        simplexTreeEllipsoids = gd.SimplexTree()
        for simplexTreeEntry in simplexTreeEllipsoidsRaw:
            simplexTreeEllipsoids.insert(simplexTreeEntry[0],simplexTreeEntry[1])
        vars['simplexTreeEllipsoids'] = simplexTreeEllipsoids

    if 'simplexTreeRips' in jsonVars:        
        simplexTreeRipsRaw = jsonVars['simplexTreeEllipsoids']
        simplexTreeRips = gd.SimplexTree()
        for simplexTreeEntry in simplexTreeRipsRaw:
            simplexTreeRips.insert(simplexTreeEntry[0],simplexTreeEntry[1])
        vars['simplexTreeRips'] = simplexTreeRips

    if 'barcodeEllipsoids' in jsonVars:
        vars['barcodeEllipsoids'] = jsonVars['barcodeEllipsoids']

    if 'barcodeRips' in jsonVars:
        vars['barcodeRips'] = jsonVars['barcodeRips']

    return vars
           
def visualisation(**kwargs):
    print('Generating plots...')

    points = kwargs['points']
    print('size of 1 point = ' + str(getsizeof(points[1])))
    dim = len(points[0,:])
    expansionDim = 3 #
    simplexTreeEllipsoids = kwargs['simplexTreeEllipsoids']
    simplexTreeEllipsoids.expansion(expansionDim)
    # simplexTreeRips = kwargs['simplexTreeRips']
    ripsComplex = gd.RipsComplex(points=points)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=expansionDim)
    barcodeEllipsoids = kwargs['barcodeEllipsoids']
    barcodeRips = kwargs['barcodeRips']

    # plotting ellipsoids (only if 2- or 3-dimensional)
    # if len(points) < 50 and (dim == 2 or dim == 3):
    plotPoints = False
    if plotPoints and (dim == 2 or dim == 3):
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
        if ('ellipsoidList' in kwargs or 'ellipseList' in kwargs) and ('rPlot' in kwargs) and (kwargs['rPlot'] != 0):
            ellipsoidList = kwargs['ellipsoidList'] if 'ellipsoidList' in kwargs else kwargs['ellipseList']
            rPlot = kwargs['rPlot']
            if dim == 2:
                plotEllipses(ellipsoidList, rPlot, axes=axData)
            if dim == 3:
                plotEllipsoids(ellipsoidList, rPlot, axes=axData)
            plotSimplexTree(points, simplexTreeEllipsoids, rPlot, axes=axData)
            axData.set_title('Data and the ellipsoid simplex tree for r = %0.2f' %(rPlot), fontsize=12)
        axData.set_aspect('equal', adjustable='box')

    else:
        fig = plt.figure(figsize=(14,7))
        gs = fig.add_gridspec(2,1)
        axBarE = fig.add_subplot(gs[0, 0])
        axBarR = fig.add_subplot(gs[1, 0])


    if 'filename' in kwargs:
        filename = kwargs['filename']
    else: filename = 'data/plotTest.png'

    # plotting barcodes
    # xAxisEnd = max(maxFiltration(simplexTreeEllipsoids)*2, maxFiltration(simplexTreeRips)) + 0.1
    # barcodePlotting.plot_persistence_barcode(barcodeEllipsoids, inf_delta=0.5, axes=axBarE, fontsize=12,\
    #                                          axis_start = -0.1, infinity = xAxisEnd / 2) #(0.1 + xAxisEnd) /2 ) # todo: put this back to xAxisEnd (without /2)
    # axBarE.set_title('Ellipsoid barcode', fontsize=12)
    # barcodePlotting.plot_persistence_barcode(barcodeRips, inf_delta=0.5, axes=axBarR, fontsize=12,\
    #                                          axis_start = -0.1, infinity = xAxisEnd) #(0.1 + xAxisEnd))
    # axBarR.set_title('Rips barcode', fontsize=12)

    # alt20230927_2: ellipsoid simplex filtration is 2r so that it's comparable to Rips
    # plotting barcodes
    print('Determining the xAxisLength... ', end='', flush=True)
    # xAxisEnd = max(maxFiltration(simplexTreeEllipsoids), maxFiltration(simplexTreeRips)) + 0.1
    xAxisEnd = 5
    print('Done')
    print('Plotting ellipsoid barcode... ', end='', flush=True)
    barcodePlotting.plot_persistence_barcode(barcodeEllipsoids, inf_delta=0.5, axes=axBarE, fontsize=12,\
                                             axis_start = -0.1, infinity = xAxisEnd) #(0.1 + xAxisEnd) /2 ) # todo: put this back to xAxisEnd (without /2)
    print('Done.')
    axBarE.set_title('Ellipsoid barcode', fontsize=12)

    print('Plotting Rips barcode... ', end='', flush=True)
    barcodePlotting.plot_persistence_barcode(barcodeRips, inf_delta=0.5, axes=axBarR, fontsize=12,\
                                             axis_start = -0.1, infinity = xAxisEnd) #(0.1 + xAxisEnd))
    print('Done.')
    axBarR.set_title('Rips barcode', fontsize=12)
    # /alt


    # plotting the vertical lines at r
    # if 'rValues' in kwargs:
    #     rValues = kwargs['rValues']
    #     for rValue in rValues:
    #         axBarE.axvline(x = rValue, color='gray', linewidth=0.5, linestyle='dashed')
    
    if 'savePlot' in kwargs and kwargs['savePlot'] is True:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print('Plot saved to file.')

    if 'showPlot' in kwargs and kwargs['showPlot'] is True:
        plt.show()

def visualisationFromFile(filename, rPlot=0.6, plotEllipsoids=False):
    vars = loadVarsFromFile(filename)

    if plotEllipsoids is True:
        visualisation(points = vars['points'],\
                    ellipsoidList = vars['ellipsoidList'], rPlot = rPlot, \
                    simplexTreeEllipsoids = vars['simplexTreeEllipsoids'], \
                    barcodeEllipsoids = vars['barcodeEllipsoids'], \
                    barcodeRips = vars['barcodeRips'], \
                    showPlot = True, \
                    savePlot = False, \
                    rValues = vars['rValues']
                    )
    else:
        visualisation(points = vars['points'],\
                    #ellipsoidList = vars['ellipsoidList'], rPlot = rPlot, \
                    simplexTreeEllipsoids = vars['simplexTreeEllipsoids'], \
                    barcodeEllipsoids = vars['barcodeEllipsoids'], \
                    barcodeRips = vars['barcodeRips'], \
                    showPlot = True, \
                    savePlot = False, \
                    rValues = vars['rValues']
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

def calculateBottleeckDistance(barcode1, barcode2, dim):
    # create the numpy array in dimension dim
    npBarcode1 = np.array()
    npBarcode2 = np.array()
    for line in barcode1:
        npBarcode1[line[0]].append(line[1])
    for line in barcode2:
        npBarcode2[line[0]].append(line[1])

    bottleneckDistance = [gd.bottleneck_distance(i,j) for i,j in zip(npBarcode1, npBarcode2)]

    return bottleneckDistance
    # TODO: finish this
    pass


def main():
    tStart = time.time()

    ###### User input ######
    boolSaveData = True
    boolShowPlot = True
    boolSavePlot = True
    rPlot = 0.7            # if rPlot = 0, ellipses won't be plotted
    expansionDim = 3
    # -------------------- #
    dim = 2              # dimension of the ambient space
    rStart = 0.1
    rEnd = 4
    rStep = 0.1
    rValues = np.arange(rStart, rEnd, rStep)
    nbhdSize = 5         # number of points for doing PCA
    nPts = 500            # number of data points
    axesRatios = np.array([2,1,1])
    # --------------------- #
    #   Specifying points   #

    dataType = 'unknown_'

    # 2d circle, ellipse:
    # variation = 0.2
    # points = createData(nPts,'circle', variation = variation)
    # dataType = 'circle-var' + str(variation) + '_'

    # points = createData(nPts,'circle', outliers = False)

    # torus in R3:
    # points = createData(nPts, 'torus')
    # dataType = 'torus_'

    # more advanced circle, annulus (also in higher dimensions):
    # points = shapes.sample_from_sphere(n=nPts, d=(dim-1), seed=0)

    # figure eight:
    # points = figure_eight.figure_eight(nPts, 1, 0.1)
    # dataType = 'figureEight_'

    # to read in a mesh from an OFF file:
    # meshFileName = '61.off'
    # points = readOFF('data' + meshFileName) # warning: 61.off is a mesh with 1k+ vertices.
    # dataType = 'mesh_' + str(meshFileName)
    
    # to read in the CycloOctane dataset:
    points = np.asarray(\
        scipy.io.loadmat('pointsCycloOctane.mat')['pointsCycloOctane']) # size: 6040 x 24
    dataType = 'cyclooctane'
    points = points[0:nPts,:]
    # scipy.io.savemat('data/cyclooctane.mat', {'points': points})
    dim = 24
    nbhdSize = 26
    axesRatios = np.pad(axesRatios, (0,dim - len(axesRatios)), constant_values=1) # creates an array of length dim
    ########################

    rStart = rValues[0]
    rEnd = rValues[-1]
    rStep = 0 if len(rValues) == 1 else abs(rValues[1] - rValues[0])
    nPts = len(points)
    dim = len(points[0])
    if len(axesRatios) != dim:
        raise ValueError('The number of ellipsoid axes does not correspond to the ambient dimension.')
    
    importantParameters = f'{dataType=}' + '_' + f'{nPts=}' + '_' + f'{nbhdSize=}' + f'{axesRatios}'
    dataPlotFilename='data/ellipsoids_'+importantParameters+datetime.now().strftime("_%Y%m%d_%H%M%S")
    
    if (boolShowPlot or boolSavePlot) and (rPlot not in rValues):
        print('\nWarning: the simplex tree plot may be inaccurate since the calculations are ' \
            +'not performed for the chosen value of rPlot. To fix this, make sure that ' \
            +'rPlot is in np.arange(rStart,rEnd,rStep).')

    tStartEllipsoids = time.time()    # for testing performace

    print('Creating KD tree... ', end='', flush=True)
    kdTree = spatial.KDTree(points)
    print('Done.')
    print('Fitting ellipsoids... ', end='', flush=True)
    ellipsoidList = fitEllipsoids(dim, kdTree, nbhdSize, axesRatios)
    print('Done.')
    longestEllipsoidAxis = max(ellipsoid.axesLengths[0] for ellipsoid in ellipsoidList)
    queryRadius = 2*longestEllipsoidAxis
    # original:
    # try:
    #     # simplexTreeEllipsoids = generateEllipoidSimplexTree(kdTree, ellipsoidList, queryRadius, \
    #     #                                               filtrationValues = rValues)
    #     simplexTreeEllipsoids = generateEllipoidSimplexTree2(kdTree, ellipsoidList, axesRatios)
    # except np.linalg.LinAlgError:
    #     print(f'{points=}')
    #     print("\nError: Attempting to add an edge in Ellipsoid Simplex from a vertex to itself. \n" + 
    #             "Reason: Two of the ellipsoids are the same because two points from " + \
    #             "the data set have the same neighbourhood. Please change the " + \
    #             "neighbourhood size (nbhdSize) and try again.")
    #     return 0
    #debug:
    print('Generating ellipsoid simplex tree... ', end='', flush=True)
    simplexTreeEllipsoids = generateEllipoidSimplexTree2(kdTree, ellipsoidList, axesRatios)
    print('Done.')
    # ----

    simplexTreeEllipsoidsExpanded = simplexTreeEllipsoids.copy()
    print('Expanding the ellipsoid simplex tree... ', end='', flush=True)
    simplexTreeEllipsoidsExpanded.expansion(expansionDim) # expands the simplicial complex to include 
                                        # dim-dimensional simplices whose 1-skeleton is in simplexTree
    print('Done.')
    print('Calculating the persistence of the expanded simplex tree... ', end='', flush=True)
    barcodeEllipsoids = simplexTreeEllipsoidsExpanded.persistence()
    print('Done.')
    
    tEndEllipsoids = time.time()
    tStartRips = time.time()

    print('Creating the Rips complex... ', end='', flush=True)
    ripsComplex = gd.RipsComplex(points=points)
    print('Done.')
    print('Creating the Rips simplex tree... ', end='', flush=True)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=expansionDim)
    print('Done.')
    print('Calculating the Rips barcode... ', end='', flush=True)
    barcodeRips = simplexTreeRips.persistence()
    print('Done.')

    tEndRips = time.time()

    tEllipsoids = tEndEllipsoids - tStartEllipsoids
    tRips = tEndRips - tStartRips
    tRatio = tEllipsoids / tRips

    print(f'{nPts=}')
    print(f'{tRatio=}')
    print(f'{nbhdSize}')
    #print(simplexTreeRips.persistence_intervals_in_dimension(1))

    tEnd = time.time()

    print('\nThe total execution time is ' + str(tEnd-tStart) + '\n')
    
    if boolSaveData is True:
        print('Saving data to file... ', end='', flush=True)
        # debug:
        # barcodeRips1 = barcodeRips
        # barcodeRips = [0,[0,1]]

        dictOfVars = {
            'dim': dim,
            'rStart': rStart,
            'rEnd': rEnd,
            'rStep': rStep,
            'rValues': rValues,
            'nbhdSize': nbhdSize,
            'nPts': nPts,
            'points': points,
            'ellipsoidList': ellipsoidList,
            'simplexTreeEllipsoids': simplexTreeEllipsoids,
            # 'simplexTreeRips': simplexTreeRips,
            'barcodeEllipsoids': barcodeEllipsoids,
            'barcodeRips': barcodeRips
        }
        saveVarsToFile(dictOfVars, \
                       filename = dataPlotFilename + '.json')
        
        # debug:
        # barcodeRips = barcodeRips1

    if (boolShowPlot or boolSavePlot) is True:
        visualisation(points = points,\
                        ellipsoidList = ellipsoidList, rPlot = rPlot, \
                        simplexTreeEllipsoids = simplexTreeEllipsoids, \
                        simplexTreeRips = simplexTreeRips, \
                        barcodeEllipsoids = barcodeEllipsoids, \
                        barcodeRips = barcodeRips, \
                        showPlot = boolShowPlot, \
                        savePlot = boolSavePlot, \
                        filename = dataPlotFilename + '.png', \
                        rValues = rValues)
        

if __name__ == "__main__":
    main()




    


