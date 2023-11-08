from ellipsoidSimplexTree import Ellipsoid
import numpy as np
import matplotlib.pyplot as plt
import barcodePlotting
import gudhi as gd
from readWriteData import loadVarsFromFile
from utils import reduceBarcode


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

def plotDataPoints(points, axes=None):
    if axes is None:
        plt.scatter(points[:,0],points[:,1])
    else:
        axes.scatter(points[:,0],points[:,1])

def visualisation(**kwargs):
    print('Generating plots...')

    xAxisEnd = kwargs['xAxisEnd']
    barcodeEllipsoids = kwargs['barcodeEllipsoids']
    barcodeRips = kwargs['barcodeRips']

    listOfPlotPointsVars = \
    ['points', 'xAxisEnd',
                        'expansionDim',
                        'ellipsoidList',
                        'simplexTreeEllipsoids',
                        'simplexTreeRips',
                        'barcodeEllipsoids',
                        'barcodeRips',
                        'plotPoints',
                        'plotBarcodes',
                        'showPlot',
                        'savePlot',
                        'filename']
    
    if set(listOfPlotPointsVars).issubset(kwargs): 
        points = kwargs['points']
        dim = len(points[0])
        expansionDim = kwargs['expansionDim']
        simplexTreeEllipsoids = kwargs['simplexTreeEllipsoids']
        simplexTreeEllipsoids.expansion(expansionDim)

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
    else: filename = '../data/plotTest.png'

    print('Plotting ellipsoid barcode... ', end='', flush=True)
    barcodePlotting.plot_persistence_barcode(barcodeEllipsoids, inf_delta=0.5, axes=axBarE, fontsize=12,\
                                             axis_start = -0.1, infinity = xAxisEnd, max_intervals=100)
    print('Done.')
    axBarE.set_title('Ellipsoid barcode', fontsize=12)

    print('Plotting Rips barcode... ', end='', flush=True)
    barcodePlotting.plot_persistence_barcode(barcodeRips, inf_delta=0.5, axes=axBarR, fontsize=12,\
                                             axis_start = -0.1, infinity = xAxisEnd, max_intervals=100) #(0.1 + xAxisEnd))
    print('Done.')
    axBarR.set_title('Rips barcode', fontsize=12)
    
    if 'savePlot' in kwargs and kwargs['savePlot'] is True:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print('Plot saved to file.')

    if 'showPlot' in kwargs and kwargs['showPlot'] is True:
        plt.show()

def visualisationFromFile(\
        filename, \
        nBarsDim0=1, nBarsDim1=0, nBarsDim2=0, \
        rPlot=0.6, \
        plotEllipsoids=False):

    print('Reading in the variables... ', end='', flush=True)
    vars = loadVarsFromFile(filename)
    barcodeEllipsoids = vars['barcodeEllipsoids']
    barcodeRips = vars['barcodeRips']
    print('Done.')

    print('Calculating the reduced barcodes... ', end='', flush=True)
    reducedBarcodeEllipsoids, maxBarEndEllipsoids = reduceBarcode( \
                                barcodeEllipsoids, \
                                nBarsDim0=nBarsDim0, \
                                nBarsDim1=nBarsDim1, \
                                nBarsDim2=nBarsDim2)
    reducedBarcodeRips, maxBarEndRips = reduceBarcode( \
                                barcodeRips, \
                                nBarsDim0=nBarsDim0, \
                                nBarsDim1=nBarsDim1, \
                                nBarsDim2=nBarsDim2)
    print('Done.')

    barcodeEllipsoids = reducedBarcodeEllipsoids
    barcodeRips = reducedBarcodeRips
    xAxisEnd = max(maxBarEndEllipsoids, maxBarEndRips)

    print('Plotting...')
    if plotEllipsoids is True:
        simplexTreeEllipsoids = vars['simplexTreeEllipsoids']
        simplexTreeRips = vars['simplexTreeRips']
        points = vars['points']
        visualisation(points = points,\
                    ellipsoidList = vars['ellipsoidList'], \
                    rPlot = rPlot, \
                    simplexTreeEllipsoids = simplexTreeEllipsoids, \
                    barcodeEllipsoids = barcodeEllipsoids, \
                    barcodeRips = barcodeRips, \
                    xAxisEnd = xAxisEnd, \
                    showPlot = True, \
                    savePlot = False
                    )
    else:
        visualisation( \
                    xAxisEnd = xAxisEnd, \
                    barcodeEllipsoids = barcodeEllipsoids, \
                    barcodeRips = barcodeRips, \
                    showPlot = True, \
                    savePlot = False, \
                    )
