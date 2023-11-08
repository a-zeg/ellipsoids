import numpy as np
import time
from datetime import datetime
from ellipsoidSimplexTree import generateEllipoidSimplexTree2
from importPoints import importPoints
from ripsSimplexTree import generateRipsSimplexTree
from readWriteData import saveVarsToFile
from utils import expandTreeAndCalculateBarcode
from utils import maxFiltration
from utils import padAxesRatios
from visualisation import visualisation

def main(**kwargs):
    ###### Default parameter values ######
    saveData = True
    showPlot = False
    savePlot = True
    plotPoints = False
    plotBarcodes = True
    collapseEdges = True
    rPlot = 1.2                                 # if rPlot = 0, ellipses won't be plotted
    dim = 3                                     # dimension of the ambient space
    nbhdSize = 5                                # number of points for doing PCA
    nPts = 100                                  # number of data points
    expansionDim = 3        
    axesRatios = np.array([2,1])
    dataType = 'circle'
    #######################################
    if 'saveData' in kwargs:
        saveData = kwargs['saveData']
    if 'showPlot' in kwargs:
        showPlot = kwargs['showPlot']
    if 'savePlot' in kwargs:
        savePlot = kwargs['savePlot']
    if 'plotPoints' in kwargs:
        plotPoints = kwargs['plotPoints']
    if 'plotBarcodes' in kwargs:
        plotBarcodes = kwargs['plotBarcodes']
    if 'collapseEdges' in kwargs:
        collapseEdges = kwargs['collapseEdges']
    if 'rPlot' in kwargs:
        rPlot = kwargs['rPlot']
    if 'dim' in kwargs:
        dim = kwargs['dim']
    if 'nbhdSize' in kwargs:
        nbhdSize = kwargs['nbhdSize']
    if 'nPts' in kwargs:
        nPts = kwargs['nPts']
    if 'expansionDim' in kwargs:
        expansionDim = kwargs['expansionDim']
    if 'axesRatios' in kwargs:
        axesRatios = kwargs['axesRatios']
    if 'dataType' in kwargs:
        dataType = kwargs['dataType']
    # ----------------------------------- #

    if dataType == 'cyclooctane':
        nbhdSize = 26

    points = importPoints(nPts, dataType=dataType) 

    # depending on the type of data used, it might be necessary to read these in again
    nPts = len(points)
    dim = len(points[0])

    # filling out axesRatios in case only the first few terms are specified
    axesRatios = padAxesRatios(axesRatios, dim)

    tStart = time.time()
    tStartEllipsoids = time.time()    # for testing performace
    [simplexTreeEllipsoids, ellipsoidList] = generateEllipoidSimplexTree2(points, nbhdSize, axesRatios)
    barcodeEllipsoids = expandTreeAndCalculateBarcode(simplexTreeEllipsoids, expansionDim, collapseEdges=collapseEdges)
    tEndEllipsoids = time.time()

    tStartRips = time.time()
    simplexTreeRips = generateRipsSimplexTree(points)
    barcodeRips = expandTreeAndCalculateBarcode(simplexTreeRips, expansionDim, collapseEdges=collapseEdges)
    tEndRips = time.time()

    tEllipsoids = tEndEllipsoids - tStartEllipsoids
    tRips = tEndRips - tStartRips
    tRatio = tEllipsoids / tRips
    tEnd = time.time()
    tTotal = tEnd - tStart

    print('\nThe total execution time is ' + str(tEnd-tStart) + '\n')
    importantParameters = 'dataType=' + f'{dataType}' + '_' + f'{nPts=}' + '_' + f'{nbhdSize=}' + f'{axesRatios}'
    dataPlotFilename='data/ellipsoids_'+importantParameters+datetime.now().strftime("_%Y%m%d_%H%M%S")

    if saveData is True:
        print('Saving data to file... ', end='', flush=True)
        dictOfVars = {
            'dim': dim,
            'expansionDim': expansionDim,
            'nbhdSize': nbhdSize,
            'nPts': nPts,
            'totalExecutionTime': tTotal,
            'timeEllipsoids/tRips': tRatio,
            'points': points,
            'ellipsoidList': ellipsoidList,
            'simplexTreeEllipsoids': simplexTreeEllipsoids,
            'simplexTreeRips': simplexTreeRips,
            'barcodeEllipsoids': barcodeEllipsoids,
            'barcodeRips': barcodeRips
        }
        saveVarsToFile(dictOfVars, \
                       filename = dataPlotFilename + '.json')
        
    if (showPlot or savePlot) is True:
        xAxisEnd = max(maxFiltration(simplexTreeEllipsoids), maxFiltration(simplexTreeRips)) + 0.1
        if dim > 3: plotPoints = False
        
        if plotPoints is True:                  # plotting points, ellipsoids, and the barcodes
            visualisation(
                        points = points,\
                        xAxisEnd = xAxisEnd,\
                        expansionDim = expansionDim, \
                        ellipsoidList = ellipsoidList, rPlot = rPlot, \
                        simplexTreeEllipsoids = simplexTreeEllipsoids, \
                        simplexTreeRips = simplexTreeRips, \
                        barcodeEllipsoids = barcodeEllipsoids, \
                        barcodeRips = barcodeRips, \
                        plotPoints = plotPoints, \
                        plotBarcodes = plotBarcodes, \
                        showPlot = showPlot, \
                        savePlot = savePlot, \
                        filename = dataPlotFilename + '.png')
        else:                                       # plotting just the barcodes
            visualisation(  
                        points = points,\
                        xAxisEnd = xAxisEnd,\
                        barcodeEllipsoids = barcodeEllipsoids, \
                        barcodeRips = barcodeRips, \
                        showPlot = showPlot, \
                        savePlot = savePlot, \
                        filename = dataPlotFilename + '.png') 
        

if __name__ == "__main__":
    nPtsValues = np.array([100])#, 200, 300, 400, 500, 800, 1000])
    dataType = 'circle'
    for nPts in nPtsValues:
        main(nPts=nPts, dataType=dataType, collapseEdges=True)




    


