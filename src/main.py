import numpy as np
import time
from datetime import datetime
from ellipsoidSimplexTree import generateEllipsoidSimplexTree2
from ellipsoidSimplexTree import generateEllipsoidSimplexTree3
from ellipsoidSimplexTree import generateEllipsoidSimplexTree4
from importPoints import importPoints
from ripsSimplexTree import generateRipsSimplexTree
from readWriteData import saveVarsToFile
from utils import expandTreeAndCalculateBarcode
from utils import maxFiltration
from utils import padAxesRatios
from visualisation import visualisation
import gudhi as gd

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
    expansionDim = 2                            # if expansionDim=3, it takes ages for 300+ points
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

    if 'dim' in kwargs:
        points = importPoints(nPts,dataType=dataType,dim=dim)
    else:
        points = importPoints(nPts, dataType=dataType) 

    # depending on the type of data used, it might be necessary to read these in again
    nPts = len(points)
    dim = len(points[0])

    # filling out axesRatios in case only the first few terms are specified
    axesRatios = padAxesRatios(axesRatios, dim)

    importantParameters = 'dataType=' + f'{dataType}' + '_' + f'{nPts=}' + '_' + f'{nbhdSize=}' + f'{axesRatios=}'
    dataPlotFilename='data/ellipsoids_'+importantParameters+datetime.now().strftime("_%Y%m%d_%H%M%S")

    ######### Ellipsoids #########
    tStart = time.time()
    tStartEllipsoidsST = time.time()    # for testing performace
    [simplexTreeEllipsoids, ellipsoidList] = generateEllipsoidSimplexTree4(points, nbhdSize, axesRatios)
    tEndEllipsoidsST = time.time()

    if saveData is True:
        print('Saving data to file... ', end='', flush=True)
        dictOfVars = {
            'dim': dim,
            'expansionDim': expansionDim,
            'nbhdSize': nbhdSize,
            'nPts': nPts,
            'points': points,
            'ellipsoidList': ellipsoidList,
            'simplexTreeEllipsoids': simplexTreeEllipsoids
        }
        saveVarsToFile(dictOfVars, \
                        filename = dataPlotFilename + '.json')

    tStartEllipsoidsBarcode = time.time()
    barcodeEllipsoids = expandTreeAndCalculateBarcode(simplexTreeEllipsoids, expansionDim, collapseEdges=collapseEdges)
    tEndEllipsoidsBarcode = time.time()

    if saveData is True:
        print('Saving data to file... ', end='', flush=True)
        dictOfVars = {
            'dim': dim,
            'expansionDim': expansionDim,
            'nbhdSize': nbhdSize,
            'nPts': nPts,
            'points': points,
            'ellipsoidList': ellipsoidList,
            'simplexTreeEllipsoids': simplexTreeEllipsoids,
            'barcodeEllipsoids': barcodeEllipsoids 
        }
        saveVarsToFile(dictOfVars, filename = dataPlotFilename + '.json')

    ########### Rips ###########

    tStartRipsST = time.time()
    simplexTreeRips = generateRipsSimplexTree(points)
    tEndRipsST = time.time()

    if saveData is True:
        print('Saving data to file... ', end='', flush=True)
        dictOfVars = {
            'dim': dim,
            'expansionDim': expansionDim,
            'nbhdSize': nbhdSize,
            'nPts': nPts,
            'points': points,
            'ellipsoidList': ellipsoidList,
            'simplexTreeEllipsoids': simplexTreeEllipsoids,
            'barcodeEllipsoids': barcodeEllipsoids,
            'simplexTreeRips': simplexTreeRips 
        }
        saveVarsToFile(dictOfVars, filename = dataPlotFilename + '.json')

    tStartRipsBarcode = time.time()
    barcodeRips = expandTreeAndCalculateBarcode(simplexTreeRips, expansionDim, collapseEdges=collapseEdges)
    tEndRipsBarcode = time.time()

    if saveData is True:
        print('Saving data to file... ', end='', flush=True)
        dictOfVars = {
            'dim': dim,
            'expansionDim': expansionDim,
            'nbhdSize': nbhdSize,
            'nPts': nPts,
            'points': points,
            'ellipsoidList': ellipsoidList,
            'simplexTreeEllipsoids': simplexTreeEllipsoids,
            'barcodeEllipsoids': barcodeEllipsoids,
            'simplexTreeRips': simplexTreeRips,
            'barcodeRips': barcodeRips 
        }
        saveVarsToFile(dictOfVars, filename = dataPlotFilename + '.json')

    tEllipsoids = tEndEllipsoidsST - tStartEllipsoidsST + tEndEllipsoidsBarcode - tStartEllipsoidsBarcode
    tRips = tEndRipsST - tStartRipsST + tEndRipsBarcode - tStartRipsBarcode
    tRatio = tEllipsoids / tRips
    tEnd = time.time()
    tTotal = tEnd - tStart

    print('\nThe total execution time is ' + str(tEnd-tStart) + '\n')

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

    axesRatios = np.array([2,1])
    expansionDim = 2
    nPtsValues = np.array([500, 200])

    # circles
    # dim = 3
    # nPtsValues = np.array([100, 300, 500, 1000, 1500])
    # dataType = 'sphere'

    # for nPts in nPtsValues:
    #     main(nPts=nPts, dataType=dataType, collapseEdges=True, axesRatios=axesRatios, dim=3)

    # cyclooctane
    # nPtsValues = np.array([100, 300, 500, 1000, 1500, 2000])
    dataType = 'cyclooctane'
    for nPts in nPtsValues:
        main(nPts=nPts, dataType=dataType, collapseEdges=True, axesRatios=axesRatios, expansionDim=expansionDim)




    


