#
# reduce barcode
# recalculate barcodes from file - maybe for this one should separate it into calculate barcodes from simplex tree and read from file?

from readWriteData import loadVarsFromFile
from readWriteData import saveVarsToFile
import numpy as np
import gudhi as gd
from datetime import datetime


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
    # ###
    # print('Number of simplices is ' + str(simplexTreeExpanded.num_simplices()))
    # ###
    print('Calculating the barcode of the expanded tree... ', end='', flush=True)
    barcode = simplexTreeExpanded.persistence()
    print('Done.')

    return barcode

def maxFiltration(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    return max(splx[1] for splx in simplexList)

def reduceBarcode(barcode, nBarsDim0 = 1, nBarsDim1 = 0, nBarsDim2 = 0):
    # return only nBarsDimk longest bars in each dimension k
    reducedBarcode = []
    maxBarEnd = 0
    for bar in barcode:
        if bar[0] == 0 and nBarsDim0 > 0:
            reducedBarcode.append(bar)
            nBarsDim0 = nBarsDim0 - 1
            barEnd = bar[1][1]
            if barEnd != float('inf') and barEnd > maxBarEnd:
                maxBarEnd = barEnd
        elif bar[0] == 1 and nBarsDim1 > 0:
            reducedBarcode.append(bar)
            nBarsDim1 = nBarsDim1 - 1
            barEnd = bar[1][1]
            if barEnd != float('inf') and barEnd > maxBarEnd:
                maxBarEnd = barEnd
        elif bar[0] == 2 and nBarsDim2 > 0:
            reducedBarcode.append(bar)
            nBarsDim2 = nBarsDim2 - 1
            barEnd = bar[1][1]
            if barEnd != float('inf') and barEnd > maxBarEnd:
                maxBarEnd = barEnd
    
    return reducedBarcode, maxBarEnd

def calculateBottleeckDistance(barcode1, barcode2, dim):
    npBarcode1 = np.array()
    npBarcode2 = np.array()
    for line in barcode1:
        npBarcode1[line[0]].append(line[1])
    for line in barcode2:
        npBarcode2[line[0]].append(line[1])

    bottleneckDistance = [gd.bottleneck_distance(i,j) for i,j in zip(npBarcode1, npBarcode2)]

def recalculateBarcodesFromFile(filename, expansionDim=2, collapseEdges=False):
    print('Reading in the variables... ', end='', flush=True)
    vars = loadVarsFromFile(filename)
    if 'expansionDim' in vars and vars['expansionDim'] == expansionDim:
        print('The original barcode is already expanded to the specified dimension.')
        return None
    if 'simplexTreeEllipsoids' in vars:
        simplexTreeEllipsoids = vars['simplexTreeEllipsoids']
    else: simplexTreeEllipsoids = gd.SimplexTree()
    if 'simplexTreeRips' in vars:
        simplexTreeRips = vars['simplexTreeRips']
    else: simplexTreeRips = gd.SimplexTree()
    print('Done.')
    
    barcodeEllipsoids = expandTreeAndCalculateBarcode(simplexTreeEllipsoids, expansionDim, collapseEdges=collapseEdges)
    barcodeRips = expandTreeAndCalculateBarcode(simplexTreeRips, expansionDim, collapseEdges=collapseEdges)
    
    dictOfVars = {
        'originalFile': filename, 
        'expansionDim': expansionDim,
        'barcodeEllipsoids': barcodeEllipsoids,
        'barcodeRips': barcodeRips
    }

    filename =  filename[:filename.rfind('.')] + '-barcodes_expansionDim=' + f'{expansionDim}' + datetime.now().strftime("_%Y%m%d_%H%M%S") + '.json'
    saveVarsToFile(dictOfVars, filename=filename)

def padAxesRatios(axesRatios, dim):
    ''' For high dimensional ellipsoids, it is enough for the user to specify 
    the first few axes. This function will set the remaining axes to 1.'''
    if dim > len(axesRatios):
        return np.pad(axesRatios, (0,dim - len(axesRatios)), constant_values=1) # creates an array of length dim
    else: 
        return axesRatios[0:dim]

