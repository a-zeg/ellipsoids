from utils import recalculateBarcodesFromFile

filenameLoad = \
    "data/ellipsoids_dataType='cyclooctane'_nPts=500_nbhdSize=26[2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]_20231025_011741.json"

recalculateBarcodesFromFile(filenameLoad, \
                            expansionDim=3, collapseEdges=True)
