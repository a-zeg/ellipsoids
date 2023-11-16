from utils import recalculateBarcodesFromFile

filenameLoad = \
    "data/300pts_cyclooctane_sphere_expansionDim3/ellipsoids_dataType=sphere_nPts=300_nbhdSize=5axesRatios=array([2, 1, 1])_20231116_000450.json"

recalculateBarcodesFromFile(filenameLoad, \
                            expansionDim=3, collapseEdges=False)
