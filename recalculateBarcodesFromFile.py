from ellipsoids import recalculateBarcodesFromFile

###### User input ######
filenameLoad = \
    "data/ellipsoids_dataType='cyclooctane'_nPts=400_nbhdSize=26[2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]_20231025_005532.json"
    # 'data/ellipsoids_dataType=cyclooctane_nPts=199_nbhdSize=26[2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]_20231018_135436.json'
    #'data/ellipsoids_nPts=50_rStep=0.49999999999999994_nbhdSize=5_20230816_124230.json'
########################

recalculateBarcodesFromFile(filenameLoad, \
                            expansionDim=3)
