from ellipsoids import visualisationFromFile

###### User input ######
filenameLoad = \
    "data/ellipsoids_dataType='cyclooctane'_nPts=100_nbhdSize=26[2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]_20231025_142315.json"
    # 'data/ellipsoids_dataType=cyclooctane_nPts=199_nbhdSize=26[2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]_20231018_135436.json'
    #'data/ellipsoids_nPts=50_rStep=0.49999999999999994_nbhdSize=5_20230816_124230.json'

rPlot = 2
nBarsDim0 = 5
nBarsDim1 = 10
nBarsDim2 = 4

########################

visualisationFromFile(filenameLoad, \
                      nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2, \
                      rPlot=rPlot, plotEllipsoids=False)
