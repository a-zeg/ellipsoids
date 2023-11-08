from visualisation import visualisationFromFile

###### User input ######
filenameLoad = \
    "../data/ellipsoids_dataType='cyclooctane'_nPts=100_nbhdSize=26[2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]_20231025_142315-barcodes_expansionDim=3_20231030_204010.json"

rPlot = 2
nBarsDim0 = 5
nBarsDim1 = 10
nBarsDim2 = 4

########################

visualisationFromFile(filenameLoad, \
                      nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2, \
                      rPlot=rPlot, plotEllipsoids=False)
