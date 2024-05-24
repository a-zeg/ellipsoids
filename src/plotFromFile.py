from visualisation import visualisationFromFile
from datetime import datetime
from os import listdir
from os.path import isfile, join

###### User input ######
filename = \
    "data/ellipsoids_dataType=sphere_nPts=100_nbhdSize=10axesRatios=array([3, 3, 1])_seed=0_20240320_020017.json"

folder = "data/shapes"

rPlot = 0.2
nBarsDim0 = 10
nBarsDim1 = 10
nBarsDim2 = 10

savePlot = True
drawEllipsoidsSimplexTree = False
drawEllipsoids = False

plotDensity = False
persistenceDim = None # if None, bBarsDimn will be ignored and bars in all dimensions will be plotted

try:
    folder
except NameError:
    visualisationFromFile(filename, \
                      nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2, \
                      rPlot=rPlot, 
                      drawEllipsoids=drawEllipsoids,
                      drawEllipsoidsSimplexTree=drawEllipsoidsSimplexTree,
                      savePlot=savePlot,
                      plotDensity=plotDensity,
                      persistenceDim=persistenceDim)
else:
    filenames = [f for f in listdir(folder) if isfile(join(folder, f)) if f.endswith(".json")]
    # [print(folder + '/' + f) for f in filenames]
    # exit()
    [visualisationFromFile(folder + '/' + f, \
                          nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2, \
                          rPlot=rPlot,
                          drawEllipsoids=drawEllipsoids,
                          drawEllipsoidsSimplexTree=drawEllipsoidsSimplexTree,
                          savePlot=savePlot,
                          plotDensity=plotDensity,
                          persistenceDim=persistenceDim) for f in filenames]


########################

                    #   filename=filenameLoad)
