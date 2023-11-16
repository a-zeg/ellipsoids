from visualisation import visualisationFromFile
from datetime import datetime

###### User input ######
filename = \
    "data/ellipsoids_dataType=sphere_nPts=500_nbhdSize=5axesRatios=array([2, 1, 1])_20231116_000655.json"

rPlot = 2
nBarsDim0 = 10
nBarsDim1 = 10
nBarsDim2 = 10

########################

visualisationFromFile(filename, \
                      nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2, \
                      rPlot=rPlot, 
                      plotEllipsoids=False,
                      savePlot=True)
                    #   filename=filenameLoad)
