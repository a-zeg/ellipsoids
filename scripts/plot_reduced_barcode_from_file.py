import os
import sys

from datetime import datetime

import matplotlib.pyplot as plt

sys.path.append(os.path.abspath('.'))

from ellipsoids.visualisation import plot_barcode
from ellipsoids.topological_computations import reduceBarcode
from ellipsoids.data_handling import read_variables

###### User input ######
filename = \
    "data/shapes/ellipsoids_data_type=figure_eight_n_pts=150_nbhd_size=5_axes_ratios=[3 1]_seed=0__20240430_162531.json"
    # "data/shapes/ellipsoids_data_type=annulus_n_pts=10_nbhd_size=5_axes_ratios=[3 1]_seed=0__20240430_154603.json"

folder = "data/shapes"

rPlot = 0.2
nBarsDim0 = 10
nBarsDim1 = 10
nBarsDim2 = 0

savePlot = True
drawEllipsoidsSimplexTree = False
drawEllipsoids = False

plotDensity = False
persistenceDim = None # if None, bBarsDimn will be ignored and bars in all dimensions will be plotted

vars = read_variables(filename)

barcodeEllipsoids = vars['barcode_ellipsoids']
barcodeRips = vars['barcode_rips']

reducedBarcodeEllipsoids, maxBarEndEllipsoids = reduceBarcode( \
                            barcodeEllipsoids, \
                            nBarsDim0=nBarsDim0, \
                            nBarsDim1=nBarsDim1, \
                            nBarsDim2=nBarsDim2)
reducedBarcodeRips, maxBarEndRips = reduceBarcode( \
                            barcodeRips, \
                            nBarsDim0=nBarsDim0, \
                            nBarsDim1=nBarsDim1, \
                            nBarsDim2=nBarsDim2)

barcodeEllipsoids = reducedBarcodeEllipsoids
barcodeRips = reducedBarcodeRips
xAxisEnd = max(maxBarEndEllipsoids, maxBarEndRips) * 1.1

fig, axes = plt.subplots()
# fig.set_size_inches(10,2)
# axes.set_ylim([-1,10])

plot_barcode(barcodeRips, inf_delta=0.5, axes=axes, fontsize=12,\
            axis_start = -0.1, infinity = xAxisEnd, max_intervals=100) 

plt.show()

