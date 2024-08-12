from visualisation import visualisationFromFile
from datetime import datetime
from data_handling import read_variables
from visualisation import plotDataPoints
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt

###### User input ######
filename = \
    "data/shapes/ellipsoids_data_type=figure_eight_n_pts=50_nbhd_size=5_axes_ratios=[3 1]_seed=0__20240430_154606.json"

# folder = "data/shapes"


vars = read_variables(filename)

points = vars['points']
points = np.asarray(points)

# fig, axs = plt.figure()
for point in points:
    plt.scatter(*point, c='k')
    # plt.set_title('Data (%d points)' %len(points), fontsize=12)

plt.show()

