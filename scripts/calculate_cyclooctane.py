import os
import numpy as np
import sys

sys.path.append(os.path.abspath('.'))

from ellipsoids.data_handling import get_paths_of_files_in_a_folder
from ellipsoids.data_handling import import_maxmin_mat
from ellipsoids.data_handling import set_filename_parameters
from ellipsoids.data_handling import generate_filename
from ellipsoids.data_handling import calculate_and_save_ellipsoids_and_rips_data

data_type = 'cyclooctaneMaxmin'

axes_ratios_all = np.array([[3,3,1]])
expansion_dim = 3
nbhd_sizes = [26]
folder = 'datasets/cyclooctane/pca_analysis'
paths = get_paths_of_files_in_a_folder(folder, extension='.mat')

for path in paths:
    for axes_ratios in axes_ratios_all:
        for nbhd_size in nbhd_sizes:

            # Import points
            points = import_maxmin_mat(path=path)

            # Generate the filename for storing the variables 
            n_pts = len(points)
            filename_parameters = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
            save_filename = generate_filename(filename_parameters)

            calculate_and_save_ellipsoids_and_rips_data(points, nbhd_size, axes_ratios, expansion_dim, save_filename)
