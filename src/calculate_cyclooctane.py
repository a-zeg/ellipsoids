from main import set_filename_parameters
from main import generate_filename
from main import calculate_and_save_ellipsoids_and_rips_data
import data_handling
import numpy as np

data_type = 'cyclooctaneMaxmin'

axes_ratios_all = np.array([[3,3,1]])
expansion_dim = 3
nbhd_sizes = [26]
folder = 'datasets/cyclooctane/pca_analysis'
paths = data_handling.get_paths_of_files_in_a_folder(folder, extension='.mat')

for path in paths:
    for axes_ratios in axes_ratios_all:
        for nbhd_size in nbhd_sizes:

            # Import points
            points, data_type_params = data_handling.importPoints(path=path)

            # Generate the filename for storing the variables 
            n_pts = len(points)
            filename_parameters = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
            save_filename = generate_filename(filename_parameters)

            calculate_and_save_ellipsoids_and_rips_data(points, nbhd_size, axes_ratios, expansion_dim, save_filename)
