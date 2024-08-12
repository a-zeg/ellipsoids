
from main import set_filename_parameters
from main import generate_filename
from main import calculate_and_save_ellipsoids_and_rips_data
import data_handling
import numpy as np
import os


def calculate_shapes():
    folder = os.path.join('datasets', 'shapes', 'shapes_maxmin')
    paths = data_handling.get_paths_of_files_in_a_folder(folder, extension='.mat')

    for path in paths:

        if 'seed=0' in path:
            # Import points
            points = data_handling.import_maxmin_mat(path=path)
            data_type_params = {'seed': 0}

            if 'circle' in path:
                data_type = 'circle'
            elif 'sphere' in path:
                data_type = 'sphere'
            elif 'figure_eight' in path:
                data_type = 'figure_eight'
            elif 'annulus' in path:
                data_type = 'annulus'
            else:
                print('No valid data type found in the filename.')
                exit()
            
            # Generate the filename for storing the variables 
            n_pts = len(points)
            ambient_dim = len(points[0])

            if ambient_dim == 2:
                axes_ratios_all = np.array([[3,1]])
            elif ambient_dim == 3:
                axes_ratios_all = np.array([[3,3,1]])

            expansion_dim = 2
            nbhd_sizes = [5]

            for axes_ratios in axes_ratios_all:

                for nbhd_size in nbhd_sizes:
                    filename_parameters = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
                    save_filename = generate_filename(filename_parameters)

                    calculate_and_save_ellipsoids_and_rips_data(points, nbhd_size, axes_ratios, expansion_dim, save_filename)


if __name__ == '__main__':
    calculate_shapes()