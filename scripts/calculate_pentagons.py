import os 
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath('.'))

from ellipsoids.topological_computations import calculate_ellipsoid_barcode
from ellipsoids.topological_computations import calculate_rips_barcode
from ellipsoids.topological_computations import reduceBarcode
from ellipsoids.data_handling import set_filename_parameters
from ellipsoids.data_handling import generate_filename
from ellipsoids.data_handling import save_variables
from ellipsoids.data_handling import read_variables
from ellipsoids.data_handling import get_paths_of_files_in_a_folder
from ellipsoids.data_handling import import_maxmin_mat
from ellipsoids.visualisation import plot_barcode



def calculate_ellipsoid_complex(datasets_folder, output_folder):

    paths = get_paths_of_files_in_a_folder(datasets_folder)

    for path in paths:

        points = import_maxmin_mat(path)

        data_type_params={}
        n_pts = len(points)

        axes_ratios_all = [np.array([3,3,1])]

        expansion_dim = 2
        nbhd_sizes = [6, 10, 15]

        for axes_ratios in axes_ratios_all:

            for nbhd_size in nbhd_sizes:
                filename_parameters = set_filename_parameters('pentagons', n_pts, nbhd_size, axes_ratios, data_type_params)
                save_filename = generate_filename(filename_parameters, folder=output_folder)

                barcode_ellipsoids, simplex_tree_ellipsoids, ellipsoid_list, _ = calculate_ellipsoid_barcode(points, nbhd_size, axes_ratios, expansion_dim=expansion_dim)
                barcode_rips, simplex_tree_rips, _ = calculate_rips_barcode(points, expansion_dim=expansion_dim)

                dict_vars = {
                    'points' : points,
                    'barcode_ellipsoids' : barcode_ellipsoids,
                    'barcode_rips' : barcode_rips,
                    'simplex_tree_ellipsoids' : simplex_tree_ellipsoids,
                    'simplex_tree_rips' : simplex_tree_rips,
                    'ellipsoid_list' : ellipsoid_list
                }

                save_variables(dict_vars, save_filename, timestamp=True)



def plot_from_file(folder):

    n_bars_dim0 = 5
    n_bars_dim1 = 30

    paths = get_paths_of_files_in_a_folder(folder, extension='.json')

    for path in paths:
        dict_vars = read_variables(path)

        points = dict_vars['points']
        points = np.asarray(points)
        barcode_ellipsoids = dict_vars['barcode_ellipsoids']
        barcode_rips = dict_vars['barcode_rips']

        fig = plt.figure()
        ax_bar_E = fig.add_subplot(121)
        ax_bar_R = fig.add_subplot(122) 

        reduced_barcode_ellipsoids, maxBarEndEllipsoids = reduceBarcode( \
                            barcode_ellipsoids, \
                            nBarsDim0=n_bars_dim0, \
                            nBarsDim1=n_bars_dim1)
        reduced_barcode_rips, maxBarEndRips = reduceBarcode( \
                                    barcode_rips, \
                                    nBarsDim0=n_bars_dim0, \
                                    nBarsDim1=n_bars_dim1)
        
        xAxisEnd = max(maxBarEndEllipsoids, maxBarEndRips) * 1.1
        plot_barcode(reduced_barcode_ellipsoids, inf_delta=0.5, axes=ax_bar_E, fontsize=12,\
                                        axis_start = -0.1, infinity = xAxisEnd, max_intervals=100)
        ax_bar_E.set_title('Ellipsoids barcode')

        plot_barcode(reduced_barcode_rips, inf_delta=0.5, axes=ax_bar_R, fontsize=12,\
                                        axis_start = -0.1, infinity = xAxisEnd, max_intervals=100)
        ax_bar_R.set_title('Rips barcode')
        
        fig.set_size_inches(12, 4)
        plt.tight_layout()
        output_name = os.path.splitext(path)[0] + '.png'
        fig.savefig(output_name, dpi=300)

        print('File saved to ' + output_name)



if __name__ == '__main__':

    output_folder = 'data/pentagons'
    datasets_folder = 'datasets/pentagons/maxmin'

    ## uncomment for calculation:
    # calculate_ellipsoid_complex(datasets_folder, output_folder)

    ## uncomment for plotting:
    plot_from_file(output_folder) 