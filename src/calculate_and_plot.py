
from topological_computations import calculate_ellipsoid_barcode
from topological_computations import calculate_rips_barcode
from topological_computations import reduceBarcode
from data_handling import set_filename_parameters
from data_handling import generate_filename
import data_handling
import numpy as np
import visualisation
import os
import random
import matplotlib.pyplot as plt
import barcodePlotting
random.seed(10)


def calculate_ellipsoid_complex(data_type: str):

    output_folder='data/for_paper'

    n = 20

    # Import points
    if data_type == 'figure_eight':
        a = 2
        b = 0.5
        points = data_handling.figure_eight(n,a,b)
    elif data_type == 'circle':
        points = data_handling.sample_from_circle(n_pts=n)


    data_type_params={}
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
            save_filename = generate_filename(filename_parameters, folder=output_folder)

            barcode_ellipsoids, simplex_tree_ellipsoids, ellipsoid_list, t_ellipsoids = calculate_ellipsoid_barcode(points, nbhd_size, axes_ratios, expansion_dim=expansion_dim)
            barcode_rips, simplex_tree_rips, t_rips = calculate_rips_barcode(points, expansion_dim=expansion_dim)

            dict_vars = {
                'points' : points,
                'barcode_ellipsoids' : barcode_ellipsoids,
                'barcode_rips' : barcode_rips,
                'simplex_tree_ellipsoids' : simplex_tree_ellipsoids,
                'ellipsoid_list' : ellipsoid_list
            }

            data_handling.saveVarsToFile(dict_vars, save_filename, timestamp=True)

def plot_from_file(filename, r_plot=0.2):

    n_bars_dim0 = 5
    n_bars_dim1 = 3

    dict_vars = data_handling.loadVarsFromFile(filename)

    ellipsoid_list = dict_vars['ellipsoid_list']
    points = dict_vars['points']
    points = np.asarray(points)
    simplex_tree_ellipsoids = dict_vars['simplex_tree_ellipsoids']
    barcode_ellipsoids = dict_vars['barcode_ellipsoids']
    barcode_rips = dict_vars['barcode_rips']

    fig = plt.figure()
    ax_bar_E = fig.add_subplot(222)
    ax_bar_R = fig.add_subplot(224)
    ax_pts = fig.add_subplot(121)

    visualisation.plotEllipses(ellipseList=ellipsoid_list, r=r_plot, axes=ax_pts)
    visualisation.plotSimplexTree(points=points,simplexTree=simplex_tree_ellipsoids,r=r_plot,axes=ax_pts)
    ax_pts.set_aspect('equal', adjustable='box')
    n_pts = len(points)
    axes_ratio = ellipsoid_list[0].axesLengths / ellipsoid_list[0].axesLengths[-1]
    axes_ratio = axes_ratio.astype(int)
    axes_ratio = list(axes_ratio)
    plt.rcParams['text.usetex'] = True

    q = axes_ratio[0]
    ax_pts.set_title(f'{n_pts} points, ' + r'$q$=' + f'{q}, ' + r'$\varepsilon$' + f'={r_plot}')

    reduced_barcode_ellipsoids, maxBarEndEllipsoids = reduceBarcode( \
                        barcode_ellipsoids, \
                        nBarsDim0=n_bars_dim0, \
                        nBarsDim1=n_bars_dim1)
    reduced_barcode_rips, maxBarEndRips = reduceBarcode( \
                                barcode_rips, \
                                nBarsDim0=n_bars_dim0, \
                                nBarsDim1=n_bars_dim1)
    
    xAxisEnd = max(maxBarEndEllipsoids, maxBarEndRips) * 1.1
    
    barcodePlotting.plot_persistence_barcode(reduced_barcode_ellipsoids, inf_delta=0.5, axes=ax_bar_E, fontsize=12,\
                                    axis_start = -0.1, infinity = xAxisEnd, max_intervals=100)
    ax_bar_E.set_title('Ellipsoids barcode')

    barcodePlotting.plot_persistence_barcode(reduced_barcode_rips, inf_delta=0.5, axes=ax_bar_R, fontsize=12,\
                                    axis_start = -0.1, infinity = xAxisEnd, max_intervals=100)
    ax_bar_R.set_title('Rips barcode')
    
    plt.tight_layout()
    fig.set_size_inches(12, 8)
    output_name = os.path.splitext(filename)[0] + f'_{r_plot=}.png'
    fig.savefig(output_name, dpi=300)
    # plt.show()

    print('File saved to ' + output_name)


if __name__ == '__main__':
    # for calculation:
    # data_type = 'circle'
    # calculate_ellipsoid_complex(data_type)

    # # for plotting:
    filename = 'data/for_paper/ellipsoids_data_type=figure_eight_n_pts=200_nbhd_size=5_axes_ratios=[3 1]__20240702_123409.json'
    r_plot = 0.2
    plot_from_file(filename, r_plot=r_plot)

