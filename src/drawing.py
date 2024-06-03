
from main import set_filename_parameters
from main import generate_filename
from main import calculate_and_save_ellipsoids_and_rips_data
import data_handling
import numpy as np
import os
# from visualisation import visualisation
from scipy import spatial
from topological_computations import fitEllipsoids
from topological_computations import generateEllipsoidSimplexTree2
import gudhi as gd
import matplotlib.pyplot as plt
import visualisation


def draw_figure(r_plot):

    # Import points
    data_type = 'circle'
    n_pts = 10
    nbhd_size = 3
    axes_ratios = np.asarray([2,1])
    expansion_dim = 3
    points, data_type_params = data_handling.importPoints(data_type=data_type, n_pts = n_pts, variation = 0.1)

    simplexTreeEllipsoids, ellipsoidsList = generateEllipsoidSimplexTree2(points, nbhd_size, axes_ratios)
    # for simplex in simplexTreeEllipsoids.get_simplices():
    #     print(simplex)

    # exit()
    simplexTreeEllipsoidsExpanded = gd.SimplexTree()
    simplexTreeEllipsoidsExpanded = simplexTreeEllipsoids.copy()
    simplexTreeEllipsoidsExpanded.expansion(expansion_dim) # expands the simplicial complex to include 
    barcodeEllipsoids = simplexTreeEllipsoidsExpanded.persistence()
    
    ripsComplex = gd.RipsComplex(points=points)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=expansion_dim)
    barcodeRips = simplexTreeRips.persistence()

    x_axis_end = 0
    for bar in barcodeRips:
        if not(np.isinf(bar[1][1])) and bar[1][1] > x_axis_end:
            x_axis_end = bar[1][1]
    
    for bar in barcodeEllipsoids:
        if not(np.isinf(bar[1][1])) and bar[1][1] > x_axis_end:
            x_axis_end = bar[1][1]

    x_axis_end = x_axis_end * 1.1


    for r in r_plot:
        fig, axes = plt.subplots()
        visualisation.plotEllipses(ellipsoidsList, r, axes)
        visualisation.plotSimplexTree(points, simplexTreeEllipsoidsExpanded, r, axes)
        axes.set_aspect('equal')
        axes.axis('equal')
        axes.set_xlim(xmin=-2.5, xmax=2.5)
        axes.set_ylim(ymin=-2, ymax=2)
        axes.set_axis_off()
        filename = f'{r=}' + '.png'
        fig.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0)
        print('Plot saved to file.')
        fig.clf()




    



    


if __name__ == '__main__':

    r_plot = [0.0, 0.2, 0.32, 1.2]
    draw_figure(r_plot)