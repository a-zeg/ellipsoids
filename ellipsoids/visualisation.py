import gudhi as gd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from typing import Type

from datetime import datetime

# from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.data_handling import read_variables
from ellipsoids.topological_computations import reduce_barcode
# from ellipsoids.visualisation.barcodePlotting import plot_persistence_barcode, plot_persistence_density
from gudhi.persistence_graphical_tools import _limit_to_max_intervals, __min_birth_max_death
from ellipsoids.topological_computations import spherisize_axes
from ellipsoids.common import Ellipsoid
from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import PlotParameters
from ellipsoids.common import EllipsoidResults
from ellipsoids.common import Results
from ellipsoids.common import Dataset
from ellipsoids.common import Experiment
from typing import Optional
from dataclasses import dataclass


def plot_ellipse(ellipse: Ellipsoid, color='grey', r:float=1, axes=None, r_spherisize:float=np.inf):
    sampleRate = 100
    t = np.linspace(0, 2*np.pi, sampleRate)
    spherisized_axes = spherisize_axes(ellipse.axes_lengths, r=r, r_spherisize=r_spherisize)
    xTemp = r*spherisized_axes[0]*np.cos(t)
    yTemp = r*spherisized_axes[1]*np.sin(t)
    x = ellipse.center[0] + ellipse.axes[0,0]*xTemp + ellipse.axes[1,0]*yTemp
    y = ellipse.center[1] + ellipse.axes[0,1]*xTemp + ellipse.axes[1,1]*yTemp
    if axes is None:
        plt.plot(x,y,c=color)
    else:
        axes.plot(x,y,c=color)



def plot_ellipsoid(ellipsoid: Ellipsoid, color='grey', r:float=1, axes=None):
    # see https://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
    sampleRate = 100

    rx = r * ellipsoid.axes_lengths[0]
    ry = r * ellipsoid.axes_lengths[1]
    rz = r * ellipsoid.axes_lengths[2]
    
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, sampleRate)
    v = np.linspace(0, np.pi, sampleRate)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    all = np.concatenate((np.reshape(x, [-1,1]), np.reshape(y, [-1,1]), np.reshape(z, [-1,1])), axis=1)
    allTransformed = ellipsoid.axes @ np.transpose(all)
    
    x = np.reshape(allTransformed[0,:],(100,100)) + ellipsoid.center[0]
    y = np.reshape(allTransformed[1,:],(100,100)) + ellipsoid.center[1]
    z = np.reshape(allTransformed[2,:],(100,100)) + ellipsoid.center[2]

    if axes is None:
        plt.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)
    else:
        axes.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)



def plot_ellipses(ellipse_list: list[Ellipsoid], r: float, axes=None, r_spherisize:float=np.inf):
    for ellipse in ellipse_list:
        plot_ellipse(ellipse, r=r, axes=axes, r_spherisize=r_spherisize)



def plot_ellipsoids(ellipsoid_list, r, axes=None):
    for ellipsoid in ellipsoid_list:
        plot_ellipsoid(ellipsoid, r = r, axes=axes)



def plotCircle(point, r=1, color='grey', axes=None):
    sample_rate = 100
    t = np.linspace(0, 2*np.pi, sample_rate)
    x = point[0] + r*np.cos(t)
    y = point[1] + r*np.sin(t)
    if axes is None:
        plt.plot(x,y,c=color)
    else:
        axes.plot(x,y,c=color)



def plotCircles(points, r=1, axes=None):
    for point in points:
        plotCircle(point, r=r, axes=axes)



def plot_simplex_tree(points, simplexTree, r, axes):
    dim = len(points[0])
    if dim > 3:
        raise Exception('Error: Attempting to plot simplex tree in dimension higher than 3.')
    generator = simplexTree.get_filtration()
    simplexList = list(generator)

    if axes is None:
        for splx in simplexList:
            if splx[1] <= r:
                vertices = splx[0]
                match len(vertices):
                    case 1:
                        plt.scatter(*np.transpose(points[vertices]), c='k', zorder=100)
                        # points[vertices] gives us points forming the vertices of the simplex
                        # transposing them and taking the * operator returns x-, y-. and z-coords separately
                    case 2:
                        plt.plot(*np.transpose(points[vertices]), c='r')
                    case 3:
                        if dim == 2:
                            plt.fill(*np.transpose(points[vertices]), c='r', alpha=0.1)
    else:
        idx = 1
        for splx in simplexList:
            idx = idx + 1
            
            # if splx[1] <= r:
            if splx[1] <= 2*r: # alt20230927_2: 2r so that it's comparable to Rips
                vertices = splx[0]
                match len(vertices):
                    case 1:
                        axes.scatter(*np.transpose(points[vertices]), c='k', zorder=100)
                    case 2:
                        axes.plot(*np.transpose(points[vertices]), c='r')
                    case 3:
                        if dim == 2:
                            axes.fill(*np.transpose(points[vertices]), c='r', alpha=0.1)



def plot_data_points(points, axes=None):
    if axes is None:
        plt.scatter(points[:,0],points[:,1])
    else:
        axes.scatter(points[:,0],points[:,1])



# from typing import List, Dict, Any
import json

def from_json(filename: str) -> list[Results]:
    with open(filename, 'r') as f:
        data = json.load(f)

    results_list = []

    # Check if any key contains the word "ellipsoid" in its name
    if any("ellipsoid" in key for key in data.keys()):
        # If any "ellipsoid"-related data is present, create an EllipsoidResults object
        ellipsoid_results = EllipsoidResults.from_dict({
            "barcode": data.get("barcode_ellipsoids", []),
            "simplex_tree": data.get("simplex_tree_ellipsoids", []),
            "execution_time": data.get("execution_time"),
            "ellipsoid_list": data.get("ellipsoid_list", [])
        })
        results_list.append(ellipsoid_results)
    else:
        # Otherwise, create a Results object
        results = Results.from_dict({
            "barcode": data.get("barcode_rips", []),
            "simplex_tree": data.get("simplex_tree_rips", []),
            "execution_time": data.get("execution_time")
        })
        results_list.append(results)

    return results_list




def plot_barcode(
    barcode=[], # rename to 'barcode
    alpha=0.6,
    max_intervals=20000,
    inf_delta=0.1,
    colormap=None,
    axes=None,
    fontsize=16,
    axis_start=None,
    infinity=None,
    bar_height=0.8,
):
    '''
    Adapted from GUDHI's plot_persistence_barcode.

    To have a realistic comparison between barcodes for ellipsoids 
    and Rips complexes, it is necessary to have the axes start and
    end at the same values. Since GUDHI's plot_persistence_barcode 
    only provides a way to set a scaling factor and not the absolute 
    start and end values, this adaptaion was created.
    '''


    barcode = _limit_to_max_intervals(
        barcode, max_intervals, key=lambda life_time: life_time[1][1] - life_time[1][0]
    )
    (min_birth, max_death) = __min_birth_max_death(barcode)
    # barcode = sorted(barcode, key=lambda life_time: life_time[1][1] - life_time[1][0])
    barcode = sorted(barcode, key=lambda birth: birth[1][0])

    delta = (max_death - min_birth) * inf_delta
    if infinity is None:
        infinity = max_death + delta
    if axis_start is None:
        axis_start = min_birth - delta
    if axes is None:
        _, axes = plt.subplots(1, 1)
    if colormap is None:
        colormap = plt.cm.Set1.colors

    x = [birth for (dim, (birth, death)) in barcode]
    y = [(death - birth) if death != float("inf") else (infinity - birth) for (dim, (birth, death)) in barcode]
    c = [colormap[dim] for (dim, (birth, death)) in barcode]
 
    axes.barh(range(len(x)), y, left=x, alpha=alpha, color=c, height=bar_height)

    dimensions = {item[0] for item in barcode}
    axes.legend(
        handles=[mpatches.Patch(color=colormap[dim], label=str(dim)) for dim in dimensions],
        loc="best",
    )

    axes.set_title("Persistence barcode", fontsize=fontsize)
    axes.set_yticks([])

    if len(x) != 0:
        axes.set_xlim((axis_start, infinity))

    # -------------------------------- 
    # the next part fixes the scaling 
    margin = 0.4
    n_bars = len(barcode)
    padding = 1.5

    fig = axes.get_figure()
    axes.set_ylim([-padding, (n_bars-1) + padding])
    fig_height = n_bars * (bar_height + margin)# + padding
    axes.invert_yaxis() # temp changing this
    
    fig.set_figheight(fig_height)

    width_inches = 10
    height_inches = n_bars * margin + padding
    fig.set_size_inches(width_inches, height_inches)
    fig.subplots_adjust(top=1)
    # --------------------------------          

    return axes



def calculate_barcodes(results_list: list[Results]):
    list_barcodes = []
    max_lengths = []

    for results in results_list:
        barcode, max_length = reduce_barcode(results.barcode)
        list_barcodes.append(barcode)
        max_lengths.append(max_length)

    return list_barcodes, max(max_lengths) * 1.1



def n_bars_to_dict(n_bars: dict):
    return {"nBarsDim0": n_bars[0], "nBarsDim1": n_bars[1], "nBarsDim2": n_bars[2]}



def quick_dim_to_bars(n_bars_list: list[int]) -> dict:
    """
    Given a list of integers (number of bars per dimension),
    returns a list of BarsInDim objects, with dimensions inferred from the index.
    """
    return {i: n for i, n in enumerate(n_bars_list)}



def reduce_barcode_descending(barcode: list[tuple], dim_to_bars: dict):
    """
    barcode is a barcode
    dim_to_bars is a dictionary with key dim and bars the number of bars in this dimension

    returns: in each dimension n, dim_to_bars[n] longest bars in that dimension
    """
    reduced_barcode = []

    for dim, n_bars in dim_to_bars.items():
        reduced_barcode_dim = [bar for bar in barcode if bar[0] == dim and n_bars>0]
        reduced_barcode_dim.sort(key=lambda bar: bar[1][1]-bar[1][0], reverse=True)
        reduced_barcode.extend(reduced_barcode_dim[:n_bars])

    return reduced_barcode



def find_max_end(barcode, length_tolerance=0.1):

    max = -np.inf

    for bar in barcode:
        bar_end = bar[1][1]
        bar_length = bar_end - bar[1][0]
        if bar_end != np.inf and bar_length > length_tolerance and bar_end > max:
            max = bar_end

    return max



def axis_end_experiments(experiments: list[Experiment]):
    reduced_barcodes \
        = [reduce_barcode_descending(experiment.results.barcode, experiment.plot_parameters.n_bars) \
           for experiment in experiments]
    max_value = max([find_max_end(reduced_barcode) for reduced_barcode in reduced_barcodes])
    return 1.1 * max_value



def ax_title_barcode(experiment: Experiment):
    parameters = experiment.parameters
    dataset = experiment.dataset

    if isinstance(parameters, EllipsoidParameters):
        return f"{dataset.data_type} [n={dataset.n_points()}]: {parameters.complex_type}-{parameters.complex_subtype} r_s={parameters.r_spherisize}"
    else:
        return f"{dataset.data_type} [n={dataset.n_points()}]: {parameters.complex_type}-{parameters.complex_subtype}"



def ax_title_plot(experiment: Experiment):
    plot_parameters = experiment.plot_parameters

    return f"r={plot_parameters.r}"



def plot_experiment(experiment: Experiment,
         ax_barcode: Optional[plt.Axes] = None,
         ax_plot: Optional[plt.Axes] = None):

    results = experiment.results
    parameters = experiment.parameters
    plot_parameters = experiment.plot_parameters
    dataset = experiment.dataset

    if ax_barcode == None and ax_plot == None:
        fig, [ax_barcode, ax_plot] = plt.subplots(1, 2, figsize=(15, 7))

    # Plot the barcode
    reduced_barcode = reduce_barcode_descending(results.barcode, plot_parameters.n_bars)
    plot_barcode(reduced_barcode,
                 axes=ax_barcode,
                 infinity=plot_parameters.x_axis_end,
                 axis_start=plot_parameters.x_axis_start)
    ax_barcode.set_title(ax_title_barcode(experiment))
        # f"Barcode of {parameters.complex_type} - {parameters.complex_subtype}")

    # maybe plot the points
    if plot_parameters.draw_points \
       or plot_parameters.draw_ellipsoids \
       or plot_parameters.draw_simplex_tree:

        plot_data_points(dataset.points, axes=ax_plot)
        ax_plot.set_aspect('equal')
        ax_plot.set_title(ax_title_plot(experiment))

        if plot_parameters.draw_ellipsoids and isinstance(results, EllipsoidResults):
            # scale axes appropriately, depending on r_spherisize
            plot_ellipses(results.ellipsoid_list, plot_parameters.r, axes=ax_plot, r_spherisize=experiment.parameters.r_spherisize)

        if plot_parameters.draw_simplex_tree:
            plot_simplex_tree(dataset.points, results.simplex_tree, plot_parameters.r, axes=ax_plot)




def should_draw_experiments(experiments: list[Experiment]):
    draw_experiment = False
    for experiment in experiments:
        plot_parameters = experiment.plot_parameters
        if plot_parameters.draw_points \
           or plot_parameters.draw_simplex_tree \
           or plot_parameters.draw_ellipsoids:
            draw_experiment = True
            break

    return draw_experiment



def plot_experiments(experiments: list[Experiment]):

    n_experiments = len(experiments)
    n_axes_per_experiment = 1 + should_draw_experiments(experiments)
    fig, axes = plt.subplots(n_experiments,
                             n_axes_per_experiment)
                             # figsize=(n_axes_per_experiment*5, n_experiments * 5))

    x_axis_end = axis_end_experiments(experiments)

    for i, experiment in enumerate(experiments):
        print(f"Generating plot for experiment {i} of {n_experiments}... ", end='', flush=True)

        if n_axes_per_experiment == 1:
            ax_barcode = axes[i]
            ax_draw = None
        else:
            if n_experiments == 1:
                ax_barcode = axes[-1]
                ax_draw = axes[0]
            else:
                ax_barcode = axes[i][-1]
                ax_draw = axes[i][0]

        experiment.plot_parameters.x_axis_end = x_axis_end
        plot_experiment(experiment, ax_barcode=ax_barcode, ax_plot=ax_draw)

        print("Done.")

    print("Done.")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.6)
    fig.set_size_inches(4*n_axes_per_experiment, 2*n_experiments)
    plt.show()



def plot_results(results_list: list[Results],
                 params_list: list[Parameters],
                 dataset: Dataset,
                 r: float = 1,
                 draw_ellipsoids: bool = False,
                 draw_simplex_tree: bool = False):

    print("Plotting results...")

    num_datasets = len(results_list)
    num_plots_per_dataset = 1 + (draw_ellipsoids or draw_simplex_tree)
    fig, axes = plt.subplots(num_datasets,
                             num_plots_per_dataset,
                             figsize=(15, num_datasets * 5))  # 3 subplots per dataset

    list_barcodes, max_length = calculate_barcodes(results_list)

    for i, results in enumerate(results_list):
        print(f"Generating plot for dataset {i} of {num_datasets}... ", end='', flush=True)

        # Plot the barcode
        if num_plots_per_dataset == 1:
            ax = axes[i]
        else:
            ax = axes[i][-1]
        plot_barcode(list_barcodes[i], axes=ax, infinity=max_length, axis_start=-0.1)
        ax.set_title(f"Barcode of {params_list[i].complex_type} - {params_list[i].complex_subtype}")

        # maybe plot the points
        if num_plots_per_dataset > 1:
            plot_data_points(dataset.points, axes=axes[i,0])
            axes[i][0].set_aspect('equal')

            print(type(results))

            if draw_ellipsoids and isinstance(results, EllipsoidResults):
                plot_ellipses(results.ellipsoid_list, r, axes=axes[i][0])

            if draw_simplex_tree:
                plot_simplex_tree(dataset.points, results.simplex_tree, r, axes=axes[i,0])

        print("Done.")

    print("Done.")

    plt.tight_layout()
    plt.show()


