import gudhi as gd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from datetime import datetime

from scipy.spatial import Delaunay

# from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.data_handling import read_from_json
from ellipsoids.data_handling import ensure_folder_exists

# from ellipsoids.topological_computations import reduce_barcode
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
from ellipsoids.common import ComplexType
from ellipsoids.common import ComplexSubtype
from typing import Optional
from dataclasses import dataclass
from ellipsoids.common import ConversionType
from ellipsoids.common import convert


def plot_ellipse(ellipse: Ellipsoid, color='grey', r:float=1, axes=None, r_spherisize:float=np.inf):
    sampleRate = 100
    t = np.linspace(0, 2*np.pi, sampleRate)
    spherisized_axes = spherisize_axes(ellipse.axes_lengths, r=r, r_spherisize=r_spherisize)
    xTemp = r*spherisized_axes[0]*np.cos(t)
    yTemp = r*spherisized_axes[1]*np.sin(t)
    x = ellipse.center[0] + ellipse.axes[0,0]*xTemp + ellipse.axes[1,0]*yTemp
    y = ellipse.center[1] + ellipse.axes[0,1]*xTemp + ellipse.axes[1,1]*yTemp

    plot_context = axes if axes is not None else plt
    plot_context.plot(x,y,c=color, alpha=0.5)



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

    plot_context = axes if axes is not None else plt
    plot_context.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)



def plot_ellipses(ellipse_list: list[Ellipsoid], r: float, axes=None, r_spherisize:float=np.inf):
    for ellipse in ellipse_list:
        plot_ellipse(ellipse, r=r, axes=axes, r_spherisize=r_spherisize)



def plot_ellipsoids(ellipsoid_list, r, axes=None):
    for ellipsoid in ellipsoid_list:
        plot_ellipsoid(ellipsoid, r = r, axes=axes)



def plot_circle(point, r=1, color='grey', axes=None):
    sample_rate = 100
    t = np.linspace(0, 2*np.pi, sample_rate)
    x = point[0] + r*np.cos(t)
    y = point[1] + r*np.sin(t)

    plot_context = axes if axes is not None else plt
    plot_context.plot(x,y,c=color, alpha=0.5)



def plot_circles(points, r=1, axes=None):
    for point in points:
        plot_circle(point, r=r, axes=axes)



def plot_simplex_tree(
        points: np.ndarray,
        simplexTree: gd.SimplexTree,
        filtration: float,
        axes: plt.Axes,
        simplex_color='r'
        ):
    dim = len(points[0])
    if dim > 3:
        raise ValueError('Error: Attempting to plot simplex tree in dimension higher than 3.')

    filtered_simplex_list = list(simplexTree.get_filtration())
    plot_context = axes if axes is not None else plt

    for simplex, simplex_filtration in filtered_simplex_list:
        if simplex_filtration > filtration: continue

        coordinates = np.transpose(points[simplex])
        match len(simplex):
            case 1:
                plot_context.scatter(*coordinates, c='k', zorder=100)
            case 2:
                plot_context.plot(*coordinates, c=simplex_color)
            case 3:
                if dim == 2:
                    plot_context.fill(*coordinates, c=simplex_color, alpha=0.1)



def plot_data_points(points, axes=None):
    plot_context = axes if axes is not None else plt
    if len(points[0]) == 2:
        plot_context.scatter(points[:,0],points[:,1])
    elif len(points[0]) == 3:
        plot_context.scatter(points[:,0],points[:,1],points[:,2])
        # plot_context.remove()
        # TODO figure this out!




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



# def calculate_barcodes(results_list: list[Results]):
#     list_barcodes = []
#     max_lengths = []

#     for results in results_list:
#         barcode, max_length = reduce_barcode(results.barcode)
#         list_barcodes.append(barcode)
#         max_lengths.append(max_length)

#     return list_barcodes, max(max_lengths) * 1.1



def n_bars_to_dict(n_bars: dict):
    return {"nBarsDim0": n_bars[0], "nBarsDim1": n_bars[1], "nBarsDim2": n_bars[2]}



def quick_dim_to_bars(n_bars_list: list[int]) -> dict:
    """
    Given a list of integers (number of bars per dimension),
    returns a list of BarsInDim objects, with dimensions inferred from the index.
    """
    return {i: n for i, n in enumerate(n_bars_list)}



def reduce_barcode_descending(barcode: list[tuple], dim_to_bars: dict):
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
        return f"{dataset.data_type} [n={dataset.n_points}]: " \
            + f"{parameters.complex_type}-{parameters.complex_subtype}; \n"\
            + f"r_s={parameters.r_spherisize}; " \
            + f"axes_ratios={parameters.axes_ratios_to_str()}; " \
            + f"nbhd_size={parameters.nbhd_size}"
    else:
        return f"{dataset.data_type} [n={dataset.n_points}]: {parameters.complex_type}-{parameters.complex_subtype}"



def ax_title_plot(experiment: Experiment):
    plot_parameters = experiment.plot_parameters

    return f"filtration={plot_parameters.filtration}"



def plot_delaunay_triangulation(points, axes):
    plot_context = axes if axes is not None else plt
    triangulation = Delaunay(points)
    plot_context.triplot(points[:,0], points[:,1], triangulation.simplices, color='b', alpha=0.3)


def plot_spatial_data(experiment, ax_plot):
    """Helper function to plot data points, ellipsoids, and simplex tree."""
    plot_parameters = experiment.plot_parameters
    parameters = experiment.parameters
    dataset = experiment.dataset
    results = experiment.results

    if ax_plot is not None:
        ax_plot.set_aspect('equal')
        ax_plot.set_title(ax_title_plot(experiment))

    if plot_parameters.draw_points:
        plot_data_points(dataset.points, axes=ax_plot)
    if plot_parameters.draw_ellipsoids:

        filtration = plot_parameters.filtration
        if parameters.complex_subtype == ComplexSubtype.RIPS:
            radius = convert(filtration, ConversionType.RIPS_TO_RADIUS)
        elif parameters.complex_subtype == ComplexSubtype.ALPHA:
            radius = convert(filtration, ConversionType.ALPHA_TO_RADIUS)
        else:
            radius = filtration

        if isinstance(results, EllipsoidResults):
            plot_ellipses(results.ellipsoid_list,
                          radius,
                          axes=ax_plot,
                          r_spherisize=parameters.r_spherisize)
        else:
            plot_circles(dataset.points, r=radius, axes=ax_plot)
    if plot_parameters.draw_simplex_tree:
        plot_simplex_tree(dataset.points, results.simplex_tree, plot_parameters.filtration, axes=ax_plot)

    if parameters.complex_type == ComplexType.ELLIPSOID \
       and parameters.complex_subtype == ComplexSubtype.ALPHA:
        plot_delaunay_triangulation(dataset.points, ax_plot)



def plot_experiment(experiment: Experiment,
                    ax_barcode: Optional[plt.Axes] = None,
                    ax_plot: Optional[plt.Axes] = None,
                    show=True):

    results = experiment.results
    parameters = experiment.parameters
    plot_parameters = experiment.plot_parameters
    dataset = experiment.dataset

    fig = None

    if ax_barcode == None and ax_plot == None:
        n_axes_per_experiment = 1 + should_plot_spatial_data([experiment])
        fig, axes = plt.subplots(1, n_axes_per_experiment)
        fig.set_size_inches(4*n_axes_per_experiment, 2)

        if n_axes_per_experiment == 1:
            ax_barcode = axes
            ax_plot = None
        else:
            ax_barcode = axes[-1]
            ax_plot = axes[0]

    reduced_barcode = reduce_barcode_descending(results.barcode, plot_parameters.n_bars)
    plot_barcode(reduced_barcode,
                 axes=ax_barcode,
                 infinity=plot_parameters.x_axis_end,
                 axis_start=plot_parameters.x_axis_start)
    ax_barcode.set_title(ax_title_barcode(experiment))
    ax_barcode.set_xlabel("Filtration")
    ax_barcode.set_ylabel(f"(Filtered) barcode \n dim : n_bars = {plot_parameters.n_bars}")


    if should_plot_spatial_data([experiment]):
        plot_spatial_data(experiment, ax_plot=ax_plot)

    if show:
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.6)
        plt.show()

    if fig is None and ax_barcode is not None:
        fig = ax_barcode.get_figure()

    return fig



def should_plot_spatial_data(experiments: list[Experiment]):
    draw_experiment = False
    for experiment in experiments:
        plot_parameters = experiment.plot_parameters
        if plot_parameters.draw_points \
           or plot_parameters.draw_simplex_tree \
           or plot_parameters.draw_ellipsoids:
            draw_experiment = True
            break

    return draw_experiment



def save_figure(figure, path: Optional[str] = None):

    if path == None:
        ensure_folder_exists("data")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = os.path.join("data", f"figure_{timestamp}.png")

    figure.savefig(path, dpi=300)
    print(f"Figure saved to {path}.")



def plot_experiments(experiments: list[Experiment], save_plot: bool = False):

    n_experiments = len(experiments)
    n_axes_per_experiment = 1 + should_plot_spatial_data(experiments)
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
        plot_experiment(experiment, ax_barcode=ax_barcode, ax_plot=ax_draw, show=False)

        print("Done.")


    print("Done.")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.6)
    fig.set_size_inches(4*n_axes_per_experiment, 2*n_experiments)
    plt.show()

    if save_plot:
        save_figure(fig)
