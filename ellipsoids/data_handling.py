import os
import argparse
from os import listdir
from os.path import isfile, join
from typing import Any, Type, TypeVar, cast
import sys
import json
from datetime import datetime
import logging
import hashlib

import gudhi as gd
import numpy as np
from scipy.io import loadmat

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Ellipsoid, TurkevsDatasetInfo, ExperimentSummary


logger = logging.getLogger(__name__)



def sample_from_circle(n_pts: int = 100, variation: float = 0.1, outlier: bool = False):
    if outlier is True: 
        n_pts = n_pts - 1

    r = 1
    t = np.linspace(0, 2*np.pi * (n_pts-1)/n_pts, n_pts)
    x = r*np.cos(t) + variation * np.random.rand(n_pts)
    y = r*np.sin(t) + variation * np.random.rand(n_pts)
    output = np.vstack((x,y)).transpose()

    if outlier is True:
        output = np.append(output,[[0,0]],axis=0)
    
    return output



def sample_from_ellipse(n_pts: int = 100, a: float = 2, b: float = 1,
                        variation: float = 0.1, shift: float = 0):
    t = np.linspace(0, 2*np.pi * (n_pts-1)/n_pts, n_pts) + shift
    x = a * np.cos(t) + variation * np.random.rand(n_pts)
    y = b * np.sin(t) + variation * np.random.rand(n_pts)
    return np.vstack((x,y)).transpose()



def sample_from_cassini_oval(n_pts: int = 100, variation: float = 0.1):
    t = np.linspace(-1, 1, int(n_pts/2))
    #t = np.sign(t)*np.abs(t)**(1/4)
    x = np.concatenate((t,t)) + variation * np.random.rand(n_pts)
    yh = (t**2 + 0.5) * np.sqrt(1 - t**2)
    y = np.concatenate((-yh, yh)) + variation * np.random.rand(n_pts)
    
    return np.vstack((x,y)).transpose()



def sample_from_sphere(n_pts: int = 100, ambient_dim: int = 3, r: float = 1):
    vec = np.random.randn(ambient_dim, n_pts)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.transpose()



def sample_from_torus(n_pts: int = 100, R: float = 2, r: float = 1):
    nPtsSampled = 0
    theta = np.zeros([n_pts])
    phi = np.zeros([n_pts])

    # rejection sampling
    while nPtsSampled < n_pts:
        thetaSample = 2 * np.pi * np.random.rand()
        phiSample = 2 * np.pi * np.random.rand()
        W = 2 * np.pi * np.random.rand()

        if W <= (R + r * np.cos(thetaSample))/(R+r):
            theta[nPtsSampled] = thetaSample
            phi[nPtsSampled] = phiSample
            nPtsSampled += 1


        x = (R + r * np.cos(theta)) * np.cos(phi)
        y = (R + r * np.cos(theta)) * np.sin(phi)
        z = r * np.sin(theta)

    return np.vstack((x,y,z)).transpose() 



def sample_from_figure_eight(n: int, a: float = 1, b: float = 0.5, variation: float = 0, seed: int = None):
    # adapted from Bastian Rieck
    """Sample a set of points from a figure eight curve.

    Parameters
    ----------
    n : int
        Number of points to sample

    a : float
        Controls extents of the curve. A larger `a` parameter will
        result in larger scaling.

    b : float
        Controls neck size of the curve. A larger `b` parameter will
        result in an increased neck size.

    Returns
    -------
    np.array
        Array of shape (n, 2). Will contain the sampled points.
    """
    rng = np.random.default_rng(seed)

    start_T = -np.pi
    end_T = np.pi - (2*np.pi / n) # correcting for periodicity (don't want two points at np.pi (i.e. -np.pi))
    T = np.linspace(start_T, end_T, num=n)

    X = a * np.sin(T)
    Y = a * np.sin(T)**2 * np.cos(T) + b * np.cos(T)

    X = np.column_stack((X, Y))
    X += rng.uniform(0, variation, size=(n, 2))
    return X



def sample_from_annulus(n: int, r: float = 1, R: float = 2, seed=None):
    # taken from Bastian Rieck
    """Sample points from a 2D annulus.

    This function samples `N` points from an annulus with inner radius `r`
    and outer radius `R`.

    Parameters
    ----------
    n : int
        Number of points to sample

    r : float
        Inner radius of annulus

    R : float
        Outer radius of annulus

    seed : int, instance of `np.random.Generator`, or `None`
        Seed for the random number generator, or an instance of such
        a generator. If set to `None`, the default random number
        generator will be used.

    Returns
    -------
    torch.tensor of shape `(n, 2)`
        Tensor containing sampled coordinates.
    """
    if r >= R:
        raise RuntimeError(
            'Inner radius must be less than or equal to outer radius'
        )

    rng = np.random.default_rng(seed)
    thetas = rng.uniform(0, 2 * np.pi, n)

    # Need to sample based on squared radii to account for density
    # differences.
    radii = np.sqrt(rng.uniform(r ** 2, R ** 2, n))

    X = np.column_stack((radii * np.cos(thetas), radii * np.sin(thetas)))
    return X



def sample_spiral(n, turns=3, inner_radius=0.1, outer_radius=1.0):
    """
    Open 1D spiral in R^2.
    No holes.
    """
    t = np.random.uniform(0, 1, size=n)

    theta = 2 * np.pi * turns * t
    r = inner_radius + (outer_radius - inner_radius) * t

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    X = np.array([x, y]).T

    return X


def sample_3d_bottleneck(n, gap=0.12):
    """
    Closed non-planar 1D curve in R^3.

    In xy projection it resembles a figure eight, but the two
    branches near the crossing are separated in z.

    Smaller gap -> stronger bottleneck / near-self-intersection.

    Topologically this is a circle (1 hole).
    """
    t = np.random.uniform(0, 2 * np.pi, size=n)

    x = np.sin(t)
    y = np.sin(2 * t)
    z = gap * np.cos(t)

    X = np.array([x, y, z]).T

    return X




def sample_braided_theta(
        n,
        holes=2,
        radius=0.25,
        twists=1.5,
):
    """
    Non-planar generalized theta graph in R^3.

    holes + 1 strands connect the same two endpoints, so

        beta_1 = holes.

    The strands twist around one another without intersecting.

    Smaller radius makes distinct strands closer in Euclidean
    distance.
    """

    number_of_strands = holes + 1

    curves = []

    for j in range(number_of_strands):

        phase = (
            2 * np.pi
            * j
            / number_of_strands
        )

        def curve(t, phase=phase):

            x = 2 * t - 1

            # Vanishes at endpoints so all strands meet there.
            r = radius * np.sin(np.pi * t)

            angle = (
                phase
                + 2 * np.pi * twists * t
            )

            y = r * np.cos(angle)
            z = r * np.sin(angle)

            return np.array([
                x,
                y,
                z
            ]).T

        curves.append(curve)

    return _sample_curve_family(
        curves,
        n
    )


# ============================================================
# 2. Folded ladder
# ============================================================

def sample_folded_ladder(
        n,
        holes=4,
        width=0.25,
        bend=1.8 * np.pi,
        radius=1.0,
):
    """
    Ladder graph with the requested number of holes, bent into
    a strong U / horseshoe shape in R^3.

    Topologically:
        beta_1 = holes.

    As bend approaches 2*pi, distant parts of the ladder become
    close in Euclidean space even though they are far apart along
    the graph.
    """

    curves = []

    def embed(u, v):
        """
        Bend the original (u,v) ladder into R^3.
        u ranges from 0 to holes.
        """

        q = u / holes

        angle = (
            -bend / 2
            + bend * q
        )

        x = radius * np.sin(angle)
        y = v
        z = radius * np.cos(angle)

        return np.array([
            x,
            y,
            z
        ]).T

    # Two long rails
    for side in [-1, 1]:

        v = side * width

        for i in range(holes):

            def rail(t, i=i, v=v):

                u = i + t

                return embed(
                    u,
                    np.full_like(t, v)
                )

            curves.append(rail)

    # Rungs
    for i in range(holes + 1):

        def rung(t, i=i):

            u = np.full_like(t, i)

            v = (
                -width
                + 2 * width * t
            )

            return embed(u, v)

        curves.append(rung)

    return _sample_curve_family(
        curves,
        n
    )


# ============================================================
# 3. Multiscale necklace
# ============================================================
#



import numpy as np


def _sample_curve_family(
        curves,
        n,
        weights=None,
        sampling_randomness=0.5,
        resolution=3000,
):
    """
    Sample from a collection of parametric curves.

    sampling_randomness:
        0.0 = perfectly regular spacing along each curve
        1.0 = maximally jittered within regular spacing cells

    Even at 1.0 this is more uniform than iid random sampling.
    """

    curve_data = []
    lengths = []

    for curve in curves:

        t_dense = np.linspace(0, 1, resolution)
        X_dense = curve(t_dense)

        ds = np.linalg.norm(
            np.diff(X_dense, axis=0),
            axis=1
        )

        s_dense = np.concatenate([
            [0],
            np.cumsum(ds)
        ])

        curve_data.append(
            (curve, t_dense, s_dense)
        )

        lengths.append(
            s_dense[-1]
        )

    lengths = np.asarray(lengths)

    if weights is None:
        probabilities = lengths / lengths.sum()

    else:
        weights = np.asarray(
            weights,
            dtype=float
        )

        probabilities = (
            weights / weights.sum()
        )

    # --------------------------------------------------------
    # Allocate approximately the correct number of points
    # to each curve, but without multinomial fluctuations
    # --------------------------------------------------------

    expected = n * probabilities

    counts = np.floor(
        expected
    ).astype(int)

    remainder = n - counts.sum()

    if remainder > 0:

        fractional = (
            expected - counts
        )

        order = np.argsort(
            fractional
        )[::-1]

        counts[
            order[:remainder]
        ] += 1

    # --------------------------------------------------------
    # Sample each curve
    # --------------------------------------------------------

    dimension = curves[0](
        np.array([0.5])
    ).shape[1]

    X_parts = []

    for count, (
        curve,
        t_dense,
        s_dense
    ) in zip(counts, curve_data):

        if count == 0:
            continue

        L = s_dense[-1]

        # Regular positions centred inside equal-length cells
        regular = (
            np.arange(count) + 0.5
        ) / count

        # Random displacement of at most half a cell
        jitter = np.random.uniform(
            -0.5,
            0.5,
            size=count
        ) / count

        u = (
            regular
            + sampling_randomness * jitter
        )

        u = np.clip(
            u,
            0,
            1
        )

        s = u * L

        t = np.interp(
            s,
            s_dense,
            t_dense
        )

        X_parts.append(
            curve(t)
        )

    X = np.concatenate(
        X_parts,
        axis=0
    )

    # Shuffle so points aren't grouped by curve
    np.random.shuffle(X)

    return X


def sample_multiscale_loops(
        n,
        holes=4,
        largest_radius=0.7,
        radius_ratio=0.65,
        stem_length=0.25,
        spacing=1.8,
        loop_mass=0.78,
        stem_mass=0.14,
        backbone_mass=0.08,
        loop_length_power=0.8,
        sampling_randomness=0.35,
):
    """
    Connected collection of loops of progressively different size.

    beta_1 = holes

    Improvements over the previous version:

    - larger loops receive more points than smaller loops
    - sampling along each piece is more spatially uniform
    - the backbone still receives relatively few points

    Parameters
    ----------
    loop_length_power:

        0 -> every loop gets the same number of points

        1 -> points assigned fully proportional to circumference

        values around 0.7--0.9 are usually a good compromise

    sampling_randomness:

        0 -> almost perfectly equally spaced

        1 -> strongly jittered, but still less clumpy than iid sampling
    """

    if holes < 1:
        raise ValueError(
            "holes must be >= 1"
        )

    radii = np.array([
        largest_radius
        * radius_ratio**i
        for i in range(holes)
    ])

    x_positions = (
        np.arange(holes)
        * spacing
    )

    backbone_curves = []
    stem_curves = []
    loop_curves = []

    # ========================================================
    # Backbone
    # ========================================================

    if holes > 1:

        for i in range(
            holes - 1
        ):

            x0 = x_positions[i]
            x1 = x_positions[i + 1]

            def backbone(
                t,
                x0=x0,
                x1=x1,
            ):

                x = (
                    x0
                    + (x1 - x0) * t
                )

                y = np.zeros_like(t)

                return np.array([
                    x,
                    y
                ]).T

            backbone_curves.append(
                backbone
            )

    # ========================================================
    # Stems and loops
    # ========================================================

    for x0, r in zip(
        x_positions,
        radii
    ):

        def stem(
            t,
            x0=x0
        ):

            x = np.full_like(
                t,
                x0
            )

            y = (
                stem_length * t
            )

            return np.array([
                x,
                y
            ]).T

        stem_curves.append(
            stem
        )

        def loop(
            t,
            x0=x0,
            r=r
        ):

            theta = (
                2 * np.pi * t
            )

            x = (
                x0
                + r * np.sin(theta)
            )

            y = (
                stem_length
                + r
                - r * np.cos(theta)
            )

            return np.array([
                x,
                y
            ]).T

        loop_curves.append(
            loop
        )

    # ========================================================
    # Weights
    # ========================================================

    curves = (
        backbone_curves
        + stem_curves
        + loop_curves
    )

    weights = []

    # Backbone
    if len(backbone_curves) > 0:

        weights += [
            backbone_mass
            / len(backbone_curves)
        ] * len(backbone_curves)

    # Stems
    weights += [
        stem_mass
        / len(stem_curves)
    ] * len(stem_curves)

    # --------------------------------------------------------
    # Loops:
    #
    # larger radius -> larger weight
    #
    # power = 0 gives equal loop weights
    # power = 1 gives circumference-proportional weights
    # --------------------------------------------------------

    loop_weights = (
        radii ** loop_length_power
    )

    loop_weights /= (
        loop_weights.sum()
    )

    loop_weights *= loop_mass

    weights += list(
        loop_weights
    )

    # ========================================================
    # Sample
    # ========================================================

    X = _sample_curve_family(
        curves,
        n,
        weights=weights,
        sampling_randomness=sampling_randomness,
    )

    X -= X.mean(
        axis=0
    )

    return X





















def import_maxmin_mat(path: str):
    return np.asarray(loadmat(path)['maxmin_points']) # size: 6040 x 24 



def get_timestamp():
    return datetime.now().strftime("_%Y%m%d_%H%M%S")



class CustomEncoder(json.JSONEncoder):
    def default(self, obj):

        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
        elif isinstance(obj, Ellipsoid):
            obj_data = {
                "center": obj.center,
                "axes": obj.axes.tolist(),
                "axes_lengths": obj.axes_lengths.tolist()
            }
            return obj_data
        
        elif isinstance(obj,gd.SimplexTree):
            return list(obj.get_filtration())
        
        elif isinstance(obj, np.integer):
            return int(obj)
        
        return json.JSONEncoder.default(self, obj)



def ensure_folder_exists(filename):
    """
    Checks if the folder for the given filename exists,
    and creates it if it does not.
    """
    folder = os.path.dirname(filename)
    if folder and not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


def save_to_json(
        data: Any,
        filename: str = os.path.join("data","test"),
        add_timestamp: bool = True):
    
    logger.info("Saving data...")

    if not filename.endswith(".json"):
        filename = filename + ".json"

    if add_timestamp:
        filename = filename.replace(".json", f"_{get_timestamp()}.json")

    ensure_folder_exists(filename);

    json_string = json.dumps(data, cls=CustomEncoder, indent=2)
    with open(filename, 'w') as outfile:
        outfile.write(json_string)
    logger.info(f"Data saved to {filename}.")

    return filename



def read_from_json(filename):
    with open(filename, "r") as f:
        json_vars = json.load(f)
    return json_vars
    # return json_process_variables(json_vars)



def _json_process_ellipsoid_list(ellipsoid_list_raw):
    ellipsoid_list = []
    for ellipsoid in ellipsoid_list_raw:
        ellipsoid_list.append(
            Ellipsoid(ellipsoid["center"],
                      np.asarray(ellipsoid["axes"]), \
                      np.asarray(ellipsoid["axes_lengths"]))
            )
    return ellipsoid_list



def _json_process_simplex_tree(simplex_tree_raw):
    simplex_tree = gd.SimplexTree()
    for simplex_tree_entry in simplex_tree_raw:
        simplex_tree.insert(simplex_tree_entry[0],simplex_tree_entry[1])
    return simplex_tree



def json_process_variables(json_vars: dict):
    vars = json_vars
    if "ellipsoid_list" in json_vars:
        vars["ellipsoid_list"] = _json_process_ellipsoid_list(json_vars["ellipsoid_list"])

    if "simplex_tree" in json_vars:
        vars["simplex_tree"] = _json_process_simplex_tree(json_vars["simplex_tree"])
    return vars



def filter_barcode(barcode, dim):
    filtered_barcode = []
    for bar in barcode:
        if bar[0] == dim:
            filtered_barcode.append(bar[1])
    return filtered_barcode



def print_list_of_simplices(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    for splx in simplexList:
        print(splx)



def get_paths_of_files_in_a_folder(folder: str, extension=None):
    if extension is not None:
        filenames = [f for f in listdir(folder) if isfile(join(folder, f)) if f.endswith(extension)]
    else:
        filenames = [f for f in listdir(folder) if isfile(join(folder, f))]

    paths = [os.path.join(folder, f) for f in filenames] 
    return paths



def remove_dim_from_barcode(barcode):
    result_barcode = []
    for bar in barcode:
            result_barcode.append(bar[1])
    return result_barcode



# T = TypeVar("T")

# def check_type(variable: Any, expected_type: Type[T]) -> T:
#     if isinstance(variable, expected_type):
#         return cast(T, variable) # signaling to the type checker (it does nothing at runtime).
#     raise TypeError(
#         f"Expected the variable {variable.__name__} to be of type {expected_type.__name__}, "
#         f"got {type(variable).__name__} instead."
#     )



def generate_signature_from_dataset_summary_parameters(dataset_summary, parameters) -> str:
    # dataset_summary_dict = cast(TurkevsDatasetInfo, dataset_summary.additional_info).to_dict()
    dataset_summary_dict = dataset_summary.additional_info.to_dict()
    dataset_summary_json = json.dumps(dataset_summary_dict, sort_keys=True, cls=CustomEncoder)
    parameters_dict = parameters.to_dict()
    parameters_json = json.dumps(parameters_dict, sort_keys=True, cls=CustomEncoder)
    full_str = f"{dataset_summary_json}_{parameters_json}"
    return hashlib.sha256(full_str.encode()).hexdigest()



def build_signature_index_from_jsonl(summaries_jsonl_path: str) -> dict:
    from ellipsoids.turkevs.turkevs_utils import iter_experiment_summaries_from_compressed_jsonl

    summary_index = {}

    for summary in iter_experiment_summaries_from_compressed_jsonl(summaries_jsonl_path):
        sig = generate_signature_from_dataset_summary_parameters(
            summary.dataset_summary,
            summary.parameters
        )
        summary_index[sig] = {
            "parameters": summary.parameters.to_dict(),
            "dataset_summary": summary.dataset_summary.to_dict()
        }
    return summary_index



def parse_args():
    parser = argparse.ArgumentParser(description="Calculate ellipsoid barcodes on a dataset folder.")
    parser.add_argument(
        "--folder",
        type=str,
        default=None,
        help="Path to the folder containing the datasets and where the results will be stored."
        )
    parser.add_argument(
        "--datasets_start",
        type=int,
        default=None,
        help="Start index for dataset slicing."
    )
    parser.add_argument(
        "--datasets_end",
        type=int,
        default=None,
        help="End index for dataset slicing."
    )
    return parser.parse_args()

