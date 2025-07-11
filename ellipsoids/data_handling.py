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



def sample_from_figure_eight(n: int, a: float = 1, b: float = 0.5, variation: float = 0):
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
    start_T = -np.pi
    end_T = np.pi - (2*np.pi / n) # correcting for periodicity (don't want two points at np.pi (i.e. -np.pi))
    T = np.linspace(start_T, end_T, num=n)

    X = a * np.sin(T)
    Y = a * np.sin(T)**2 * np.cos(T) + b * np.cos(T)

    X = np.column_stack((X, Y))
    X += np.random.default_rng().uniform(0, variation, size=(n, 2))
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
    dataset_summary_dict = cast(TurkevsDatasetInfo, dataset_summary.additional_info).to_dict()
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

