import numpy as np
import os
import sys
from typing import Optional

sys.path.append(os.path.abspath('.'))

from ellipsoids.data_handling import get_turkevs_dataset_id
from ellipsoids.data_handling import read_from_json

from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import Dataset
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import Experiment



def find_turkevs_dataset_path(id: str, folder: str = os.path.join("data/turkevs")):
    folder = os.path.join(folder, f"turkevs_id={id}")
    matches = []
    for dirpath, _, filenames in os.walk(folder):
        for filename in filenames:
            if "datasets" in filename and id in filename:
                full_path = os.path.join(dirpath, filename)
                matches.append(full_path)
    if not matches:
        raise FileNotFoundError(f"No dataset found with id={id} within the folder {folder}.")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple datasets found with id={id} within the folder {folder}: {len(matches)} matches found.")
    return matches[0]



def calculate_turkevs(datasets_path: str,
                      complexes_to_calculate: list[Parameters],
                      save_folder: Optional[str]=None):

    all_dicts = read_from_json(datasets_path)
    all_datasets = [Dataset.from_dict(d) for d in all_dicts]

    id = get_turkevs_dataset_id(datasets_path)
    sample_info = all_datasets[0].additional_info
    dataset_id = sample_info.dataset_id if sample_info else "unknown"

    if id != dataset_id:
        raise ValueError(f"Error: there is a mismatch between the dataset ID {dataset_id} \
            and the path-derived ID {id}.")

    # alternative:
    # if id != dataset_id:
    # print(f"Warning: ID mismatch! Using dataset ID '{dataset_id}' instead of path-derived ID '{id}'")
    # id = dataset_id

    for dataset in all_datasets:
        turkevs_dataset_info: TurkevsDatasetInfo = dataset.additional_info
        assert turkevs_dataset_info is not None, "Missing TurkevsDatasetInfo in a dataset!"

        print(f"Performing calculations for transformation: {turkevs_dataset_info.transformation.fullname}, mesh index: {turkevs_dataset_info.mesh_index}.")

        for params in complexes_to_calculate:
            experiment = Experiment(dataset, params)
            experiment.run()
            experiment.save_summary(folder=save_folder)
            experiment.print_execution_time()



if __name__ == '__main__':

    id = "0000"
    path = find_turkevs_dataset_path(id=id)

    # complexes_to_calculate = [EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS),
    #                           EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA),
    #                           Parameters(complex_subtype=ComplexSubtype.RIPS),
    #                           Parameters(complex_subtype=ComplexSubtype.ALPHA)]

    complexes_to_calculate = [EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS),
                              Parameters(complex_subtype=ComplexSubtype.RIPS)]

    for parameters in complexes_to_calculate:
        parameters.save_simplex_tree = False
        if isinstance(parameters, EllipsoidParameters):
            parameters.nbhd_size = 5
            parameters.axes_ratios = np.array([3,1])
            parameters.save_ellipsoid_list = False

    calculate_turkevs(datasets_path=path, complexes_to_calculate=complexes_to_calculate)
