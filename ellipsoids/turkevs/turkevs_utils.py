#!/usr/bin/env python3

import os
import sys
import re
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional, Any, cast
from collections import defaultdict
import logging
import hashlib
import gzip
import json

sys.path.append(os.path.abspath('.'))

import matplotlib.pyplot as plt
import numpy as np
from ellipsoids.common import TurkevsTransformation
from ellipsoids.common import Dataset
from ellipsoids.common import Parameters
from ellipsoids.common import ExperimentSummary
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import TurkevsTransformation
from ellipsoids.common import Parameters
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import ExperimentSummary
from ellipsoids.common import Dataset
from ellipsoids.data_handling import CustomEncoder
from ellipsoids.data_handling import read_from_json
from ellipsoids.data_handling import get_paths_of_files_in_a_folder
from ellipsoids.data_handling import check_type
from ellipsoids.data_handling import filter_barcode
from ellipsoids.data_handling import filter_barcode
from ellipsoids.data_handling import get_paths_of_files_in_a_folder
from ellipsoids.data_handling import read_from_json
import ellipsoids.turkevs.model as model
from ellipsoids.turkevs.data_construction import build_dataset_holes
from ellipsoids.turkevs.data_construction import calculate_point_clouds_under_trnsf
from ellipsoids.turkevs.model import get_score
import ellipsoids.turkevs.data_construction as data_construction
from ellipsoids.turkevs.ph import extend_pds_to_length
import ellipsoids.turkevs.ph as ph
import ellipsoids.turkevs.model as model
import ellipsoids.turkevs.ml as ml
import ellipsoids.turkevs.ph_ml as ph_ml
import ellipsoids.turkevs.nn_shallow as nn_shallow
import ellipsoids.turkevs.nn_deep as nn_deep
import ellipsoids.turkevs.point_net as point_net
import ellipsoids.turkevs.plots as plots



logger = logging.getLogger(__name__)



def generate_turkevs_results_folder_name(n_point_clouds: int, n_points: int, seed: int):
    return f"npc={n_point_clouds}_np={n_points}_s={seed}"



def generate_turkevs_datasets(n_point_clouds: int,
                              n_points: int,
                              seed: int,
                              dataset_id=""):
    np.random.seed(seed)
    logger.info("Constructing Turkevs datasets...")

    initial_point_clouds, labels, _ = build_dataset_holes(n_point_clouds, n_points)
    all_datasets: list[Dataset] = []

    if n_point_clouds != len(initial_point_clouds):
        raise ValueError("Mismatch in the number of point clouds.")

    for transformation in TurkevsTransformation:
        logger.info(f"Applying transformation: {transformation.fullname}.")
        transformed_point_clouds = (
            initial_point_clouds if transformation is TurkevsTransformation.STANDARD
            else calculate_point_clouds_under_trnsf(initial_point_clouds,
                                                    transformation=transformation.fullname)
        )
        for mesh_index, (points,label) in enumerate(zip(transformed_point_clouds, labels)):
            dataset = Dataset(
                points=points,
                data_type="turkevs",
                additional_info=TurkevsDatasetInfo(
                    dataset_id=dataset_id,
                    seed=seed,
                    point_cloud_index=mesh_index,
                    transformation=transformation,
                    label=label,
                    n_point_clouds=n_point_clouds
                )
            )
            all_datasets.append(dataset)
    logger.info("Datasets constructed.")
    return all_datasets



def add_dataset_id(datasets: list[Dataset], dataset_id: str):
    for dataset in datasets:
        info = cast(TurkevsDatasetInfo, dataset.additional_info)
        info.dataset_id = dataset_id
    return datasets



def generate_datasets_hash(datasets: list[Dataset], hash_length: int = 6):
    hasher = hashlib.sha256()

    for dataset in datasets:
        hasher.update(dataset.points.astype(np.float32).tobytes())

        point_cloud_index = getattr(dataset.additional_info, "point_cloud_index", None)
        if point_cloud_index is not None:
            hasher.update(str(point_cloud_index).encode())

        transformation = getattr(dataset.additional_info, "transformation", None)
        if transformation is not None:
            hasher.update(transformation.name.encode())

        label = getattr(dataset.additional_info, "label", None)
        if label is not None:
            hasher.update(str(label).encode())

    full_hash = hasher.hexdigest()
    return full_hash[:hash_length]



# def format_dataset_number_id(id: int):
#     return str(id).zfill(4)


# def check_id_exists(parent_folder: str, folder: str):
#     paths = os.listdir(parent_folder)
#     if folder in paths:
#         raise FileExistsError(f"A folder {folder} already exists in folder {parent_folder}.")



# def generate_dataset_number_id(folder):
#     '''
#     The datasets used in the Turkevs holes tests are generated using the code from the paper and saved in a json file.
#     To keep track of which datasets correspond to which data / graphs, 'id' is introduced.

#     This function checks the ids of all the dataset files in the given folder and returns the next available one.
#     '''
#     if not os.path.exists(folder):
#         return format_dataset_number_id(0)
#     paths = os.listdir(folder)
#     existing_ids = set()

#     # add all ids to existing_ids
#     for path in paths:
#         match = re.search(r'id=(\d+)', path)
#         if match:
#             try:
#                 current_id = int(match.group(1))
#                 existing_ids.add(current_id)
#             except ValueError:
#                 continue

#     # find a new id
#     new_id = 0
#     while new_id in existing_ids:
#         new_id += 1

#     return format_dataset_number_id(new_id)



def get_turkevs_dataset_id(path: str):
    match = re.search(r'id=([a-f0-9]+)', path)
    if match:
        return str(match.group(1))
    raise ValueError(f"Dataset ID could not be found in the path: {path}")



def read_turkevs_datasets(datasets_path):
    logger.info(f"Reading in Turkevs datasets from {datasets_path}...")
    datasets_dicts = read_from_json(datasets_path)
    datasets = [Dataset.from_dict(d) for d in datasets_dicts]
    id = get_turkevs_dataset_id(datasets_path)
    datasets = ensure_turkevs_datasets_consistency(datasets, id)
    logger.info("Datasets read.")
    return datasets



def ensure_turkevs_datasets_consistency(datasets: list[Dataset], dataset_id: Optional[str] = None):
    for dataset in datasets:
        dataset = check_type(dataset, Dataset)
        dataset.additional_info = check_type(dataset.additional_info, TurkevsDatasetInfo)
        if dataset_id is not None and dataset.additional_info.dataset_id != dataset_id:
            raise ValueError(f"Error: there is a mismatch between the dataset ID {dataset_id} \
                and the provided ID {id}.")
    return datasets



def load_all_experiment_summaries(folder: str) -> list[ExperimentSummary]:
    experiment_summaries = []
    paths = get_paths_of_files_in_a_folder(folder=folder, extension=".json")
    for path in paths:
        summary_dict = read_from_json(path)
        summary = ExperimentSummary.from_dict(summary_dict)
        experiment_summaries.append(summary)
    return experiment_summaries



def find_files_by_keyword(folder, keyword="", extension=".json"):
    paths = get_paths_of_files_in_a_folder(folder=folder, extension=extension)
    keyword_paths = [path for path in paths if keyword in os.path.basename(path)]
    if len(keyword_paths) == 0:
        raise FileNotFoundError(f"No file found for id={id} in folder {folder} and containing keyword: {keyword}.")
    return keyword_paths



def filter_turkevs_dict(data: dict[str, list],
                        indices: np.ndarray,
                        transformations: list[TurkevsTransformation] = [t for t in TurkevsTransformation]):
    '''
    Data is a dictionary with keys 'std', 'trns', etc and values
    all point clouds without transformations, with translation, etc.
    This function only keeps the point clouds at indices given by the variable indices.
    '''
    filtered_data = {}
    for transformation in transformations:
        data_transformation = data[transformation.shortname]
        filtered_data[transformation.shortname] = [data_transformation[i] for i in indices]
    return filtered_data



def divide_indices_into_train_and_test(indices, train_percentage=0.8, test_percentage=0.2):
    n_indices = len(indices)
    train_size = int(train_percentage * n_indices)
    test_size = int(test_percentage * n_indices)
    train_indices = np.random.choice(indices, size=train_size, replace=False)
    non_train_indices = np.setdiff1d(indices, train_indices)
    test_indices = np.random.choice(non_train_indices, size=test_size, replace=False)
    return train_indices, test_indices



@dataclass
class EvaluationResults:
    accuracy: dict
    pipeline_name: str
    parameters: Optional[Parameters] = None

    def to_dict(self):
        result = {
            "accuracy": self.accuracy,
            "pipeline_name": self.pipeline_name,
        }
        if self.parameters is not None:
            result["parameters"] = self.parameters.to_dict()
        return result

    @classmethod
    def from_dict(cls, data: dict):
        parameters = None
        if "parameters" in data and data["parameters"] is not None:
            parameters = Parameters.from_any_dict(data["parameters"])
        return cls(
            accuracy=data["accuracy"],
            pipeline_name=data["pipeline_name"],
            parameters=parameters
        )



def get_datasets_per_transformation_sorted(datasets: list[Dataset]) -> dict[str, list[np.ndarray]]:
    """
    Returns a dictionary with keys ["std", "rot", ...] and values lists of point clouds, sorted by point_cloud_index.
    """
    datasets_per_transformation: dict[str, list[Dataset]] = defaultdict(list)

    for d in datasets:
        if isinstance(d.additional_info, TurkevsDatasetInfo):
            key = d.additional_info.transformation.shortname
            datasets_per_transformation[key].append(d)

    sorted_datasets_per_transformation: dict[str, list[np.ndarray]] = {}

    for t, datasets in datasets_per_transformation.items():
        # validate_point_cloud_indices(datasets)
        sorted_datasets = sorted(datasets,
                                 key=lambda d: cast(TurkevsDatasetInfo, d.additional_info).point_cloud_index)
        sorted_datasets_per_transformation[t] = [d.points for d in sorted_datasets]

    return sorted_datasets_per_transformation



def get_labels_sorted(datasets: list[Dataset]) -> np.ndarray:
    """
    Returns a list of labels sorted by Dataset.additional_info.point_cloud_index.
    The labels are obtained from the non-transformed datasets, i.e. those with
    Dataset.additional_info.transformation = TurkevsTransformation.STANDARD.
    """
    logger.info("Extracting and sorting labels... ")
    standard_datasets = [ d for d in datasets
                          if isinstance(d.additional_info, TurkevsDatasetInfo)
                          and d.additional_info.transformation == TurkevsTransformation.STANDARD ]
    standard_datasets_sorted = sorted(standard_datasets,
                                      key=lambda d: cast(TurkevsDatasetInfo, d.additional_info).point_cloud_index)
    labels = [cast(TurkevsDatasetInfo, d.additional_info).label for d in standard_datasets_sorted]
    logger.info("Labels extracted and sorted.")
    return np.asarray(labels)



def get_barcodes_dim_1_per_parameters_per_transformation_sorted(experiment_summaries: list[ExperimentSummary]) \
        -> dict[Parameters, dict[str, list[np.ndarray]]]:
    """
    Returns a dictionary with:
    - keys: parameters
    - values: dictionary with:
    -- keys: TurkevsTransformation.shortname
    -- values: list of barcodes, sorted by point_cloud_index
    """

    # group by parameters and then by transformations
    grouped: dict[Any, dict[str, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    for summary in experiment_summaries:
        info = summary.dataset_summary.additional_info
        if isinstance(info, TurkevsDatasetInfo):
            t = info.transformation.shortname
            grouped[summary.parameters][t].append((info.point_cloud_index, summary.barcode))

    # sort by point_cloud_index
    sorted_and_filtered: dict[Any, dict[str, list[np.ndarray]]] = {}
    for params, barcodes_by_transformation in grouped.items():
        sorted_and_filtered[params] = {}
        for t, point_cloud_index_and_barcodes in barcodes_by_transformation.items():
            sorted_barcodes = [np.asarray(filter_barcode(barcode=barcode, dim=1))
                               for _, barcode in sorted(point_cloud_index_and_barcodes,
                                                        key=lambda x: x[0])]
            padded_barcodes = pad_to_max_length(sorted_barcodes)
            sorted_and_filtered[params][t] = padded_barcodes

    return sorted_and_filtered



def pad_to_max_length(barcodes: list[np.ndarray]):
    max_length = max([len(barcode) for barcode in barcodes])
    return extend_pds_to_length(barcodes, max_length)



def train_model(ml_module, data_train, labels_train):
    logger.info(f"Training with module {ml_module.__name__}...")
    logger.debug("Tuning the hyperparameters...")
    model_ = ml_module.tune_hyperparameters(data_train, labels_train)
    logger.debug("Hyperparameters tuned.")
    logger.debug("Training the model...")
    model_trained, _ = model.fit(data_train, labels_train, model_)
    logger.info("Training completed.")
    return model_trained



def get_scores_under_trnsfs(data_trnsfs, trnsfs, labels, model):
    accs = {}
    for t in trnsfs:
        data = data_trnsfs[t]
        accs[t] = get_score(data, labels, model)
    return accs



def test_model(model_trained,
               data_test,
               labels_test):
    logger.info("Testing the model by evaluating noise robustness...")
    trnsfs = [transformation.shortname for transformation in TurkevsTransformation]
    accs_trnsfs = get_scores_under_trnsfs(data_test,trnsfs,labels_test,model_trained)
    logger.info("Testing completed.")
    return accs_trnsfs



def evaluate_model(
        data_per_transformations,
        labels,
        train_indices,
        test_indices,
        ml_module,
        pipeline_name
        ) -> EvaluationResults:
    logger.info(f"Running pipeline {pipeline_name}...")
    labels_train = labels[train_indices]
    labels_test = labels[test_indices]
    data_train = [data_per_transformations[TurkevsTransformation.STANDARD.shortname][i] for i in train_indices]
    data_test_per_transformations = filter_turkevs_dict(data_per_transformations, test_indices)

    if isinstance(data_per_transformations["std"], np.ndarray):
        data_train = np.asarray(data_train)
        data_test_per_transformations = {t: np.asarray(data_test_per_transformations[t]) for t in data_test_per_transformations}

    model_trained = train_model(ml_module, data_train, labels_train)
    accs_trnsfs = test_model(model_trained, data_test_per_transformations, labels_test)
    return EvaluationResults(pipeline_name=pipeline_name, accuracy=accs_trnsfs)



def get_train_and_test_indices(datasets_per_transformation: dict["str", list]):
    logger.info("Generating train and test indices...")
    n_std_point_clouds = len(datasets_per_transformation[TurkevsTransformation.STANDARD.shortname])
    indices = np.arange(n_std_point_clouds)
    logger.info("Train and test indices generated.")
    train_indices, test_indices = divide_indices_into_train_and_test(indices)
    return train_indices, test_indices



def point_clouds_from_scratch(n_point_clouds, n_points) -> tuple[dict["str", list], np.ndarray]:
    from ellipsoids.turkevs.data_construction import build_dataset_holes
    data_pc, labels, _ = build_dataset_holes(n_point_clouds,n_points)
    data_pc_trnsfs: dict["str", list] = {t.shortname: data_construction.calculate_point_clouds_under_trnsf(data_pc, t.fullname)
                      for t in TurkevsTransformation if t is not TurkevsTransformation.STANDARD}
    data_pc_trnsfs[TurkevsTransformation.STANDARD.shortname] = data_pc
    return data_pc_trnsfs, labels



def check_consistent_barcode_counts(data: dict) -> bool:
    """
    Verifies that all transformations across all parameter sets have
    the same number of barcodes.
    """
    reference_counts = None
    inconsistencies = []

    for params, barcodes_per_transformation in data.items():
        for t, barcodes in barcodes_per_transformation.items():
            count = len(barcodes)
            key = (params,t)

            if reference_counts is None:
                reference_counts = count
            elif count != reference_counts:
                inconsistencies.append((key, count))

    if inconsistencies:
        logger.warning("Inconsistencies found in the number of barcodes per parameter sets and transformations:")
        for (params, t), count in inconsistencies:
            logger.warning(f"Parameters: {params}, Transformation: {t}, Count: {count} (Expected: {reference_counts})")
        return False

    return True



def barcodes_equal_unordered(b1: np.ndarray, b2: np.ndarray, rtol=1e-5, atol=1e-8) -> bool:
    if b1.shape != b2.shape:
        return False
    # Sort rows for both barcodes (by birth then death)
    b1_sorted = np.array(sorted(b1.tolist()))
    b2_sorted = np.array(sorted(b2.tolist()))
    return np.allclose(b1_sorted, b2_sorted, rtol=rtol, atol=atol)



def compare_barcodes_per_transformations(d1, d2, rtol=1e-3, atol=1e-3, verbose=True):
    if d1.keys() != d2.keys():
        if verbose:
            print("Different keys: ", d1.keys() ^ d2.keys())
        return False

    equal = True
    for key in d1:
        barcodes1 = d1[key]
        barcodes2 = d2[key]

        if len(barcodes1) != len(barcodes2):
            if verbose:
                print(f"Different number of barcodes for transformation '{key}': {len(barcodes1)} vs {len(barcodes2)}")
            equal = False
            continue

        for i, (b1, b2) in enumerate(zip(barcodes1, barcodes2)):
            if not barcodes_equal_unordered(b1, b2, rtol=rtol, atol=atol):
                if verbose:
                    print(f"Difference at transformation '{key}', barcode index {i}")
                equal = False
                print(f"b1: \n{b1}")
                print(f"\nb2: \n{b2}")
                exit()

    return equal




class PreprocessingCache:
    def __init__(self):
        self.cache = {}

    def get_or_compute_from_fn(self, fn, *args, **kwargs):
        if fn is None:
            return None
        key = fn.__name__
        if key not in self.cache:
            self.cache[key] = fn(*args, **kwargs)
        return self.cache[key]



def calculate_alpha_pde1_per_trnsf(datasets_per_transformation_sorted, trnsfs = TurkevsTransformation.shortnames()):
    return {
        t: ph.calculate_pds_point_clouds(datasets_per_transformation_sorted[t],
                                        fil_complex="alpha",
                                        fil_fun="dtm",
                                        m=0.03,
                                        p=1)[1]
        for t in trnsfs
    }



def calculate_simple_pde1_per_trnsf(pde1_per_trnsf):
    if pde1_per_trnsf is None:
        raise ValueError("Data error: alpha_pde1 data is None.")

    return {
        t: ph.sorted_lifespans_pds(pde1_per_trnsf[t], size=10) for t in pde1_per_trnsf
        }



def calculate_distance_matrix_per_trnsf(datasets_per_transformation_sorted, trnsfs = TurkevsTransformation.shortnames()):
    return data_construction.calculate_distance_matrices_flat_under_trnsfs(
        point_clouds_trnsfs=datasets_per_transformation_sorted,
        trnsfs = trnsfs)



def calculate_3d_point_clouds_per_trnsf(datasets_per_transformation_sorted, trnsfs = TurkevsTransformation.shortnames()):
    return data_construction.calculate_3d_point_clouds_under_trnsfs(
        point_clouds_trnsfs = datasets_per_transformation_sorted,
        trnsfs = trnsfs)



def make_consecutive_labels(labels):
    return data_construction.encode_labels(labels)



class CModels(Enum):
    PHE = auto()
    PH = auto()
    PH_simple = auto()
    ML = auto()
    NN_shallow = auto()
    NN_deep = auto()
    PointNet = auto()



pipeline_registry = {
    CModels.PHE: {
        "data_fn": None,
        "labels_fn": None,
        "model_module": ph_ml,
    },
    CModels.PH: {
        "data_fn": calculate_alpha_pde1_per_trnsf,
        "labels_fn": None,
        "model_module": ph_ml,
    },
    CModels.PH_simple: {
        "data_fn": calculate_alpha_pde1_per_trnsf,
        "data_postprocessing_fn": calculate_simple_pde1_per_trnsf,
        "labels_fn": None,
        "model_module": ml,
    },
    CModels.ML: {
        "data_fn": calculate_distance_matrix_per_trnsf,
        "labels_fn": make_consecutive_labels,
        "model_module": ml,
    },
    CModels.NN_shallow: {
        "data_fn": calculate_distance_matrix_per_trnsf,
        "labels_fn": make_consecutive_labels,
        "model_module": nn_shallow,
    },
    CModels.NN_deep: {
        "data_fn": calculate_distance_matrix_per_trnsf,
        "labels_fn": make_consecutive_labels,
        "model_module": nn_deep,
    },
    CModels.PointNet: {
        "data_fn": calculate_3d_point_clouds_per_trnsf,
        "labels_fn": make_consecutive_labels,
        "model_module": point_net,
    },
}



def make_parameter_filter(**criteria):
    def filter_fn(params):
        for attr, expected_value in criteria.items():
            value = getattr(params, attr, None)
            if isinstance(expected_value, (list, tuple, set)):
                if value not in expected_value:
                    return False
            else:
                if value != expected_value:
                    return False
            return True
    return filter_fn



# def run_classification_experiments(experiment_summaries: list[ExperimentSummary],
#                                    datasets: list[Dataset],
#                                    pipelines: list = [m for m in CModels]):
#     datasets_per_transformation_sorted = get_datasets_per_transformation_sorted(datasets)
#     labels = get_labels_sorted(datasets)
#     train_indices, test_indices = get_train_and_test_indices(datasets_per_transformation_sorted)

#     evaluation_results = []

#     if CModels.PHE in pipelines:
#         barcodes_dim_1_per_parameters = get_barcodes_dim_1_per_parameters_per_transformation_sorted(experiment_summaries)
#         summaries_consistent = check_consistent_barcode_counts(barcodes_dim_1_per_parameters)
#         if not summaries_consistent:
#             raise ValueError("Inconsistent barcode counts in experiment summaries.")

#         for parameters, pde_1_per_transformation in barcodes_dim_1_per_parameters.items():
#             pde_accuracy = evaluate_model(pde_1_per_transformation, labels, train_indices, test_indices, ph_ml, CModels.PHE.name)
#             pde_accuracy.parameters = parameters
#             evaluation_results.append(pde_accuracy)

#         # rips_from_ellipsoids = next(
#         #     (
#         #         barcodes_per_transformations
#         #         for params, barcodes_per_transformations in barcodes_dim_1_per_parameters.items()
#         #         if params.complex_type == ComplexType.BALL and params.complex_subtype == ComplexSubtype.RIPS
#         #     ),
#         #     None
#         #     )

#     trnsfs = [t.shortname for t in TurkevsTransformation]

#     if any(m in pipelines for m in [CModels.PH, CModels.PH_simple]):
#         data_pd1_per_transformation = {
#             t: ph.calculate_pds_point_clouds(datasets_per_transformation_sorted[t],
#                                             fil_complex="alpha",
#                                             fil_fun="dtm",
#                                             m=0.03,
#                                             p=1)[1]
#             for t in trnsfs
#         }

#         if CModels.PH in pipelines:
#             pd_accuracy = evaluate_model(data_pd1_per_transformation, labels, train_indices, test_indices, ph_ml, CModels.PH.name)
#             evaluation_results.append(pd_accuracy)



#         if CModels.PH_simple in pipelines:
#             data_ph_simple_per_transformation = {
#                 t: ph.sorted_lifespans_pds(data_pd1_per_transformation[t], size=10) for t in trnsfs }
#             ph_simple_accuracy = evaluate_model(data_ph_simple_per_transformation, labels, train_indices, test_indices, ml, CModels.PH_simple.name)
#             evaluation_results.append(ph_simple_accuracy)

#     if any(m in pipelines for m in [CModels.ML, CModels.NN_shallow, CModels.NN_deep]):
#         data_dis_mat_flat_trnsfs = data_construction.calculate_distance_matrices_flat_under_trnsfs(
#             point_clouds_trnsfs=datasets_per_transformation_sorted,
#             trnsfs=trnsfs)
#         labels_con, _ = data_construction.encode_labels(labels)

#         if CModels.ML in pipelines:
#             ml_accuracy = evaluate_model(data_dis_mat_flat_trnsfs, labels_con, train_indices, test_indices, ml, CModels.ML.name)
#             evaluation_results.append(ml_accuracy)
#         if CModels.NN_shallow in pipelines:
#             nn_shallow_accuracy = evaluate_model(data_dis_mat_flat_trnsfs, labels_con, train_indices, test_indices, nn_shallow, CModels.NN_shallow.name)
#             evaluation_results.append(nn_shallow_accuracy)
#         if CModels.NN_deep in pipelines:
#             nn_deep_accuracy = evaluate_model(data_dis_mat_flat_trnsfs, labels_con, train_indices, test_indices, nn_deep, CModels.NN_deep.name)
#             evaluation_results.append(nn_deep_accuracy)

#     if CModels.PointNet in pipelines:
#         data_pc_3d_trnsfs = data_construction.calculate_3d_point_clouds_under_trnsfs(point_clouds_trnsfs = datasets_per_transformation_sorted, trnsfs = trnsfs)
#         labels_con, _ = data_construction.encode_labels(labels)
#         point_net_accuracy = evaluate_model(data_pc_3d_trnsfs, labels_con, train_indices, test_indices, point_net, CModels.PointNet.name)
#         evaluation_results.append(point_net_accuracy)

#     return evaluation_results



@dataclass
class AggregatedResults():
    mean_accuracy_per_trnsf: dict
    sd_accuracy_per_trnsf: dict
    n_runs: int
    pipeline_name: str
    results_label: str
    parameters: Optional[Parameters] = None
    # dataset_id: Optional[str] = None

    def to_dict(self):
            return {
                "parameters": self.parameters.to_dict() if self.parameters is not None else None,
                "mean_accuracy_per_trnsf": self.mean_accuracy_per_trnsf,
                "sd_accuracy_per_trnsf": self.sd_accuracy_per_trnsf,
                "n_runs": self.n_runs,
                "pipeline_name": self.pipeline_name,
                "results_label": self.results_label
            }



def calculate_aggregated_classification_results(results_folder: str):
    paths = get_paths_of_files_in_a_folder(results_folder, ".json")
    if len(paths) == 0:
        raise ValueError("Error: no valid paths found.")

    evaluation_results = []
    for path in paths:
        ev_res_dicts = read_from_json(path)
        evaluation_results.extend(EvaluationResults.from_dict(ev_res_dict) for ev_res_dict in ev_res_dicts)

    grouped_results = defaultdict(list)
    for ev_res in evaluation_results:
        key = (ev_res.pipeline_name, ev_res.parameters)
        grouped_results[key].append(ev_res.accuracy)

    aggregated_results = []
    for (pipeline_name, params), acc_list in grouped_results.items():
        trnsf = acc_list[0].keys()
        acc_mean = {t: np.mean([acc[t] for acc in acc_list]) for t in trnsf}
        acc_sd = {t: np.std([acc[t] for acc in acc_list]) for t in trnsf}
        aggregated_results.append(AggregatedResults(parameters=params,
                                                                 mean_accuracy_per_trnsf=acc_mean,
                                                                 sd_accuracy_per_trnsf=acc_sd,
                                                                 n_runs=len(acc_list),
                                                                 pipeline_name=pipeline_name,
                                                                 results_label=pipeline_name))
    return sort_aggregated_results(aggregated_results)


def sort_aggregated_results(aggregated_results, pipeline_order=[CModels.PHE, CModels.PH, CModels.PH_simple, CModels.NN_shallow, CModels.NN_deep, CModels.PointNet]):
    pipeline_order_names = [p.name for p in pipeline_order]
    def sort_key(res):
        pipeline_priority = pipeline_order_names.index(res.pipeline_name) if res.pipeline_name in pipeline_order_names else len(pipeline_order)
        return (pipeline_priority, str(res.parameters))  # Secondary sort by parameters
    return sorted(aggregated_results, key=sort_key)



def add_aggregate_results_labels(aggregated_results):
    name_counters = defaultdict(int)
    for agg_result in aggregated_results:
        name_counters[agg_result.pipeline_name] += 1

    name_indices = defaultdict(int)
    for agg_result in aggregated_results:
        pipeline_name = agg_result.pipeline_name
        name_indices[pipeline_name] += 1
        index = name_indices[pipeline_name]
        results_label = f"{pipeline_name} [{index}]" if name_counters[pipeline_name] > 1 else pipeline_name
        agg_result.results_label = results_label

    return aggregated_results



def plot_aggregated_results(aggregated_results_per_datatype, output_folder, print_parameters=True, base_filename="aggregated_results"):
    transformations = [t.fullname for t in TurkevsTransformation]
    plot_data = {}
    err_data = {}
    params_for_label = {}
    min_n_runs = np.inf
    for agg_res in aggregated_results_per_datatype:
        min_n_runs = agg_res.n_runs if agg_res.n_runs < min_n_runs else min_n_runs
        plot_data[agg_res.results_label] = [agg_res.mean_accuracy_per_trnsf[t.shortname]
                                            for t in TurkevsTransformation]  # y-values
        err_data[agg_res.results_label] = [agg_res.sd_accuracy_per_trnsf[t.shortname]
                                            for t in TurkevsTransformation]  # y-values
        params_for_label[agg_res.results_label] = agg_res.parameters  # for legend
    average_accs = plot_data
    results_labels = list(plot_data.keys())

    title = f"Mean classification accuracy per transformation (≥{min_n_runs} runs)\n\n"
    plots.plot_bar_chart(transformations, average_accs, results_labels, err_data, title)
    plt.savefig(os.path.join(output_folder, base_filename), bbox_inches = "tight")

    # plot_legend_file = "accs_trnsfs_averages_legend.txt"
    plot_legend_file = f"{base_filename}_legend.txt"
    plot_parameters = f"{title} legend:\n\n"
    for plot_label, parameters in params_for_label.items():
        if parameters is not None:
            plot_parameters += f"{plot_label}:\n{parameters.to_dict()}\n\n"

    with open (os.path.join(output_folder, plot_legend_file), "w") as f:
        f.write(plot_parameters)

    if print_parameters:
        print(plot_parameters)


from typing import Iterator


def iter_experiment_summaries_from_compressed_jsonl(filepath: str) -> Iterator[ExperimentSummary]:
    if not os.path.exists(filepath):
        if os.path.exists(filepath + ".gz"):
            filepath += ".gz"
        else:
            raise FileNotFoundError(f"No such file: '{filepath}'.")
    with gzip.open(filepath, "rt", encoding="utf-8") as f:
        for line in f:
            try:
                data = json.loads(line)
                yield ExperimentSummary.from_dict(data)
            except json.JSONDecodeError as e:
                logger.warning(f"Skipping invalid JSON line: {e}")



def read_experiment_summaries_from_jsonl(filepath: str) -> list[ExperimentSummary]:
    return list(iter_experiment_summaries_from_compressed_jsonl(filepath))



    # if not os.path.exists(filepath):
    #     if os.path.exists(filepath + ".gz"):
    #         filepath += ".gz"
    #     else:
    #         raise FileNotFoundError(f"No such file: '{filepath}' or '{filepath}.gz'")
    # experiment_summaries = []
    # with gzip.open(filepath, "rt", encoding="utf-8") as f:
    #     for line in f:
    #         try:
    #             data = json.loads(line)
    #             experiment_summaries.append(ExperimentSummary.from_dict(data))
    #         except json.JSONDecodeError as e:
    #             logger.warning(f"Skipping invalid JSON line: {e}")
    # return experiment_summaries



def flush_buffer(buffer: list, output_path: str, compress: bool = True):
    if not buffer:
        return

    open_func = gzip.open if compress else open
    mode = "at" if compress else "a"
    output_path = output_path + ".gz" if compress else output_path

    with open_func(output_path, mode, encoding="utf-8") as f:
        for entry in buffer:
            f.write(json.dumps(entry, cls=CustomEncoder, separators=(",", ":")) + "\n")
        buffer.clear()



def convert_json_files_to_jsonl(folder: str, output_path: str, compress: bool = True):
    paths = get_paths_of_files_in_a_folder(folder, "json")
    logger.info(f"Found {len(paths)} files to convert to jsonl.")

    open_func = gzip.open if compress else open
    mode = "wt" if compress else "w"  # text mode for writing strings
    output_path = output_path + ".gz" if compress else output_path

    with open_func(output_path, mode, encoding="utf-8") as f_out:
        for path in paths:
            try:
                with open(path, "r", encoding="utf-8") as in_f:
                    data = json.load(in_f)
                f_out.write(json.dumps(data, separators=(",", ":")) + "\n")
            except OSError as e:
                logger.error(f"I/O error while processing {path}: {e}")
