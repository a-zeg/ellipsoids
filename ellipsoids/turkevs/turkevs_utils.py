#!/usr/bin/env python3

import os
import sys
import re
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional, Any, cast, Iterator
from collections import defaultdict
import logging
import hashlib
import gzip
import json

sys.path.append(os.path.abspath('.'))

import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split

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
from ellipsoids.data_handling import CustomEncoder, sample_multiscale_loops
from ellipsoids.data_handling import read_from_json
from ellipsoids.data_handling import get_paths_of_files_in_a_folder
# from ellipsoids.data_handling import check_type
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
# from ellipsoids.turkevs.additional_datasets import
from ellipsoids.data_handling import (
    sample_from_figure_eight,
    sample_spiral,
    sample_3d_bottleneck,
    sample_braided_theta,
    sample_folded_ladder,
    sample_multiscale_loops,
)
from ellipsoids.turkevs.data_construction import TURKEVS_CORE_SHAPES

from functools import partial

logger = logging.getLogger(__name__)



def generate_turkevs_results_folder_name(n_point_clouds: int, n_points: int, seed: int):
    return f"npc={n_point_clouds}_np={n_points}_s={seed}"



def generate_variations_of_point_cloud_():
    pass


# ADDITIONAL_SHAPES = [
#     (sample_spiral, 0),
#     (sample_3d_bottleneck, 1),
#     (sample_from_figure_eight, 1),

#     (partial(sample_braided_theta, holes=2), 2),
#     (partial(sample_braided_theta, holes=4), 4),
#     (partial(sample_braided_theta, holes=9), 9),

#     (partial(sample_folded_ladder, holes=2), 2),
#     (partial(sample_folded_ladder, holes=4), 4),
#     (partial(sample_folded_ladder, holes=9), 9),

#     (partial(sample_multiscale_loops, holes=2), 2),
#     (partial(sample_multiscale_loops, holes=4), 4),
#     (partial(sample_multiscale_loops, holes=9), 9),
# ]
ADDITIONAL_SHAPES = [
    ("spiral", sample_spiral, 0),
    ("3d_bottleneck", sample_3d_bottleneck, 1),
    ("figure_eight", sample_from_figure_eight, 1),

    ("braided_theta_2", partial(sample_braided_theta, holes=2), 2),
    ("braided_theta_4", partial(sample_braided_theta, holes=4), 4),
    ("braided_theta_9", partial(sample_braided_theta, holes=9), 9),

    ("folded_ladder_2", partial(sample_folded_ladder, holes=2), 2),
    ("folded_ladder_4", partial(sample_folded_ladder, holes=4), 4),
    ("folded_ladder_9", partial(sample_folded_ladder, holes=9), 9),

    ("multiscale_loops_2", partial(sample_multiscale_loops, holes=2), 2),
    ("multiscale_loops_4", partial(sample_multiscale_loops, holes=4), 4),
    ("multiscale_loops_9", partial(sample_multiscale_loops, holes=9), 9),
]



SHAPE_GROUPS = {
    "turkevs": TURKEVS_CORE_SHAPES,
    "additional": ADDITIONAL_SHAPES,
}


def get_n_shapes(shape_groups):
    unknown_groups = set(shape_groups) - set(SHAPE_GROUPS)

    if unknown_groups:
        raise ValueError(
            f"Unknown shape groups: {sorted(unknown_groups)}"
        )

    if len(shape_groups) != len(set(shape_groups)):
        raise ValueError("shape_groups contains duplicates.")

    return sum(
        len(SHAPE_GROUPS[group])
        for group in shape_groups
    )



# def additional_datasets_holes(n_point_clouds_per_dataset: int, n_points: int):
#     datasets = []
#     labels = []

#     logger.info("Generating additional datasets...")

#     for sampling_function, label in ADDITIONAL_SHAPES:
#         for _ in range(n_point_clouds_per_dataset):
#             datasets.append(sampling_function(n_points))
#             labels.append(label)

#     return datasets, labels



def additional_datasets_holes(
    n_shape_variations: int,
    n_points: int,
):
    datasets = []
    labels = []
    shape_names = []
    shape_variation_indices = []

    logger.info("Generating additional datasets...")

    for shape_name, sampling_function, label in ADDITIONAL_SHAPES:
        for variation_index in range(n_shape_variations):
            datasets.append(
                sampling_function(n_points)
            )
            labels.append(label)
            shape_names.append(shape_name)
            shape_variation_indices.append(variation_index)

    return (
        datasets,
        np.asarray(labels),
        np.asarray(shape_names),
        np.asarray(shape_variation_indices),
    )


# N_TURKEVS_SHAPES = 20



# def generate_turkevs_datasets(n_point_clouds: int,
#                               n_points: int,
#                               seed: int,
#                               dataset_id="",
#                               ):




# def generate_turkevs_datasets(
#         n_shape_variations,
#         n_points,
#         seed,
#         shape_groups=("turkevs", "additional"),
#         dataset_id="",
# ):
    
#     np.random.seed(seed)
#     logger.info("Constructing Turkevs datasets...")

#     initial_point_clouds = []
#     labels = []

#     if "turkevs" in shape_groups:
#         n_core_point_clouds = (
#             len(TURKEVS_CORE_SHAPES) * n_shape_variations
#         )

#         core_point_clouds, core_labels, _ = build_dataset_holes(
#             n_core_point_clouds,
#             n_points,
#         )

#         initial_point_clouds.extend(core_point_clouds)
#         labels.extend(core_labels)

#     # initial_point_clouds, labels, _ = build_dataset_holes(n_point_clouds, n_points)

#     # n_point_clouds_per_dataset = int(n_point_clouds / N_TURKEVS_SHAPES)
#     # additional_point_clouds, additional_labels = \
#     #     additional_datasets_holes(n_point_clouds_per_dataset, n_points)

#     # initial_point_clouds.extend(additional_point_clouds)
#     # labels = np.concatenate((labels, additional_labels))


#     if "additional" in shape_groups:
#         additional_point_clouds, additional_labels = (
#             additional_datasets_holes(
#                 n_shape_variations,
#                 n_points,
#             )
#         )

#         initial_point_clouds.extend(additional_point_clouds)
#         labels.extend(additional_labels)
    
#     labels = np.asarray(labels)

#     # for i, pc in enumerate(initial_point_clouds):
#     #     print(type(pc))
#     #     print(labels[i])
#     #     print("\n")

#     all_datasets: list[Dataset] = []
#     n_shapes = get_n_shapes(shape_groups)
#     n_standard_point_clouds = n_shapes * n_shape_variations

#     if len(initial_point_clouds) != n_standard_point_clouds:
#         raise ValueError(
#             f"Expected {n_standard_point_clouds} standard point clouds "
#             f"({n_shapes} shapes × {n_shape_variations} variations), "
#             f"but generated {len(initial_point_clouds)}."
#         )

#     # if n_point_clouds != len(initial_point_clouds):
#     #     raise ValueError("Mismatch in the number of point clouds.")

#     # for transformation in TurkevsTransformation:
#     #     logger.info(f"Applying transformation: {transformation.fullname}.")
#     #     transformed_point_clouds = (
#     #         initial_point_clouds if transformation is TurkevsTransformation.STANDARD
#     #         else calculate_point_clouds_under_trnsf(initial_point_clouds,
#     #                                                 transformation=transformation.fullname)
#     #     )
#     #     for mesh_index, (points,label) in enumerate(zip(transformed_point_clouds, labels)):
#     #         dataset = Dataset(
#     #             points=points,
#     #             data_type="turkevs",
#     #             additional_info=TurkevsDatasetInfo(
#     #                 dataset_id=dataset_id,
#     #                 seed=seed,
#     #                 point_cloud_index=mesh_index,
#     #                 transformation=transformation,
#     #                 label=label,
#     #                 n_point_clouds=n_standard_point_clouds
#     #             )
#     #         )
#     #         all_datasets.append(dataset)

#     for transformation in TurkevsTransformation:
#         logger.info(f"Applying transformation: {transformation.fullname}.")

#         transformed_point_clouds = (
#             initial_point_clouds
#             if transformation is TurkevsTransformation.STANDARD
#             else calculate_point_clouds_under_trnsf(
#                 initial_point_clouds,
#                 transformation=transformation.fullname,
#             )
#         )

#         for point_cloud_index, (points, label) in enumerate(
#             zip(transformed_point_clouds, labels)
#         ):
#             dataset = Dataset(
#                 points=points,
#                 data_type="turkevs",
#                 additional_info=TurkevsDatasetInfo(
#                     dataset_id=dataset_id,
#                     seed=seed,
#                     point_cloud_index=point_cloud_index,
#                     transformation=transformation,
#                     label=label,

#                     n_shapes=n_shapes,
#                     n_shape_variations=n_shape_variations,
#                     n_point_clouds=n_standard_point_clouds,
#                 ),
#             )
#             all_datasets.append(dataset)


    
#     logger.info("Datasets constructed.")
#     return all_datasets






def generate_turkevs_datasets(
    n_shape_variations: int,
    n_points: int,
    seed: int,
    shape_groups=("turkevs", "additional"),
    dataset_id="",
):
    np.random.seed(seed)
    logger.info("Constructing Turkevs datasets...")

    shape_groups = tuple(shape_groups)

    n_shapes = get_n_shapes(shape_groups)
    n_standard_point_clouds = n_shapes * n_shape_variations

    initial_point_clouds = []
    labels = []
    shape_names = []
    shape_groups_per_cloud = []
    shape_variation_indices = []

    if "turkevs" in shape_groups:
        (
            core_point_clouds,
            core_labels,
            core_shape_names,
            core_variation_indices,
        ) = build_dataset_holes(
            n_shape_variations=n_shape_variations,
            n=n_points,
        )

        initial_point_clouds.extend(core_point_clouds)
        labels.extend(core_labels)
        shape_names.extend(core_shape_names)
        shape_groups_per_cloud.extend(
            ["turkevs"] * len(core_point_clouds)
        )
        shape_variation_indices.extend(core_variation_indices)

    if "additional" in shape_groups:
        (
            additional_point_clouds,
            additional_labels,
            additional_shape_names,
            additional_variation_indices,
        ) = additional_datasets_holes(
            n_shape_variations=n_shape_variations,
            n_points=n_points,
        )

        initial_point_clouds.extend(additional_point_clouds)
        labels.extend(additional_labels)
        shape_names.extend(additional_shape_names)
        shape_groups_per_cloud.extend(
            ["additional"] * len(additional_point_clouds)
        )
        shape_variation_indices.extend(
            additional_variation_indices
        )

    labels = np.asarray(labels)

    if len(initial_point_clouds) != n_standard_point_clouds:
        raise ValueError(
            f"Expected {n_standard_point_clouds} standard point clouds "
            f"({n_shapes} shapes x {n_shape_variations} variations), "
            f"but generated {len(initial_point_clouds)}."
        )


    all_datasets: list[Dataset] = []

    for transformation in TurkevsTransformation:
        logger.info(
            f"Applying transformation: {transformation.fullname}."
        )

        transformed_point_clouds = (
            initial_point_clouds
            if transformation is TurkevsTransformation.STANDARD
            else calculate_point_clouds_under_trnsf(
                initial_point_clouds,
                transformation=transformation.fullname,
            )
        )

        for point_cloud_index, (
            points,
            label,
            shape_name,
            shape_group,
            shape_variation_index,
        ) in enumerate(
            zip(
                transformed_point_clouds,
                labels,
                shape_names,
                shape_groups_per_cloud,
                shape_variation_indices,
            )
        ):
            dataset = Dataset(
                points=points,
                data_type="turkevs",
                additional_info=TurkevsDatasetInfo(
                    dataset_id=dataset_id,
                    seed=seed,
                    point_cloud_index=point_cloud_index,
                    transformation=transformation,
                    label=label,

                    shape_group=shape_group,
                    shape_name=shape_name,
                    shape_variation_index=int(
                        shape_variation_index
                    ),

                    n_shapes=n_shapes,
                    n_shape_variations=n_shape_variations,

                    # Keeping for backwards compatibility.
                    # Now has a clear meaning:
                    # number of STANDARD point clouds.
                    n_point_clouds=n_standard_point_clouds,
                ),
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



def read_turkevs_datasets(datasets_path):
    logger.info(f"Reading in Turkevs datasets from {datasets_path}...")
    datasets_dicts = read_from_json(datasets_path)
    datasets = [Dataset.from_dict(d) for d in datasets_dicts]
    datasets = ensure_turkevs_datasets_consistency(datasets)
    logger.info("Datasets read.")
    return datasets



def ensure_turkevs_datasets_consistency(datasets: list[Dataset]):
    # dataset_id = check_type(datasets[0].additional_info, TurkevsDatasetInfo).dataset_id
    dataset_id = cast(TurkevsDatasetInfo, datasets[0].additional_info).dataset_id
    for dataset in datasets:
        dataset_info = cast(TurkevsDatasetInfo, dataset.additional_info)
        if dataset_id is not None and dataset_info.dataset_id != dataset_id:
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



# def divide_indices_into_train_and_test(indices, train_percentage=0.8, test_percentage=0.2):
#     n_indices = len(indices)
#     train_size = int(train_percentage * n_indices)
#     test_size = int(test_percentage * n_indices)
#     train_indices = np.random.choice(indices, size=train_size, replace=False)
#     non_train_indices = np.setdiff1d(indices, train_indices)
#     test_indices = np.random.choice(non_train_indices, size=test_size, replace=False)
#     return train_indices, test_indices
def divide_indices_into_train_and_test(
    indices,
    stratify_values,
    train_percentage=0.8,
    test_percentage=0.2,
    random_state=None,
):
    train_indices, test_indices = train_test_split(
        indices,
        train_size=train_percentage,
        test_size=test_percentage,
        stratify=stratify_values,
        random_state=random_state,
        shuffle=True,
    )

    # Sorting is convenient because predictions will then correspond
    # to increasing point_cloud_index.
    train_indices = np.sort(train_indices)
    test_indices = np.sort(test_indices)

    return train_indices, test_indices


# @dataclass
# class EvaluationResults:
#     accuracy: dict
#     pipeline_name: str
#     parameters: Optional[Parameters] = None

#     def to_dict(self):
#         result = {
#             "accuracy": self.accuracy,
#             "pipeline_name": self.pipeline_name,
#         }
#         if self.parameters is not None:
#             result["parameters"] = self.parameters.to_dict()
#         return result

#     @classmethod
#     def from_dict(cls, data: dict):
#         parameters = None
#         if "parameters" in data and data["parameters"] is not None:
#             parameters = Parameters.from_any_dict(data["parameters"])
#         return cls(
#             accuracy=data["accuracy"],
#             pipeline_name=data["pipeline_name"],
#             parameters=parameters
#         )

@dataclass
class EvaluationResults:
    accuracy: dict
    pipeline_name: str
    parameters: Optional[Parameters] = None

    predictions: Optional[dict] = None
    train_indices: Optional[list[int]] = None
    test_indices: Optional[list[int]] = None
    labels_test: Optional[list] = None
    split_seed: Optional[int] = None

    def to_dict(self):
        result = {
            "accuracy": self.accuracy,
            "pipeline_name": self.pipeline_name,
            "predictions": self.predictions,
            "train_indices": self.train_indices,
            "test_indices": self.test_indices,
            "labels_test": self.labels_test,
            "split_seed": self.split_seed,
        }

        if self.parameters is not None:
            result["parameters"] = self.parameters.to_dict()

        return result

    @classmethod
    def from_dict(cls, data: dict):
        parameters = None

        if (
            "parameters" in data
            and data["parameters"] is not None
        ):
            parameters = Parameters.from_any_dict(
                data["parameters"]
            )

        return cls(
            accuracy=data["accuracy"],
            pipeline_name=data["pipeline_name"],
            parameters=parameters,

            # .get() keeps OLD classification files readable.
            predictions=data.get("predictions"),
            train_indices=data.get("train_indices"),
            test_indices=data.get("test_indices"),
            labels_test=data.get("labels_test"),
            split_seed=data.get("split_seed"),
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



def get_shape_names_sorted(datasets: list[Dataset]) -> np.ndarray:
    """
    Returns shape names for the standard datasets, sorted by point_cloud_index.
    """
    standard_datasets = [
        d for d in datasets
        if isinstance(d.additional_info, TurkevsDatasetInfo)
        and d.additional_info.transformation == TurkevsTransformation.STANDARD
    ]

    standard_datasets_sorted = sorted(
        standard_datasets,
        key=lambda d: cast(
            TurkevsDatasetInfo,
            d.additional_info
        ).point_cloud_index,
    )

    shape_names = [
        cast(TurkevsDatasetInfo, d.additional_info).shape_name
        for d in standard_datasets_sorted
    ]

    if any(shape_name is None for shape_name in shape_names):
        raise ValueError(
            "Cannot stratify by shape_name because some datasets "
            "do not contain shape_name metadata."
        )

    return np.asarray(shape_names)



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
            # sorted_barcodes = [np.asarray(filter_barcode(barcode=barcode, dim=1))
            #                    for _, barcode in sorted(point_cloud_index_and_barcodes,
            #                                             key=lambda x: x[0])]
            # padded_barcodes = pad_to_max_length(sorted_barcodes)
            # sorted_and_filtered[params][t] = padded_barcodes
            #
            sorted_barcodes = []

            for _, barcode in sorted(
                    point_cloud_index_and_barcodes,
                    key=lambda x: x[0],
            ):
                pd = np.asarray(
                    filter_barcode(barcode=barcode, dim=1),
                    dtype=float,
                )

                # Important for empty H1 diagrams:
                # keep them as valid (0, 2) persistence diagrams.
                if pd.size == 0:
                    pd = np.empty((0, 2), dtype=float)
                else:
                    pd = pd.reshape(-1, 2)

                sorted_barcodes.append(pd)

            sorted_and_filtered[params][t] = sorted_barcodes

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



# def test_model(model_trained,
#                data_test,
#                labels_test):
#     logger.info("Testing the model by evaluating noise robustness...")
#     trnsfs = [transformation.shortname for transformation in TurkevsTransformation]
#     accs_trnsfs = get_scores_under_trnsfs(data_test,trnsfs,labels_test,model_trained)
#     logger.info("Testing completed.")
#     return accs_trnsfs

def test_model(
    model_trained,
    data_test,
    labels_test,
):
    logger.info(
        "Testing the model by evaluating noise robustness..."
    )

    trnsfs = [
        transformation.shortname
        for transformation in TurkevsTransformation
    ]

    accs_trnsfs = {}
    predictions_trnsfs = {}

    for t in trnsfs:
        data = data_test[t]

        accs_trnsfs[t] = get_score(
            data,
            labels_test,
            model_trained,
        )

        predictions = model_trained.predict(data)

        predictions_trnsfs[t] = (
            np.asarray(predictions).tolist()
        )

    logger.info("Testing completed.")

    return accs_trnsfs, predictions_trnsfs



# def evaluate_model(
#         data_per_transformations,
#         labels,
#         train_indices,
#         test_indices,
#         ml_module,
#         pipeline_name
#         ) -> EvaluationResults:
def evaluate_model(
        data_per_transformations,
        labels,
        train_indices,
        test_indices,
        ml_module,
        pipeline_name,
        split_seed=None,
) -> EvaluationResults:
    logger.info(f"Running pipeline {pipeline_name}...")
    labels_train = labels[train_indices]
    labels_test = labels[test_indices]
    data_train = [data_per_transformations[TurkevsTransformation.STANDARD.shortname][i] for i in train_indices]
    data_test_per_transformations = filter_turkevs_dict(data_per_transformations, test_indices)

    if isinstance(data_per_transformations["std"], np.ndarray):
        data_train = np.asarray(data_train)
        data_test_per_transformations = {t: np.asarray(data_test_per_transformations[t]) for t in data_test_per_transformations}

    # model_trained = train_model(ml_module, data_train, labels_train)
    # accs_trnsfs = test_model(model_trained, data_test_per_transformations, labels_test)
    # return EvaluationResults(pipeline_name=pipeline_name, accuracy=accs_trnsfs)
    model_trained = train_model(
        ml_module,
        data_train,
        labels_train,
    )

    accs_trnsfs, predictions_trnsfs = test_model(
        model_trained,
        data_test_per_transformations,
        labels_test,
    )

    return EvaluationResults(
        pipeline_name=pipeline_name,
        accuracy=accs_trnsfs,
        predictions=predictions_trnsfs,
        # Convert numpy values to normal Python values for JSON.
        train_indices=np.asarray(train_indices).astype(int).tolist(),
        test_indices=np.asarray(test_indices).astype(int).tolist(),
        labels_test=np.asarray(labels_test).tolist(),
        split_seed=split_seed,
    )


# def get_train_and_test_indices(datasets_per_transformation: dict["str", list]):
#     logger.info("Generating train and test indices...")
#     n_std_point_clouds = len(datasets_per_transformation[TurkevsTransformation.STANDARD.shortname])
#     indices = np.arange(n_std_point_clouds)
#     logger.info("Train and test indices generated.")
#     train_indices, test_indices = divide_indices_into_train_and_test(indices)
#     return train_indices, test_indices
def get_train_and_test_indices(
    datasets_per_transformation,
    stratify_values,
    random_state=None,
):
    logger.info("Generating stratified train and test indices...")

    n_std_point_clouds = len(
        datasets_per_transformation[
            TurkevsTransformation.STANDARD.shortname
        ]
    )

    if len(stratify_values) != n_std_point_clouds:
        raise ValueError(
            "stratify_values and standard point clouds "
            "have different lengths."
        )

    indices = np.arange(n_std_point_clouds)

    train_indices, test_indices = (
        divide_indices_into_train_and_test(
            indices=indices,
            stratify_values=stratify_values,
            random_state=random_state,
        )
    )

    logger.info("Train and test indices generated.")

    return train_indices, test_indices


# def point_clouds_from_scratch(n_point_clouds, n_points) -> tuple[dict["str", list], np.ndarray]:
#     from ellipsoids.turkevs.data_construction import build_dataset_holes
#     data_pc, labels, _ = build_dataset_holes(n_point_clouds,n_points)
#     data_pc_trnsfs: dict["str", list] = {t.shortname: data_construction.calculate_point_clouds_under_trnsf(data_pc, t.fullname)
#                       for t in TurkevsTransformation if t is not TurkevsTransformation.STANDARD}
#     data_pc_trnsfs[TurkevsTransformation.STANDARD.shortname] = data_pc
#     return data_pc_trnsfs, labels



def point_clouds_from_scratch(
    n_shape_variations,
    n_points,
) -> tuple[dict[str, list], np.ndarray]:

    data_pc, labels, _, _ = build_dataset_holes(
        n_shape_variations=n_shape_variations,
        n=n_points,
    )

    data_pc_trnsfs = {
        t.shortname:
            data_construction.calculate_point_clouds_under_trnsf(
                data_pc,
                t.fullname,
            )
        for t in TurkevsTransformation
        if t is not TurkevsTransformation.STANDARD
    }

    data_pc_trnsfs[
        TurkevsTransformation.STANDARD.shortname
    ] = data_pc

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



# def calculate_alpha_pde1_per_trnsf(datasets_per_transformation_sorted, trnsfs = TurkevsTransformation.shortnames()):
#     return {
#         t: ph.calculate_pds_point_clouds(datasets_per_transformation_sorted[t],
#                                         fil_complex="alpha",
#                                         fil_fun="dtm",
#                                         m=0.03,
#                                         p=1)[1]
#         for t in trnsfs
#     }


def calculate_alpha_pde1_per_trnsf(
    datasets_per_transformation_sorted,
    trnsfs=TurkevsTransformation.shortnames(),
):
    return {
        t: ph.calculate_pds_point_clouds(
            datasets_per_transformation_sorted[t],
            fil_complex="alpha",
            fil_fun="dtm",
            m=0.03,
            p=1,
            pad=False,
        )[1]
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
    return data_construction.encode_labels(labels)[0]



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



# def make_parameter_filter(**criteria):
#     def filter_fn(params):
#         for attr, expected_value in criteria.items():
#             value = getattr(params, attr, None)
#             if isinstance(expected_value, (list, tuple, set)):
#                 if value not in expected_value:
#                     return False
#             else:
#                 if value != expected_value:
#                     return False
#         return True
#     return filter_fn

# def make_parameter_filter(**criteria):
#     def filter_fn(params):
#         for attr, expected_value in criteria.items():
#             value = getattr(params, attr, None)

#             if expected_value is None:
#                 if value is not None:
#                     return False

#             elif isinstance(expected_value, (list, tuple, set)):
#                 # Handle array comparisons
#                 if isinstance(value, np.ndarray):
#                     found_match = any(np.array_equal(value, ev) for ev in expected_value)
#                     if not found_match:
#                         return False
#                 else:
#                     if value not in expected_value:
#                         return False

#             else:
#                 if isinstance(value, np.ndarray) and isinstance(expected_value, np.ndarray):
#                     if not np.array_equal(value, expected_value):
#                         return False
#                 else:
#                     if value != expected_value:
#                         return False

#         return True

#     return filter_fn

def make_parameter_filter(**criteria):
    def filter_fn(params):
        for attr, expected_value in criteria.items():
            value = getattr(params, attr, None)

            if expected_value is None:
                if value is not None:
                    return False

            elif isinstance(expected_value, (list, tuple, set)):
                if isinstance(value, np.ndarray):
                    if not any(np.array_equal(value, ev) for ev in expected_value):
                        return False
                else:
                    if value not in expected_value:
                        return False

            else:
                if isinstance(value, np.ndarray) or isinstance(expected_value, np.ndarray):
                    if not np.array_equal(value, expected_value):
                        return False
                else:
                    if value != expected_value:
                        return False

        return True
    return filter_fn




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



def calculate_aggregated_classification_results(results_folder: str) -> list[AggregatedResults]:
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



def plot_aggregated_results(aggregated_results_per_datatype,
                            output_path,
                            print_parameters=True):
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

    title = f"Mean classification accuracy per transformation ({min_n_runs} runs)\n\n"
    plots.plot_bar_chart(transformations, average_accs, results_labels, err_data, title)
    plt.savefig(output_path, bbox_inches = "tight")

    plot_legend_path = f"{output_path}_legend.txt"
    plot_parameters = f"{title} legend:\n\n"
    for plot_label, parameters in params_for_label.items():
        if parameters is not None:
            plot_parameters += f"{plot_label}:\n{parameters.to_dict()}\n\n"

    with open (plot_legend_path, "w") as f:
        f.write(plot_parameters)

    if print_parameters:
        print(plot_parameters)



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



def flush_buffer(buffer: list, output_path: str, compress: bool = True):
    if not buffer:
        return

    print(output_path)

    open_func = gzip.open if compress else open
    mode = "at" if compress else "a"
    output_path = output_path + ".json" if os.path.splitext(output_path)[1] == "" else output_path
    output_path = output_path + ".gz" if compress and os.path.splitext(output_path)[1] != ".gz" else output_path

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



import plotly.graph_objects as go
from plotly.subplots import make_subplots
        
def plot_turkevs_datasets(
    datasets: list[Dataset],
    marker_size: float = 3,
    height_per_row: int = 250,
    width_per_col: int = 800,
):
    """
    Plot Turkevs point clouds in a grid.

    Rows:
        transformations, in TurkevsTransformation enum order

    Columns:
        point_cloud_index

    The number of point clouds / columns is inferred automatically
    from the supplied Dataset objects.

    Parameters
    ----------
    datasets:
        List of Dataset objects with TurkevsDatasetInfo as additional_info.

    marker_size:
        Size of Plotly markers.

    height_per_row:
        Approximate subplot height.

    width_per_col:
        Approximate subplot width.

    Returns
    -------
    plotly.graph_objects.Figure
    """

    # ------------------------------------------------------------
    # Keep only Turkevs datasets
    # ------------------------------------------------------------

    datasets = [
        dataset
        for dataset in datasets
        if dataset.data_type == "turkevs"
        and isinstance(dataset.additional_info, TurkevsDatasetInfo)
    ]

    if not datasets:
        raise ValueError("No Turkevs datasets found.")

    transformations = list(TurkevsTransformation)

    # ------------------------------------------------------------
    # Determine point clouds / columns dynamically
    # ------------------------------------------------------------

    point_cloud_indices = sorted({
        dataset.additional_info.point_cloud_index
        for dataset in datasets
    })

    n_columns = len(point_cloud_indices)
    n_rows = len(transformations)

    # Map actual point_cloud_index -> Plotly column number
    column_for_index = {
        point_cloud_index: i + 1
        for i, point_cloud_index in enumerate(point_cloud_indices)
    }

    # ------------------------------------------------------------
    # Build lookup:
    #
    # (point_cloud_index, transformation) -> Dataset
    # ------------------------------------------------------------

    dataset_lookup = {}

    for dataset in datasets:
        info = dataset.additional_info

        key = (
            info.point_cloud_index,
            info.transformation,
        )

        if key in dataset_lookup:
            raise ValueError(
                "Multiple datasets found for "
                f"point_cloud_index={info.point_cloud_index}, "
                f"transformation={info.transformation.fullname}."
            )

        dataset_lookup[key] = dataset

    # ------------------------------------------------------------
    # Determine dimensionality of each column
    # ------------------------------------------------------------

    dimensions = {}

    for point_cloud_index in point_cloud_indices:

        column_datasets = [
            dataset
            for dataset in datasets
            if dataset.additional_info.point_cloud_index == point_cloud_index
        ]

        dims = {
            dataset.ambient_dim
            for dataset in column_datasets
        }

        if not dims <= {2, 3}:
            raise ValueError(
                f"Point cloud {point_cloud_index} has unsupported "
                f"dimension(s): {dims}"
            )

        if len(dims) != 1:
            raise ValueError(
                f"Point cloud {point_cloud_index} mixes dimensions: {dims}"
            )

        dimensions[point_cloud_index] = dims.pop()

    # ------------------------------------------------------------
    # Subplot types
    # ------------------------------------------------------------

    specs = []

    for _ in transformations:
        row_specs = []

        for point_cloud_index in point_cloud_indices:
            dimension = dimensions[point_cloud_index]

            row_specs.append(
                {"type": "xy"}
                if dimension == 2
                else {"type": "scene"}
            )

        specs.append(row_specs)

    # ------------------------------------------------------------
    # Column labels
    # ------------------------------------------------------------

    column_titles = [
        f"Point cloud {idx}"
        for idx in point_cloud_indices
    ]

    # ------------------------------------------------------------
    # Create figure
    # ------------------------------------------------------------

    fig = make_subplots(
        rows=n_rows,
        cols=n_columns,
        specs=specs,
        column_titles=column_titles,
        horizontal_spacing=0.03,
        vertical_spacing=0.03,
    )

    # ------------------------------------------------------------
    # Add traces
    # ------------------------------------------------------------

    for row_idx, transformation in enumerate(transformations, start=1):

        for point_cloud_index in point_cloud_indices:

            col_idx = column_for_index[point_cloud_index]

            key = (
                point_cloud_index,
                transformation,
            )

            # Allow missing transformations
            if key not in dataset_lookup:
                continue

            dataset = dataset_lookup[key]
            points = np.asarray(dataset.points)

            # -------------------------
            # 2D
            # -------------------------

            if dataset.ambient_dim == 2:

                fig.add_trace(
                    go.Scattergl(
                        x=points[:, 0],
                        y=points[:, 1],
                        mode="markers",
                        marker=dict(
                            size=marker_size,
                        ),
                        showlegend=False,
                        customdata=np.full(
                            len(points),
                            transformation.fullname,
                            dtype=object,
                        ),
                        hovertemplate=(
                            f"Point cloud: {point_cloud_index}<br>"
                            f"Transformation: {transformation.fullname}<br>"
                            "x=%{x:.4f}<br>"
                            "y=%{y:.4f}"
                            "<extra></extra>"
                        ),
                    ),
                    row=row_idx,
                    col=col_idx,
                )

                fig.update_xaxes(
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    title_text=None,
                    row=row_idx,
                    col=col_idx,
                )

                fig.update_yaxes(
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    title_text=None,
                    row=row_idx,
                    col=col_idx,
                )

            # -------------------------
            # 3D
            # -------------------------

            elif dataset.ambient_dim == 3:

                fig.add_trace(
                    go.Scatter3d(
                        x=points[:, 0],
                        y=points[:, 1],
                        z=points[:, 2],
                        mode="markers",
                        marker=dict(
                            size=marker_size,
                            opacity=0.8,
                        ),
                        showlegend=False,
                        hovertemplate=(
                            f"Point cloud: {point_cloud_index}<br>"
                            f"Transformation: {transformation.fullname}<br>"
                            "x=%{x:.4f}<br>"
                            "y=%{y:.4f}<br>"
                            "z=%{z:.4f}"
                            "<extra></extra>"
                        ),
                    ),
                    row=row_idx,
                    col=col_idx,
                )

    # ------------------------------------------------------------
    # Clean up 3D scenes
    # ------------------------------------------------------------

    for key in fig.layout:
        if str(key).startswith("scene"):
            fig.layout[key].update(
                aspectmode="data",
                xaxis=dict(
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    showbackground=False,
                    title="",
                ),
                yaxis=dict(
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    showbackground=False,
                    title="",
                ),
                zaxis=dict(
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    showbackground=False,
                    title="",
                ),
            )

    # ------------------------------------------------------------
    # Row labels
    # ------------------------------------------------------------

    for i, transformation in enumerate(transformations):

        y = 1 - (i + 0.5) / n_rows

        fig.add_annotation(
            x=-0.015,
            y=y,
            xref="paper",
            yref="paper",
            text=transformation.fullname.title(),
            showarrow=False,
            xanchor="right",
            yanchor="middle",
            textangle=-90,
            font=dict(size=13),
        )

    # ------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------

    fig.update_layout(
        width=max(
            700,
            width_per_col * n_columns,
        ),
        height=height_per_row * n_rows,
        margin=dict(
            l=100,
            r=20,
            t=70,
            b=30,
        ),
        template="plotly_white",
    )

    return fig
