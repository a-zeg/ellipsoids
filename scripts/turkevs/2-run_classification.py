import os
import sys
import numpy as np
import logging

sys.path.append(os.path.abspath('.'))

from ellipsoids.data_handling import save_to_json
from ellipsoids.turkevs.turkevs_utils import find_files_by_keyword
from ellipsoids.turkevs.turkevs_utils import CModels
from ellipsoids.turkevs.config import get_classification_results_path, get_experiment_summaries_path
from ellipsoids.turkevs.turkevs_utils import read_turkevs_datasets
from ellipsoids.turkevs.turkevs_utils import read_experiment_summaries_from_jsonl
from ellipsoids.logging_setup import setup_logging
from ellipsoids.data_handling import parse_args

from ellipsoids.turkevs.turkevs_utils import get_datasets_per_transformation_sorted
from ellipsoids.turkevs.turkevs_utils import get_labels_sorted
from ellipsoids.turkevs.turkevs_utils import get_train_and_test_indices
from ellipsoids.turkevs.turkevs_utils import get_barcodes_dim_1_per_parameters_per_transformation_sorted
from ellipsoids.turkevs.turkevs_utils import check_consistent_barcode_counts
from ellipsoids.turkevs.turkevs_utils import evaluate_model
from ellipsoids.turkevs.turkevs_utils import PreprocessingCache
from ellipsoids.turkevs.turkevs_utils import pipeline_registry
from ellipsoids.turkevs.config import get_datasets_path


setup_logging()
logger = logging.getLogger(__name__)


def main():

    setup_path = "data/turkevs/turkevs_npc=100_np=20_s=0"
    pipelines = [
        CModels.PHE,
        CModels.PH,
        CModels.PH_simple,
        # CModels.ML,
        # CModels.NN_shallow,
        # CModels.NN_deep,
        # CModels.PointNet,
    ]

    args = parse_args()
    setup_path = args.folder or setup_path
    logger.info(f"Calculating ellipsoids data for the folder {setup_path}...")

    datasets_path = get_datasets_path(setup_path)
    datasets = read_turkevs_datasets(datasets_path)
    datasets_per_transformation_sorted = get_datasets_per_transformation_sorted(datasets)
    labels = get_labels_sorted(datasets)
    train_indices, test_indices = get_train_and_test_indices(datasets_per_transformation_sorted)
    special_pipelines = [CModels.PHE if CModels.PHE in pipelines else None]
    regular_pipelines = [p for p in pipelines if p not in special_pipelines]

    evaluation_results = []
    if CModels.PHE in special_pipelines:
        summaries_path = get_experiment_summaries_path(setup_path)
        experiment_summaries = read_experiment_summaries_from_jsonl(summaries_path)
        barcodes_dim_1_per_parameters = get_barcodes_dim_1_per_parameters_per_transformation_sorted(experiment_summaries)
        summaries_consistent = check_consistent_barcode_counts(barcodes_dim_1_per_parameters)
        if not summaries_consistent:
            raise ValueError(f"Data incomplete: inconsistent barcode counts in experiment summaries in {setup_path}.")

        model_key = CModels.PHE
        current_pipeline_setup = pipeline_registry[model_key]
        for parameters, pde_1_per_transformation in barcodes_dim_1_per_parameters.items():
            pde_accuracy = evaluate_model(
                pde_1_per_transformation,
                labels,
                train_indices,
                test_indices,
                current_pipeline_setup["model_module"],
                model_key.name)
            pde_accuracy.parameters = parameters
            evaluation_results.append(pde_accuracy)

    cache = PreprocessingCache()
    for model_key in regular_pipelines:
        current_pipeline_setup = pipeline_registry[model_key]
        data_fn = current_pipeline_setup["data_fn"]
        data_post_fn = current_pipeline_setup.get("data_postprocessing_fn", None)
        labels_fn = current_pipeline_setup["labels_fn"]
        model_module = current_pipeline_setup["model_module"]

        raw_data = cache.get_or_compute_from_fn(data_fn, datasets_per_transformation_sorted)
        if raw_data is None:
            raise ValueError(f"Data error in {model_key} pipeline.")
        data = data_post_fn(raw_data) if data_post_fn else raw_data

        labels = labels if labels_fn is None else cache.get_or_compute_from_fn(labels_fn, labels)

        accuracy = evaluate_model(
            data,
            labels,
            train_indices,
            test_indices,
            model_module,
            model_key.name)
        evaluation_results.append(accuracy)

    classification_results_path = get_classification_results_path(setup_path)
    save_to_json([ev_res.to_dict() for ev_res in evaluation_results],
                 classification_results_path,
                 add_timestamp=True)



if __name__ == '__main__':

    n_runs = 10
    for _ in np.arange(n_runs):
        main()
