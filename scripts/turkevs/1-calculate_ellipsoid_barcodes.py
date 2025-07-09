import numpy as np
import os
import sys
from typing import cast
import logging
from math import floor, sqrt

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import Experiment
from ellipsoids.common import DatasetSummary
from ellipsoids.turkevs.turkevs_utils import read_turkevs_datasets
from ellipsoids.turkevs.config import REL_DATASETS_DIR, REL_SUMMARIES_DIR
from ellipsoids.logging_setup import setup_logging

from ellipsoids.data_handling import build_signature_index_from_jsonl, generate_signature_from_dataset_summary_parameters, parse_args
from ellipsoids.turkevs.turkevs_utils import flush_buffer



setup_logging(level=logging.DEBUG)
logger = logging.getLogger(__name__)



def main():

    # parameters
    folder = "data/turkevs/turkevs_npc=100_np=300_s=0"
    check_if_already_calculated = True
    save_results_to_one_file = True
    compress = True
    buffer_size = 10

    fixed_complexes = [
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=7, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=9, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=7, axes_ratios=np.asarray([2,1]), r_spherisize=0.5),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=9, axes_ratios=np.asarray([2,1]), r_spherisize=0.5),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=7, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=9, axes_ratios=None),
        # EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([2,1])),
        # EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1])),
        # # EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([3,1])),
        # EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([3,1])),
        # EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA, nbhd_size=3, axes_ratios=np.asarray([2,1])),
        # EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA, nbhd_size=5, axes_ratios=np.asarray([2,1])),
        # EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        Parameters(complex_subtype=ComplexSubtype.RIPS),
        Parameters(complex_subtype=ComplexSubtype.ALPHA),
        ]


    # handling user input through flags, setting paths, slicing datasets if applicable
    args = parse_args()
    folder = args.folder or folder

    logger.info(f"Calculating ellipsoids data for the folder {folder}...")

    datasets_folder = os.path.join(folder, REL_DATASETS_DIR)
    datasets_path = [os.path.join(datasets_folder, file) for file in os.listdir(datasets_folder) if "datasets" in file][0] # assuming only one file with "datasets" per folder
    datasets = read_turkevs_datasets(datasets_path)

    summaries_folder = os.path.join(folder, REL_SUMMARIES_DIR)
    summaries_basename = f"summaries_{os.path.basename(folder)}"
    summaries_jsonl_path = os.path.join(summaries_folder, f"{summaries_basename}.jsonl")

    datasets_start = args.datasets_start
    datasets_end = args.datasets_end
    if datasets_start is not None or datasets_end is not None:
        datasets_start = datasets_start if datasets_start is not None else 0
        datasets_end = datasets_end if datasets_end is not None else len(datasets)
        summaries_jsonl_path = os.path.join(summaries_folder, f"{summaries_basename}_part={datasets_start}-{datasets_end}.jsonl")
        datasets = datasets[datasets_start:datasets_end]

    # creating a list of hashes of (ellipsoid) parameters / datasets that have already been already calculated
    signature_index = {}
    if os.path.isfile(summaries_jsonl_path):
        signature_index = build_signature_index_from_jsonl(summaries_jsonl_path) if check_if_already_calculated else {}

    buffer = [] # for saving data in chunks

    for dataset in datasets:
        dataset_info = cast(TurkevsDatasetInfo, dataset.additional_info)
        logger.info(f"Dataset transformation: "
                + f"{dataset_info.transformation.fullname}, "
                + f"point cloud index: {dataset_info.point_cloud_index}.")

        # parameters for the variable complex
        nbhd_size_var = max(3, floor(sqrt(dataset.n_points)))
        variable_complex = EllipsoidParameters(
            complex_subtype=ComplexSubtype.RIPS,
            nbhd_size=nbhd_size_var,
            axes_ratios=np.asarray([2,1])
            )
        complexes_to_calculate = [variable_complex] + fixed_complexes

        dataset_summary = None
        if check_if_already_calculated:
            dataset_summary = DatasetSummary(
                    data_type = dataset.data_type,
                    n_points = dataset.n_points,
                    additional_info = dataset.additional_info
                )

        for i, params in enumerate(complexes_to_calculate):
            if check_if_already_calculated:
                signature = generate_signature_from_dataset_summary_parameters(dataset_summary, params)
                if signature in signature_index:
                    continue

            logger.info(f"Parameter set {i+1} of {len(complexes_to_calculate)}...")
            experiment = Experiment(dataset, params)
            experiment.run()

            if save_results_to_one_file:
                buffer.append(experiment.serialize_summary())
                if len(buffer) >= buffer_size:
                    flush_buffer(buffer, summaries_jsonl_path, compress=compress)
            else:
                experiment.save_summary(folder=summaries_folder)

    if save_results_to_one_file:
        flush_buffer(buffer, summaries_jsonl_path, compress=compress)
        results_path = summaries_jsonl_path + ".gz" if compress else summaries_jsonl_path
        logger.info(f"Results saved to {results_path}.")
    else:
        logger.info(f"Results saved to {summaries_folder}.")



if __name__ == '__main__':
    main()
