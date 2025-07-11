import numpy as np
import os
import sys
from typing import cast
import logging

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import Experiment
from ellipsoids.common import DatasetSummary
from ellipsoids.turkevs.turkevs_utils import read_turkevs_datasets
from ellipsoids.turkevs.config import get_datasets_path, get_experiment_summaries_path
from ellipsoids.logging_setup import setup_logging
from ellipsoids.data_handling import build_signature_index_from_jsonl, generate_signature_from_dataset_summary_parameters, parse_args
from ellipsoids.turkevs.turkevs_utils import flush_buffer



setup_logging(level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    setup_path = "data/turkevs/turkevs_npc=100_np=20_s=0"
    check_if_already_calculated = True
    buffer_size = 10

    complexes_to_calculate = [
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=7, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=9, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=7, axes_ratios=np.asarray([2,1]), r_spherisize=0.5),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=9, axes_ratios=np.asarray([2,1]), r_spherisize=0.5),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=7, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=9, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        Parameters(complex_subtype=ComplexSubtype.RIPS),
        Parameters(complex_subtype=ComplexSubtype.ALPHA),
        ]
    # nbhd_sizes = [5,12,15,18,20,25,30,40,50,75,100]
    # complexes_to_calculate = []
    # for nbhd_size in nbhd_sizes:
    #     complexes_to_calculate.append(
    #         EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=nbhd_size, axes_ratios=None),
    #     )

    args = parse_args()
    setup_path = args.folder or setup_path
    datasets_start = args.datasets_start
    datasets_end = args.datasets_end

    logger.info(f"Calculating ellipsoids data for the folder {setup_path}...")

    datasets_path = get_datasets_path(setup_path)
    datasets = read_turkevs_datasets(datasets_path)
    summaries_jsonl_path = get_experiment_summaries_path(setup_path)

    if datasets_start is not None or datasets_end is not None:
        datasets_start = datasets_start if datasets_start is not None else 0
        datasets_end = datasets_end if datasets_end is not None else len(datasets)
        summaries_jsonl_path = get_experiment_summaries_path(setup_path, suffix=f"_part={datasets_start}-{datasets_end}")
        datasets = datasets[datasets_start:datasets_end]

    signature_index = {}
    if os.path.isfile(summaries_jsonl_path):
        signature_index = build_signature_index_from_jsonl(summaries_jsonl_path) if check_if_already_calculated else {}

    buffer = [] # for saving data in chunks

    for dataset in datasets:
        dataset_info = cast(TurkevsDatasetInfo, dataset.additional_info)
        logger.info(f"Dataset transformation: "
                + f"{dataset_info.transformation.fullname}, "
                + f"point cloud index: {dataset_info.point_cloud_index}.")

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
            buffer.append(experiment.serialize_summary())
            if len(buffer) >= buffer_size:
                flush_buffer(buffer, summaries_jsonl_path, compress=True)

    flush_buffer(buffer, summaries_jsonl_path, compress=True)
    logger.info(f"Results saved to {summaries_jsonl_path}.")



if __name__ == '__main__':
    main()
