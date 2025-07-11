import os
import sys
import logging

sys.path.append(os.path.abspath('.'))

import numpy as np
from ellipsoids.data_handling import save_to_json
from ellipsoids.data_handling import parse_args
from ellipsoids.turkevs.turkevs_utils import add_aggregate_results_labels, plot_aggregated_results
from ellipsoids.turkevs.turkevs_utils import calculate_aggregated_classification_results
from ellipsoids.turkevs.config import get_aggregated_results_basepath, get_classification_results_folder, get_classification_results_path
from ellipsoids.logging_setup import setup_logging
from ellipsoids.turkevs.turkevs_utils import make_parameter_filter
from ellipsoids.common import ComplexType, ComplexSubtype


setup_logging()
logger = logging.getLogger(__name__)


def main():

    setup_path = "data/turkevs/turkevs_npc=100_np=300_s=0"

    filter_fn = make_parameter_filter(
        complex_type=[
            ComplexType.BALL,
            ComplexType.ELLIPSOID
            ],
    )

    args = parse_args()
    setup_path = args.folder or setup_path
    logger.info(f"Calculating aggregated results for the setup {setup_path}...")

    classification_results_folder = get_classification_results_folder(setup_path)
    aggregated_results = calculate_aggregated_classification_results(classification_results_folder)
    filtered_results = []
    for results in aggregated_results:
        if results.parameters is None: filtered_results.append(results)
        if filter_fn is not None and filter_fn(results.parameters):
            filtered_results.append(results)
    filtered_results = add_aggregate_results_labels(filtered_results)

    aggregated_results_path = get_aggregated_results_basepath(setup_path)
    save_to_json([a_res.to_dict() for a_res in filtered_results], aggregated_results_path, add_timestamp=False)
    plot_aggregated_results(filtered_results, aggregated_results_path, print_parameters=False)



if __name__ == '__main__':
    main()
