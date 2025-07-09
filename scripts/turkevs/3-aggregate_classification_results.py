import os
import sys
import logging

sys.path.append(os.path.abspath('.'))

import numpy as np
from ellipsoids.data_handling import save_to_json
from ellipsoids.data_handling import parse_args
from ellipsoids.turkevs.turkevs_utils import add_aggregate_results_labels, plot_aggregated_results
from ellipsoids.turkevs.turkevs_utils import calculate_aggregated_classification_results
from ellipsoids.turkevs.config import TURKEVS_DATA_DIR, REL_AGGREGATED_RESULTS_DIR, REL_CLASSIFICATION_RESULTS_DIR
from ellipsoids.logging_setup import setup_logging
from ellipsoids.turkevs.turkevs_utils import make_parameter_filter
from ellipsoids.common import ComplexType, ComplexSubtype


setup_logging()
logger = logging.getLogger(__name__)


def main(folder):
    # folder = "data/turkevs/turkevs_npc=100_np=100_s=0"
    filter_fn = make_parameter_filter(
        complex_type=[
            ComplexType.BALL,
            ComplexType.ELLIPSOID
            ],
    )

    args = parse_args()
    folder = args.folder or folder
    logger.info(f"Calculating ellipsoids data for the folder {folder}...")

    classification_results_folder = os.path.join(folder, REL_CLASSIFICATION_RESULTS_DIR)
    aggregated_results_folder = os.path.join(folder, REL_AGGREGATED_RESULTS_DIR)

    aggregated_results = calculate_aggregated_classification_results(classification_results_folder)
    filtered_results = []
    for results in aggregated_results:
        if results.parameters is None: filtered_results.append(results)
        if filter_fn is not None and filter_fn(results.parameters):
            filtered_results.append(results)
    filtered_results = add_aggregate_results_labels(filtered_results)

    base_filename = f"aggregated_results_{os.path.basename(folder)}"
    save_to_json([a_res.to_dict() for a_res in filtered_results], os.path.join(aggregated_results_folder, base_filename), add_timestamp=False)
    plot_aggregated_results(filtered_results, aggregated_results_folder, print_parameters=False, base_filename=base_filename)



if __name__ == '__main__':
    # main()
    # folders = [os.path.join(TURKEVS_DATA_DIR, folder) for folder in os.listdir(TURKEVS_DATA_DIR) if os.path.isdir(os.path.join(TURKEVS_DATA_DIR, folder))]
    folder = "data/turkevs/turkevs_npc=100_np=300_s=0"
    folders = [folder]

    for folder in folders:
        try:
            main(folder)
        except Exception as e:
            print(f"Error processing folder '{folder}': {e}")
            continue
