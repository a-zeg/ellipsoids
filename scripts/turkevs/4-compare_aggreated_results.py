#!/usr/bin/env python3

import os
import sys
import logging

sys.path.append(os.path.abspath('.'))

from ellipsoids.turkevs.config import TURKEVS_DATA_DIR, REL_AGGREGATED_RESULTS_DIR
from ellipsoids.logging_setup import setup_logging
from ellipsoids.common import Parameters
from ellipsoids.data_handling import read_from_json
from statistics import mean
from collections import defaultdict

setup_logging()
logger = logging.getLogger(__name__)


def main():

    folders = [os.path.join(TURKEVS_DATA_DIR, folder) for folder in os.listdir(TURKEVS_DATA_DIR) if os.path.isdir(os.path.join(TURKEVS_DATA_DIR, folder))]

    mean_accumulator = defaultdict(list)
    for folder in folders:
        aggregated_results_folder = os.path.join(folder, REL_AGGREGATED_RESULTS_DIR)
        if not os.path.exists(aggregated_results_folder):
            continue
        aggregated_results_path = os.path.join(aggregated_results_folder, f"aggregated_results_{os.path.basename(folder)}.json")

        aggregated_results = read_from_json(aggregated_results_path)
        for result in aggregated_results:
            if result["parameters"] is None:
                continue
            param_obj = Parameters.from_any_dict(result["parameters"])
            mean_accs = result["mean_accuracy_per_trnsf"]
            mean_val = mean(mean_accs.values())
            mean_accumulator[param_obj].append(mean_val)

    sorted_results = sorted(
        mean_accumulator.items(),
        key=lambda item: mean(item[1]),
        reverse=True
    )

    for param, acc_list in sorted_results:
        avg_accuracy = mean(acc_list)
        print(f"Parameters: {param}")
        print(f"Average mean accuracy across datasets: {avg_accuracy:.4f}")
        print("-" * 50)

if __name__ == '__main__':
    main()
