import os
import sys

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Experiment
from ellipsoids.visualisation import plot_experiments
from ellipsoids.data_handling import get_paths_of_files_in_a_folder



if __name__ == '__main__':

    folder = "data/test_read"

    experiments = []
    filepaths = get_paths_of_files_in_a_folder(folder)
    for filepath in filepaths:
        experiments.append(Experiment.read_from_json(filepath))

    plot_experiments(experiments)
