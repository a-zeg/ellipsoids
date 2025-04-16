import os
import sys

sys.path.append(os.path.abspath('.'))

import numpy as np
from ellipsoids.common import Dataset
from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import Experiment
from ellipsoids.data_handling import sample_from_figure_eight
from ellipsoids.visualisation import plot_experiments



if __name__ == '__main__':
    dataset = Dataset(sample_from_figure_eight(30, a=1, b=0.3), "figure_eight")

    all_parameters = [Parameters(complex_subtype=ComplexSubtype.RIPS),
                      EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS,
                                          axes_ratios=np.array([3,1]))]

    all_experiments = []
    for parameters in all_parameters:
        parameters.save_simplex_tree = True
        parameters.collapse_edges=False

        experiment = Experiment(dataset, parameters)
        experiment.run()
        experiment.plot_parameters.draw_points = True
        experiment.plot_parameters.draw_ellipsoids = True
        experiment.plot_parameters.draw_simplex_tree = True
        experiment.plot_parameters.filtration = 0.7

        all_experiments.append(experiment)

    plot_experiments(all_experiments)
