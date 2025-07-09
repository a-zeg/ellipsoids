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
from ellipsoids.data_handling import sample_from_cassini_oval
from ellipsoids.data_handling import sample_from_circle
from ellipsoids.visualisation import plot_experiments



if __name__ == '__main__':
    dataset = Dataset(sample_from_figure_eight(30, a=1, b=0.3), "figure_eight")
    dataset = Dataset(sample_from_figure_eight(100, a=1, b=0.3), "figure_eight")
    # dataset = Dataset(sample_from_circle(30), "circle")
    # dataset = Dataset(sample_from_circle(100), "circle")

    all_parameters = [
        # EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS,
        #                     axes_ratios=None,
        #                     nbhd_size=3),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS,
                            axes_ratios=None,
                            nbhd_size=5),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS,
                            axes_ratios=None,
                            nbhd_size=10),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS,
                            axes_ratios=np.array([2,1])),
        ]

    all_experiments = []
    for parameters in all_parameters:
        parameters.save_simplex_tree = True
        parameters.collapse_edges=False

        experiment = Experiment(dataset, parameters)
        experiment.run()
        experiment.plot_parameters.draw_points = True
        experiment.plot_parameters.draw_ellipsoids = True
        experiment.plot_parameters.draw_simplex_tree = True
        experiment.plot_parameters.filtration = 1

        all_experiments.append(experiment)

    plot_experiments(all_experiments)
