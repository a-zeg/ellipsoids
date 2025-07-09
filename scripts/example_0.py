import os
import sys

sys.path.append(os.path.abspath('.'))

import numpy as np
from ellipsoids.common import Dataset
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import Experiment
from ellipsoids.data_handling import sample_from_circle
from ellipsoids.visualisation import plot_experiment



def main():

    dataset = Dataset(sample_from_circle(20), "circle")
    parameters = EllipsoidParameters(
        save_simplex_tree=True,
        save_ellipsoid_list=True,
        axes_ratios=np.asarray([2,1]),
        collapse_edges=False,
        )
    experiment = Experiment(dataset, parameters)

    experiment.run()
    experiment.save_to_json()
    experiment.print_execution_time()

    experiment.plot_parameters.draw_points = True
    experiment.plot_parameters.draw_ellipsoids = True
    experiment.plot_parameters.draw_simplex_tree = True
    experiment.plot_parameters.filtration = 0.8

    plot_experiment(experiment)



if __name__ == '__main__':
    main()
