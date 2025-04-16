import os
import sys
sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Dataset
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import Experiment
from ellipsoids.data_handling import sample_from_circle
from ellipsoids.visualisation import plot_experiment


if __name__ == '__main__':
    dataset = Dataset(sample_from_circle(20), "circle")
    parameters = EllipsoidParameters(save_simplex_tree=True)
    experiment = Experiment(dataset, parameters)

    experiment.run()
    experiment.save_to_json()

    experiment.plot_parameters.draw_points = True
    experiment.plot_parameters.draw_ellipsoids = True
    experiment.plot_parameters.draw_simplex_tree = True
    plot_experiment(experiment)
