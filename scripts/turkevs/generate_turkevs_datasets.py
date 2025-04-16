import os
import sys

sys.path.append(os.path.abspath('.'))

from ellipsoids.data_handling import generate_turkevs_datasets



if __name__ == '__main__':

    n_point_clouds = 20
    n_points = 20

    folder = os.path.join("datasets", "turkevs")
    seed = 0

    generate_turkevs_datasets(n_point_clouds=n_point_clouds,
                              n_points=n_points,
                              seed=seed,
                              folder=folder,
                              save_to_file=True)
