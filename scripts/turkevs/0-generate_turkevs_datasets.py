import os
import sys

sys.path.append(os.path.abspath('.'))

from ellipsoids.turkevs.turkevs_utils import generate_turkevs_datasets
from ellipsoids.turkevs.turkevs_utils import generate_datasets_hash
from ellipsoids.turkevs.turkevs_utils import add_dataset_id
from ellipsoids.turkevs.turkevs_utils import generate_turkevs_results_folder_name
from ellipsoids.data_handling import save_to_json
from ellipsoids.turkevs.config import REL_DATASETS_DIR, TURKEVS_DATA_DIR

import itertools


if __name__ == '__main__':

    npcs = [100,300]
    npts = [20,50,100,200,300]

    all_tuples = list(itertools.product(npcs, npts))
    for n_point_clouds, n_points in all_tuples:

        # n_point_clouds = 100
        # n_points = 200
        seed = 0

        parent_folder = TURKEVS_DATA_DIR

        results_folder = generate_turkevs_results_folder_name(n_point_clouds=n_point_clouds, n_points=n_points, seed=seed)
        paths = os.listdir(parent_folder)
        if any(results_folder in path for path in paths):
            raise FileExistsError(f"Folder with parameters {results_folder} already exists in {parent_folder}.")

        datasets = generate_turkevs_datasets(n_point_clouds=n_point_clouds,
                                            n_points=n_points,
                                            seed=seed)
        dataset_id = generate_datasets_hash(datasets)
        datasets = add_dataset_id(datasets=datasets, dataset_id=dataset_id)

        filename_parameters = f"{n_point_clouds=}_{n_points=}_seed={seed}_id={dataset_id}"
        datasets_filename =  f"turkevs_datasets_{filename_parameters}"
        datasets_path = os.path.join(parent_folder, REL_DATASETS_DIR, f"turkevs_{results_folder}", datasets_filename)

        save_to_json([d.to_dict() for d in datasets], filename=datasets_path, add_timestamp=False)
