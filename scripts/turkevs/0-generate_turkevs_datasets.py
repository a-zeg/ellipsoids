import os
import sys
import logging

sys.path.append(os.path.abspath('.'))

from ellipsoids.turkevs.turkevs_utils import generate_turkevs_datasets
from ellipsoids.turkevs.turkevs_utils import generate_datasets_hash
from ellipsoids.turkevs.turkevs_utils import add_dataset_id
from ellipsoids.data_handling import save_to_json
from ellipsoids.turkevs.config import TURKEVS_DATA_DIR
from ellipsoids.turkevs.config import get_setup_name, get_datasets_path
from ellipsoids.logging_setup import setup_logging



setup_logging()
logger = logging.getLogger(__name__)



def main():
    n_point_clouds = 100
    n_points = 20
    seed = 0
    parent_folder = TURKEVS_DATA_DIR

    logger.info(f"Generating {n_point_clouds} point clouds with {n_points} using seed {seed}...")

    setup_name = get_setup_name(n_point_clouds=n_point_clouds,
                                n_points=n_points,
                                seed=seed)
    if any(setup_name == os.path.basename(path) for path in os.listdir(parent_folder)):
        raise FileExistsError(f"Folder {setup_name} already exists in {parent_folder}.")

    datasets = generate_turkevs_datasets(n_point_clouds=n_point_clouds,
                                        n_points=n_points,
                                        seed=seed)
    dataset_id = generate_datasets_hash(datasets)
    datasets = add_dataset_id(datasets=datasets, dataset_id=dataset_id)

    datasets_path = get_datasets_path(setup_name)
    datasets_path = save_to_json([d.to_dict() for d in datasets], filename=datasets_path, add_timestamp=False)

    logger.info(f"Datasets saved to {datasets_path}.")



if __name__ == '__main__':
    main()
