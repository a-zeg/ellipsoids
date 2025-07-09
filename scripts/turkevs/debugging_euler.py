import numpy as np
import os
import sys
from typing import cast
import argparse
import logging

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import TurkevsTransformation
from ellipsoids.common import Experiment
from ellipsoids.turkevs.turkevs_utils import read_turkevs_datasets
from ellipsoids.turkevs.config import REL_EXPERIMENT_SUMMARIES_DIR, REL_DATASETS_DIR
from ellipsoids.logging_setup import setup_logging



setup_logging()
logger = logging.getLogger(__name__)


def main():
    # folder = "data/turkevs/turkevs_npc=100_np=20_s=0"
    folder = "data/turkevs/turkevs_npc=300_np=100_s=0"
    complexes_to_calculate = [
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([3,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([3,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA, nbhd_size=3, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA, nbhd_size=5, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        Parameters(complex_subtype=ComplexSubtype.RIPS),
        Parameters(complex_subtype=ComplexSubtype.ALPHA),
        ]


    # when running as a batch job, one can also specify the folder with a --folder flag
    parser = argparse.ArgumentParser(description="Calculate ellipsoid barcodes on a dataset folder.")
    parser.add_argument('--folder', type=str, default=folder,
                        help="Path to the folder containing the datasets and where the results will be stored.")
    args = parser.parse_args()
    folder = args.folder

    logger.info(f"Calculating ellipsoids data for the folder {folder}...")

    datasets_folder = os.path.join(folder, REL_DATASETS_DIR)
    datasets_path = [os.path.join(datasets_folder, file) for file in os.listdir(datasets_folder) if "datasets" in file][0] # assuming only one file with "datasets" per folder
    datasets = read_turkevs_datasets(datasets_path)
    summaries_folder = os.path.join(folder, REL_EXPERIMENT_SUMMARIES_DIR)

    for dataset in datasets:
        dataset_info = cast(TurkevsDatasetInfo, dataset.additional_info)
        # print(f"trnsf = {dataset_info.transformation}; index = {dataset_info.point_cloud_index}")
        # print(f"trnsf = {dataset_info.transformation}; index = {dataset_info.point_cloud_index}")
        # if dataset_info.transformation == TurkevsTransformation.STANDARD and dataset_info.point_cloud_index == 15:
        if dataset_info.transformation == TurkevsTransformation.STANDARD and dataset_info.point_cloud_index == 14:
            print("HELLO")
            logger.info(f"Dataset transformation: "
                    + f"{dataset_info.transformation.fullname}, "
                    + f"point cloud index: {dataset_info.point_cloud_index}.")

            for i, params in enumerate(complexes_to_calculate):
                logger.info(f"Parameter set {i+1} of {len(complexes_to_calculate)}...")
                experiment = Experiment(dataset, params)
                experiment.run()
                experiment.save_summary(folder=summaries_folder)


    logger.info(f"Results saved to {summaries_folder}.")



def circle():

    complexes_to_calculate = [
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=None),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([3,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([3,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA, nbhd_size=3, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.ALPHA, nbhd_size=5, axes_ratios=np.asarray([2,1])),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=3, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        EllipsoidParameters(complex_subtype=ComplexSubtype.RIPS, nbhd_size=5, axes_ratios=np.asarray([2,1]), r_spherisize=1),
        Parameters(complex_subtype=ComplexSubtype.RIPS),
        Parameters(complex_subtype=ComplexSubtype.ALPHA),
        ]


    from ellipsoids.data_handling import sample_from_circle
    from ellipsoids.common import Dataset

    dataset = Dataset(sample_from_circle(200), "circle")



    for i, params in enumerate(complexes_to_calculate):
        logger.info(f"Parameter set {i+1} of {len(complexes_to_calculate)}...")
        experiment = Experiment(dataset, params)
        experiment.run()





if __name__ == '__main__':
    main()
    # circle()







    # from sklearn.decomposition import PCA

    # nbhd_pts = np.asarray([[1,0], [0,0], [0,0], [0,0]])

    # pca = PCA(n_components=len(nbhd_pts[0]))
    # pca.fit(nbhd_pts)
    # axes = pca.components_

    # singular_values_scaled = pca.singular_values_ / pca.singular_values_[0]

    # # if np.any(singular_values_scaled == np.nan):
    # #     singular_values_scaled



    # print(f"axes = {axes}")
    # print(f"singular_values_scaled = {singular_values_scaled}")
