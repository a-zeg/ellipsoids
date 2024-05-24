# from  import get_paths_of_files_in_a_folder
from main import set_filename_parameters
from main import generate_filename
from main import calculate_and_save_ellipsoids_and_rips_data
import data_handling
import numpy as np
import pickle
import data_construction
import re
import os

def generate_turkevs_datasets(N: int, n: int, filename: str):
    '''
    Generate datasets for Turkevs experiments (the original point clouds +
    transformed point clouds) and save them as pkl files.

    Parameters
        ----------
        N : int
            The number of point clouds.
            There are 20 different point cloud types: 4 each with 0, 1, 2, 3, and 4 
            holes), so N should be divisible by 20.
        n : int
            The number of points in each point cloud.
        filename : str
            The name of the file in which to save the pickled point clouds.
    '''

    if not filename.endswith('.pkl'):
        filename = filename + '.pkl'

    print("\n\nConstructing the data...")
    data_pc, _, _ = data_construction.build_dataset_holes(N, n)
    data_pc_trns = data_construction.calculate_point_clouds_under_trnsf(data_pc, transformation = "translation")
    data_pc_rot = data_construction.calculate_point_clouds_under_trnsf(data_pc, transformation = "rotation")
    data_pc_stretch = data_construction.calculate_point_clouds_under_trnsf(data_pc, transformation = "stretch")
    data_pc_shear = data_construction.calculate_point_clouds_under_trnsf(data_pc, transformation = "shear")
    data_pc_gauss = data_construction.calculate_point_clouds_under_trnsf(data_pc, transformation = "gaussian")
    data_pc_out = data_construction.calculate_point_clouds_under_trnsf(data_pc, transformation = "outliers")
    data_pc_trnsfs = {}
    data_pc_trnsfs["std"] = data_pc
    data_pc_trnsfs["trns"] = data_pc_trns
    data_pc_trnsfs["rot"] = data_pc_rot
    data_pc_trnsfs["stretch"] = data_pc_stretch
    data_pc_trnsfs["shear"] = data_pc_shear
    data_pc_trnsfs["gauss"] = data_pc_gauss
    data_pc_trnsfs["out"] = data_pc_out

    with open(filename, "wb") as f:
        pickle.dump(data_pc_trnsfs, f)

def generate_id(folder):
    '''
    The datasets used in the Turkevs holes tests are generated using the code from the paper and saved in a pkl file.
    To keep track of which datasets correspond to which data / graphs, 'id' was introduced.

    This function checks the ids of all the pkl files in the given folder and returns the next available one.
    '''
    paths = data_handling.get_paths_of_files_in_a_folder(folder, extension='pkl')
    id = 0
    for path in paths:
        match = re.search(r'id=(\d+)', path) 
        if match is not None:
            current_id = int(match.group(1))
            print(current_id)
            if current_id > id:
                id = current_id
    return str(id+1).zfill(4)


if __name__ == '__main__':

    folder = 'datasets/turkevs'
    N = 100
    n = 100
    id = generate_id(folder)
    filename = 'datasets/turkevs/pc_test_trnsfs' + '_' + f'{N=}' + '_' + f'{n=}' + '_' + 'id=' + id + '.pkl'
    generate_turkevs_datasets(N, n, filename)
