"""
UNDER CONSTRUCTION
"""


'''
This script is the streamlined version of https://github.com/renata-turkes/turkevs2022on/blob/main/experiments/holes.py
that also takes in ellipsoids data and calculates only the accuracies.

This script expects as input the name of the folder where the ellipsoids json files are located 
and the name of the subfolder / id.

The ellipsoids persistence diagrams, as well as the (transformed) point clouds are read in from these files.

These point clouds are then used to train and test the other algorithms (this removes the possibility of calculating
ellipsoids barcodes for one dataset and then accidentally reading in a different dataset to test other algorithms).
'''


print("\n\nNUMBER OF HOLES: Point clouds ordinal classification via persistent homology (PH) and/or deep learning (DL).\n\n")

import numpy as np
import os
import sys

import matplotlib.pyplot as plt

sys.path.append(os.path.abspath('.'))

import ellipsoids.turkevs.data_construction as data_construction
import ellipsoids.turkevs.ph as ph
import ellipsoids.turkevs.model as model
import ellipsoids.turkevs.ml as ml
import ellipsoids.turkevs.ph_ml as ph_ml
import ellipsoids.turkevs.nn_shallow as nn_shallow
import ellipsoids.turkevs.nn_deep as nn_deep
import ellipsoids.turkevs.point_net as point_net
import ellipsoids.turkevs.plots as plots

import ellipsoids.data_handling as data_handling
from ellipsoids.data_handling import ensure_folder_exists, read_from_json
from ellipsoids.data_handling import get_paths_of_files_in_a_folder
from ellipsoids.common import TurkevsTransformation
from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.common import ExperimentSummary
from ellipsoids.common import Dataset


from dataclasses import dataclass
from typing import Callable, Any

class Pipeline:
    name: str
    train_data_fn: Callable[[], Any]
    test_data_fn: Callable[[], dict]
    model_fn: Callable
    labels_train: np.ndarray
    labels_test: np.ndarray



def load_all_experiment_summaries(folder: str) -> list[ExperimentSummary]:
    experiment_summaries = []
    paths = get_paths_of_files_in_a_folder(folder=folder, extension=".json")
    for path in paths:
        summary_dict = read_from_json(path)
        summary = ExperimentSummary.from_dict(summary_dict)
        experiment_summaries.append(summary)
    return experiment_summaries



def extract_barcodes_and_labels(summaries: list[ExperimentSummary]):
    barcodes = []
    labels = []
    for summary in summaries:
        barcodes.append(summary.barcode)
        label = (
            summary.dataset_summary.additional_info.label
            if isinstance(summary.dataset_summary.additional_info, TurkevsDatasetInfo)
            else -1
        )
        labels.append(label)
    return np.array(barcodes), np.array(labels)



def find_turkevs_datasets_path(id: str, folder: str = os.path.join("datasets", "turkevs")):
    paths = get_paths_of_files_in_a_folder(folder=folder, extension=".json")
    for path in paths:
        if f"_id={id}" in path:
            return path
    raise FileNotFoundError(f"No dataset file found for id={id} in folder {folder}.")



def find_turkevs_experiment_summaries_path(id: str, folder: str = os.path.join("data", "turkevs")):
    subfolder = os.path.join(folder, f"turkevs_id={id}")
    if not os.path.isdir(subfolder):
        raise FileNotFoundError(f"No results folder found for id={id} at {subfolder}.")

    print(f"Found results folder: {subfolder}")
    return subfolder



def generate_turkevs_results_path(id: str, folder: str = os.path.join("data", "turkevs")):
    subfolder = os.path.join(folder, f"turkevs_id={id}_results")

    if os.path.isdir(subfolder) and os.listdir(subfolder):
        raise FileExistsError(f"The results folder with id={id} already exists and is not empty.")

    os.makedirs(subfolder, exist_ok=True)
    print(f"Created results folder: {subfolder}")
    return subfolder



class GroupingKey:
    def __init__(self, params):
        self.params = frozenset(params.items())  # making the params immutable (needed to use instance of GroupingKey as a dictionary key)

    def __eq__(self, other):
        return self.params == other.params # comparison for dictionary keys later on

    def __hash__(self):
        return hash(self.params) # also needed for comparisons ( i think )



from collections import defaultdict

def group_experiment_summaries_by_parameters(summaries):
    grouped_summaries = defaultdict(list)
    for summary in summaries:
        params = vars(summary.parameters)
        key = GroupingKey(params)
        grouped_summaries[key].append(summary)
    return grouped_summaries







#############################

def read_turkevs_datasets(filepath: str, id: str, transformation: TurkevsTransformation = TurkevsTransformation.STANDARD):
    print(f"Reading datasets from: {filepath}")
    datasets = []

    all_datasets_dict = read_from_json(filepath)
    for dataset_dict in all_datasets_dict:
        dataset = Dataset.from_dict(dataset_dict)
        info = dataset.additional_info

        is_turkevs = isinstance(info, TurkevsDatasetInfo)
        matches_transformation = info.transformation == transformation if is_turkevs else False
        matches_id = (id is None or info.dataset_id == id) if is_turkevs else False

        if is_turkevs and matches_transformation and matches_id:
            datasets.append(dataset)

    print(f"Found {len(datasets)} datasets with transformation={transformation} and id={id}")
    if not datasets:
            raise ValueError(f"No datasets found with id={id} and transformation={transformation.fullname}.")

    return datasets



# def get_points_and_labels_from_dict(points_dict_of_lists):
#     '''
#     The variable points_dict_of_lists is a dictionary of lists.
#     This function converts it into a dictionary of numpy arrays,
#     performs consistency checks and generates labels.
#     '''
#     # converting the dictionary of lists into dictionary of numpy arrays
#     data_pc_trnsfs = {}
#     for key, value in points_dict_of_lists.items():
#         data_pc_trnsfs[key] = []
#         for element in points_dict_of_lists[key]:
#             data_pc_trnsfs[key].append(np.asarray(element))

#     labels = []
#     num_point_clouds = len(data_pc_trnsfs['std'])
#     holes_numbers = [0,1,2,4,9]

#     trnsfs = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

#     # consistency checks:
#     n_meshes_per_transformation = []
#     for transformation in trnsfs:
#         n_meshes_per_transformation.append(len(data_pc_trnsfs[transformation]))
#     if n_meshes_per_transformation.count(n_meshes_per_transformation[0]) != len(n_meshes_per_transformation):
#         exit('Data inconsistency: unequal number of meshes per transformation.')
#     try: num_point_clouds % len(holes_numbers) == 0
#     except:
#         exit('Data inconsistency: the number of point clouds is not divisible by the number of different mesh types (holes).')

#     # generating labels (assuming the points were read in in correct order)
#     n_shapes = int(num_point_clouds / len(holes_numbers))
#     for i in holes_numbers:
#         labels += [i]*n_shapes # if n_shapes = 2, creates [0, 0, 1, 1, 2, 2, 4, 4, 9, 9]
#     labels = np.asarray(labels)

#     return data_pc_trnsfs, labels



def train_model(ml_module, data_train, labels_train, model_name=""):
    print(f"\n\nTraining {model_name}...")
    print("\nTuning the hyperparameters...", end="", flush=True)
    model_ = ml_module.tune_hyperparameters(data_train, labels_train)
    print("Done.")
    print("\nTraining...")
    model_trained, _ = model.fit(data_train, labels_train, model_)
    print("Done.")
    return model_trained



def test_model(model_trained, data_test, labels_test):
    print("\n\nTesting the model: evaluating the noise robustness...", end="", flush=True)
    trnsfs = [transformation.shortname for transformation in TurkevsTransformation]
    accs_trnsfs = model.get_scores_under_trnsfs(data_test,trnsfs,labels_test,model_trained)
    print("Done")
    return accs_trnsfs



def calculate_accuracy_trnsfs(ml_module, data_train, labels_train, data_test_trnsfs, labels_test, TRAIN_SIZES, results_path='', name=''):

    model_trained = train_model(ml_module, data_train, labels_train, model_name=name)
    accs_trnsfs = test_model(model_trained, data_test_trnsfs, labels_test)

    return accs_trnsfs



def select_dict_data(data, indices):
    '''
    Data is a dictionary with keys 'std', 'trns', etc and values
    all meshes without transformations, with translation, etc.
    This function only keeps the meshes at indices given by the variable indices.
    '''
    trnsfs = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]
    filtered_data = {}
    for transformation in trnsfs:
        data_transformation = data[transformation]
        filtered_data[transformation] = [data_transformation[i] for i in indices]

    return filtered_data



def divide_indices_into_train_and_test(indices, train_percentage=0.8, test_percentage=0.2):
    n_indices = len(indices)
    train_size = int(train_percentage * n_indices)
    test_size = int(test_percentage * n_indices)
    train_indices = np.random.choice(indices, size=train_size, replace=False)
    non_train_indices = np.setdiff1d(indices, train_indices)
    test_indices = np.random.choice(non_train_indices, size=test_size, replace=False)
    return train_indices, test_indices



##########


def convert_to_turkevs_barcode_type(summaries: list) -> dict:

    result = {transformation.shortname: [] for transformation in TurkevsTransformation}

    mesh_counts = defaultdict(int)
    for s in summaries:
        t = s.dataset.additional_info.transformation
        idx = s.dataset.additional_info.mesh_index
        mesh_counts[t] = max(mesh_counts[t], idx+1)

    for t in result:
        result[t] = [None] * mesh_counts[t]

    for s in summaries:
        t = s.dataset.additional_info.transformation
        mesh_index = s.dataset.additional_info.mesh_index
        result[t][mesh_index] = s.barcodes




def run_experiments(experiment_summaries_path: str, dataset_path: str, results_path: str):
     
    # PH parameters.
    FIL_COMPLEX = "alpha"
    FIL_FUN = "dtm"
    DTM_M = 0.03
    DTM_P = 1

    # DL parameters.
    NUM_EPOCHS = 25 # 25
    BATCH_SIZE = 32
    TRAIN_SIZES = np.linspace(0.1, 1, 10) # np.linspace(0.1, 1, 10)

    trnsfs = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]


    print("\n\nChoice of hyperparameters: ")
    print("FIL_COMPLEX = ", FIL_COMPLEX)
    print("FIL_FUN = ", FIL_FUN)
    print("DTM_M = ", DTM_M)
    print("DTM_P = ", DTM_P)
    print("NUM_EPOCHS = ", NUM_EPOCHS)
    print("BATCH_SIZE = ", BATCH_SIZE)
    print("TRAIN_SIZES = ", TRAIN_SIZES)

    '''
    this code should:
    1. import all (ellipsoid and rips persistence diagrams in dimension 1 for datasets of given id) into two lists, import labels into a list
    2. divide all those into train and test data
    3.
    '''



    experiment_summaries = load_all_experiment_summaries(experiment_summaries_path)
    summaries_grouped_by_parameters = group_experiment_summaries_by_parameters(experiment_summaries)

    for parameters, summaries in summaries_grouped_by_parameters:
        pde_1 = convert_to_turkevs_barcode_type(summaries)

    # pde has like pde["transformation"][dataset]


    barcodes,labels = extract_barcodes_and_labels(experiment_summaries)




    _, pde1, _, pdr1, points_dict_of_lists, labels = data_handling.import_turkevs_transformed(json_data_folder)
    data_pc_trnsfs, labels_check = get_points_and_labels_from_dict(points_dict_of_lists)
    data_pc = data_pc_trnsfs['std']
    num_point_clouds = len(data_pc)

    datasets = read_turkevs_datasets(datasetspath, id=id)

    # indices = np.arange(num_point_clouds)
    train_indices, test_indices = divide_indices_into_train_and_test(np.arange(num_point_clouds))

    vars_to_save = {}
    vars_to_save['train_indices'] = train_indices
    vars_to_save['test_indices'] = test_indices

    # Labels.
    labels_train = labels[train_indices]
    labels_test = labels[test_indices]
    labels_con, label_encoder = data_construction.encode_labels(labels)
    labels_con_train = labels_con[train_indices]
    labels_con_test = labels_con[test_indices]

    data_pc_test_trnsfs = select_dict_data(data_pc_trnsfs, test_indices)

    accs = {}
    pipelines = []




    # Ellipsoids
    print("\n\nReading in  ellipsoids PDs (input for PH)...")
    data_pde_train = pde1['std'][train_indices]
    data_pde_test_trnsfs = select_dict_data(pde1, test_indices)
    accs["PHE"] = calculate_accuracy_trnsfs(ph_ml, data_pde_train, labels_train, data_pde_test_trnsfs, labels_test, TRAIN_SIZES, results_path=results_path, name='PHE')
    pipelines.append('PHE')

    # Rips
    print("\n\nReading in  ellipsoids PDs (input for PH)...")
    data_pde_train = pdr1['std'][train_indices]
    data_pde_test_trnsfs = select_dict_data(pde1, test_indices)
    accs["PHR"] = calculate_accuracy_trnsfs(ph_ml, data_pde_train, labels_train, data_pde_test_trnsfs, labels_test, TRAIN_SIZES, results_path=results_path, name='PHR')
    pipelines.append('PHR')

    # PDs.
    print("\n\nCalculating PDs (input for PH)...")
    _, data_pd = ph.calculate_pds_point_clouds(data_pc, fil_complex = FIL_COMPLEX, fil_fun = FIL_FUN, m = DTM_M, p = DTM_P)    
    data_pd_train = data_pd[train_indices] 
    data_pd_test = data_pd[test_indices]
    data_pd_test_trnsfs = {}
    for trnsf in trnsfs:
        data_pc_test = data_pc_test_trnsfs[trnsf]
        _, data_pd_test = ph.calculate_pds_point_clouds(data_pc_test, fil_complex = FIL_COMPLEX, fil_fun = FIL_FUN, m = DTM_M, p = DTM_P)
        data_pd_test_trnsfs[trnsf] = data_pd_test
    accs["PH"] = calculate_accuracy_trnsfs(ph_ml, data_pd_train, labels_train, data_pd_test_trnsfs, labels_test, TRAIN_SIZES, results_path=results_path, name='PH')
    pipelines.append('PH')

    # PH simple.
    print("\n\nCalculating sorted lifespans from PDs (input for simple PH)...")
    data_ph = ph.sorted_lifespans_pds(data_pd, size = 10)
    data_ph_train = data_ph[train_indices] 
    data_ph_test_trnsfs = {}
    for trnsf in trnsfs:
        data_pd_test = data_pd_test_trnsfs[trnsf]
        data_ph_test = ph.sorted_lifespans_pds(data_pd_test, size = 10)
        data_ph_test_trnsfs[trnsf] = data_ph_test
    accs["PH simple"] = calculate_accuracy_trnsfs(ml, data_ph_train, labels_train, data_ph_test_trnsfs, labels_test, TRAIN_SIZES, results_path=results_path, name='PH_simple')
    pipelines.append('PH simple')

    # # [comment the next three blocks for speed]
    # Distance matrices.
    print("\n\nCalculating distance matrices (input for ML and NN)...")
    data_dis_mat_flat = data_construction.calculate_distance_matrices_flat(data_pc)
    data_dis_mat_flat_train = data_dis_mat_flat[train_indices]
    data_dis_mat_flat_test_trnsfs = data_construction.calculate_distance_matrices_flat_under_trnsfs(point_clouds_trnsfs=data_pc_test_trnsfs, trnsfs=trnsfs)

    # ML
    accs["ML"] = calculate_accuracy_trnsfs(ml, data_dis_mat_flat_train, labels_train, data_dis_mat_flat_test_trnsfs, labels_test, TRAIN_SIZES, results_path=results_path, name='ML')
    pipelines.append('ML')
    # NN shallow
    accs["NN shallow"] = calculate_accuracy_trnsfs(nn_shallow, data_dis_mat_flat_train, labels_con_train, data_dis_mat_flat_test_trnsfs, labels_con_test, TRAIN_SIZES, results_path=results_path, name='NN_shallow')
    pipelines.append('NN shallow')
    # NN deep
    accs["NN deep"] = calculate_accuracy_trnsfs(nn_deep, data_dis_mat_flat_train, labels_con_train, data_dis_mat_flat_test_trnsfs, labels_con_test, TRAIN_SIZES, results_path=results_path, name='NN_deep')
    pipelines.append('NN deep')

    # PointNet (need 3D point clouds and labels_con_test)
    print("\n\nCalculating 3D point clouds (input for PointNet)...")
    data_pc_3d = data_construction.calculate_3d_point_clouds(data_pc)
    data_pc_3d_train = data_pc_3d[train_indices] 
    data_pc_3d_test_trnsfs = data_construction.calculate_3d_point_clouds_under_trnsfs(point_clouds_trnsfs = data_pc_test_trnsfs, trnsfs = trnsfs)
    accs["PointNet"] = calculate_accuracy_trnsfs(point_net, data_pc_3d_train, labels_con_train, data_pc_3d_test_trnsfs, labels_con_test, TRAIN_SIZES, results_path=results_path, name='PointNet')
    pipelines.append('PointNet')

    vars_to_save['accs'] = accs

    unique_id = id + data_handling.get_timestamp()

    filename_save_vars = os.path.join(path_results,'turkevs_variables_' + unique_id)
    data_handling.save_to_json(vars_to_save, filename=filename_save_vars, add_timestamp=False)

    transformations = ["original", "translation", "rotation", "stretch", "shear", "gaussian", "outliers"]
    fig = plots.plot_bar_chart(transformations, accs, pipelines)
    plt.savefig(os.path.join(path_results, "accs_trnsfs_" + unique_id), bbox_inches = "tight")



if __name__ == '__main__':

    n_runs = 20

    ids = ['id=0005'] # (0001 is the first downsampled, also calculated with the prev version of the code)

    for id in ids:

        experiment_summaries_path = find_turkevs_experiment_summaries_path(id)
        dataset_path = find_turkevs_datasets_path(id)
        results_path = generate_turkevs_results_path(id)

        # experiment_summaries_path

        for i in np.arange(n_runs):
            run_experiments(experiment_summaries_path, dataset_path, results_path)
