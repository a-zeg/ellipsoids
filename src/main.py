import numpy as np
import time
from datetime import datetime
from topological_computations import calculate_ellipsoid_barcode
from topological_computations import calculate_rips_barcode
import data_handling 
# from ripsSimplexTree import generateRipsSimplexTree
# from utils import expandTreeAndCalculateBarcode
# from utils import maxFiltration
# from utils import padAxesRatios
from os import listdir
from os.path import isfile, join
import os
# import data_construction




def generate_filename(filename_parameters: dict, folder='data', timestamp = ''):
    '''
    Generates filename by creating a string from the variables in 
    filename_parameters and appends the timestamp if the value of 
    `timestamp` is True.
    '''
    filename = 'ellipsoids'
    for key, value in filename_parameters.items():
        filename = filename + '_' + key + '=' + str(value)

    if timestamp != '':
        filename = filename + '_' + timestamp

    return os.path.join(folder, filename)

def set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params: dict):

    filename_params = {
        'data_type': data_type,
        'n_pts': n_pts,
        'nbhd_size': nbhd_size,
        'axes_ratios': axes_ratios,
    }

    filename_params.update(data_type_params)

    return filename_params


def filter_dictionary(vars_to_save: list[str], dict_all_vars):

    results = {}
    for var_name in vars_to_save:
        if var_name in dict_all_vars:
            results[var_name] = dict_all_vars[var_name]
        else: 
            print('Warning: ' + var_name + ' does not exist in the local variables and will not be saved.')
    
    return results



def calculate_and_save_ellipsoids_and_rips_data(points, nbhd_size, axes_ratios, expansion_dim, filename):
    
    # Specify the names of variables to be saved
    vars_to_save = [
        'ambient_dim',
        'expansion_dim',
        'nbhd_size',
        'n_pts',
        't_total',
        't_ellipsoids_over_t_rips',
        'points',
        'barcode_ellipsoids',
        'barcode_rips'
    ]

    if 'ambient_dim' in vars_to_save:
        ambient_dim = len(points[0])
    if 'n_pts' in vars_to_save:
        n_pts = len(points)

    # Calculate barcodes for ellipsoids and Rips complexes
    barcode_ellipsoids, simplex_tree_ellipsoids, ellipsoid_list, t_ellipsoids = calculate_ellipsoid_barcode(points, nbhd_size, axes_ratios, expansion_dim=expansion_dim)
    barcode_rips, simplex_tree_rips, t_rips = calculate_rips_barcode(points, expansion_dim=expansion_dim)

    # Get the execution time
    t_ellipsoids_over_t_rips = t_ellipsoids / t_rips
    t_total = t_ellipsoids + t_rips
    print('\nThe total execution time is ' + str(t_total) + '\n')

    # Save variables to file
    params_dict = filter_dictionary(vars_to_save, locals())
    data_handling.saveVarsToFile(params_dict, filename=filename)
        

# def get_paths_of_files_in_a_folder(folder, extension='.mat'):
#     filenames = [f for f in listdir(folder) if isfile(join(folder, f)) if f.endswith(extension)]
#     paths = [(folder + '/' + f) for f in filenames] 
#     return paths
    

def main():

    # data_type = 'cyclooctaneMaxmin'
    # axes_ratios_all = np.array([[3,3,1]])
    # expansion_dim = 3
    # nbhd_sizes = [26, 30, 40]
    # folder = 'datasets/cyclooctane/maxmin_test'
    # paths = get_paths_of_files_in_a_folder(folder, extension='.mat')

    # for path in paths:
    #     for axes_ratios in axes_ratios_all:
    #         for nbhd_size in nbhd_sizes:

    #             # Import points
    #             points, data_type_params = data_handling.importPoints(path=path)

    #             # Generate the filename for storing the variables 
    #             n_pts = len(points)
    #             filename_parameters = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
    #             save_filename = generate_filename(filename_parameters)

    #             calculate_and_save_ellipsoids_and_rips_data(points, nbhd_size, axes_ratios, expansion_dim, save_filename)



    # picklepath = 'datasets/turkevs2022on-main/DATASETS/holes/point_clouds_300.pkl'
    picklepath = 'datasets/turkevs2022on-main/DATASETS/holes/point_clouds_test_trnsfs.pkl' 

    objects = []
    with (open(picklepath, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    
    transformations = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

    
    points_all = objects[0]

    # print(len(points_all['std']))
    transformation = transformations[0]
    print(len(points_all[transformation]))
    print(len(points_all[transformation][0]))


    

    # N = len(objects[0])

    # start = int(sys.argv[1])
    # N = int(sys.argv[2])

    # print(type(objects[0]))

    data_pc_transformed = objects[0]

    
    exit()

    start = 0
    N = len(points_all[transformation[0]])

    print(start)
    print(N)

    for transformation in transformations:
        for i in np.arange(N):

            print('Dataset ' + str(i) + ' of ' + str(N))
            points = objects[0][i]

            data_type = 'Turkevs-' + transformation + '-' + str(i).zfill(3)

            axes_ratios = np.array([3,1])
            expansion_dim = 2
            seed = 0
            nbhd_size = 9
            ambient_dim = len(points[0])

            if ambient_dim == 2:
                axes_ratios = np.array([3,1])
            elif ambient_dim == 3:
                axes_ratios = np.array([3,3,1])

            # Generate the filename for storing the variables 
            n_pts = len(points)
            data_type_params = {}
            filename_parameters = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
            save_filename = generate_filename(filename_parameters)

            calculate_and_save_ellipsoids_and_rips_data(
                points, 
                nbhd_size, 
                axes_ratios,
                expansion_dim,
                save_filename
            )  
        
import pickle 
import sys

def calculatePDforTurkevs():
    picklepath = 'datasets/turkevs2022on-main/DATASETS/holes/point_clouds_300.pkl'

    objects = []
    with (open(picklepath, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break

    # N = len(objects[0])

    start = int(sys.argv[1])
    N = int(sys.argv[2])

    print(start)
    print(N)

    for i in np.arange(start, start+N):

        print('Dataset ' + str(i) + ' of ' + str(N))
        points = objects[0][i]

        dataType = 'Turkevs-' + str(i).zfill(3)

        axesRatios = np.array([3,1])
        expansionDim = 2
        seed = 0
        nbhdSize = 9
        dim = len(points[0])

        if dim == 2:
            axesRatios = np.array([3,1])
        elif dim == 3:
            axesRatios = np.array([3,3,1])

        main(
            nPts=len(points), 
            dataType=dataType, 
            collapseEdges=True, 
            axesRatios=axesRatios, 
            expansionDim=expansionDim, 
            dim=dim,
            seed=int(seed),
            nbhdSize=nbhdSize,
            drawPoints=False,
            drawEllipsoidsSimplexTree=False,
            drawEllipsoids = False, 
            rPlot=0.3,
            points=points,
            saveSimplexTree = False,
            savePoints = False
        )  


if __name__ == "__main__":

    main()