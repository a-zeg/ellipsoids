# from  import get_paths_of_files_in_a_folder
import data_handling
import numpy as np
import os
import json
import sys
import re

def calculate_turkevs(path: str):
    '''
    Read in the turkevs point cloud data (the output of generate_turkevs.py), calculate the ellipsoids barcodes
    and save to files.

    The names of the files include the values of important parameters such as:
    - id - indicates which dataset is used (ids correspond to different runs of generate_turkevs.py)
    - transformation (e.g. 'std' in [...]Turkevs-std-[...])
    - index (e.g. 002 in [...]Turkevs-std-002) (note to self: this is important to get right because the labels depend on the indices)
    '''

    longest_axis = 3    # will be expanded to [3:1] or [3:3:1] depending on ambient dimension below
    expansion_dim = 2
    nbhd_size = 5

    with open(path, 'r') as f:
        data_pc_trnsfs = json.load(f)

    id = re.search(r'id=(\d+)', path).group(1)

    subfolder = os.path.join('data', 'id=' + str(id))
    if not os.path.isdir(subfolder):
        os.makedirs(subfolder)
        print('Created folder ' + subfolder)

    # objects = []
    # with (open(picklepath, "rb")) as openfile:
    #     while True:
    #         try:
    #             objects.append(pickle.load(openfile))
    #         except EOFError:
    #             break
    # data_transformed = objects[0]

    labels = data_pc_trnsfs['labels']

    transformations = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

    for transformation in transformations:
        N = len(data_pc_trnsfs[transformation])
        for i in np.arange(N):

            print(transformation + ': dataset ' + str(i) + ' of ' + str(N))
            points = np.asarray(data_pc_trnsfs[transformation][i])

            print(len(points))
            data_type = 'Turkevs-' + transformation + '-' + str(i).zfill(3)
            ambient_dim = len(points[0])

            if ambient_dim == 2:
                axes_ratios = np.array([longest_axis,1])
            elif ambient_dim == 3:
                axes_ratios = np.array([longest_axis,longest_axis,1])

            # Generate the filename for storing the variables 
            n_pts = len(points)
            data_type_params = {'id': id}
            filename_parameters = data_handling.set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
            save_filename = data_handling.generate_filename(filename_parameters, folder=subfolder) 

            vars_dict = {}
            vars_dict['label'] = labels[i]
            vars_dict['id'] = id
            vars_dict['transformation'] = transformation
            vars_dict['index'] = i

            data_handling.calculate_and_save_ellipsoids_and_rips_data(
                points,
                nbhd_size,
                axes_ratios,
                expansion_dim,
                save_filename,
                vars_dict
            )  



if __name__ == '__main__':

    folder = 'datasets/turkevs'
    # filename = 'pc_test_trnsfs_N=100_n=20_id=0008.json'
    # path = os.path.join(folder,filename)

    path = str(sys.argv[1])
    # path = 'datasets/turkevs/pc_test_trnsfs_N=20_n=20_id=0012.json'

    calculate_turkevs(path)

