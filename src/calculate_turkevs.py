# from  import get_paths_of_files_in_a_folder
from main import set_filename_parameters
from main import generate_filename
from main import calculate_and_save_ellipsoids_and_rips_data
import numpy as np
import pickle
import re
import os

def calculate_turkevs(picklepath: str):
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
    nbhd_size = 9

    objects = []
    with (open(picklepath, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break

    id = re.search(r'id=(\d+)', picklepath).group(1)

    transformations = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]
    data_transformed = objects[0]

    for transformation in transformations:
        N = len(data_transformed[transformation])
        for i in np.arange(N):

            print(transformation + ': dataset ' + str(i) + ' of ' + str(N))
            points = data_transformed[transformation][i]

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
            filename_parameters = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)
            save_filename = generate_filename(filename_parameters) 

            calculate_and_save_ellipsoids_and_rips_data(
                points,
                nbhd_size,
                axes_ratios,
                expansion_dim,
                save_filename
            )  


if __name__ == '__main__':

    folder = 'datasets/turkevs'
    filename = 'pc_test_trnsfs_N=100_n=100_id=0006.pkl'
    path = os.path.join(folder,filename)
    calculate_turkevs(path)