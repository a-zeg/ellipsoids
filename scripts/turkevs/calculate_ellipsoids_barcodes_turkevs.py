import numpy as np
import os
import json
import sys
import re
from dataclasses import dataclass

sys.path.append(os.path.abspath('.'))

from ellipsoids.data_handling import set_filename_parameters
from ellipsoids.data_handling import generate_filename
from ellipsoids.data_handling import calculate_and_save_ellipsoids_and_rips_data

from ellipsoids.common import Parameters
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import Dataset
from ellipsoids.common import ComplexType
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import TurkevsParameters
from ellipsoids.common import Experiment


def get_turkevs_dataset_id(path: str):
    return re.search(r'id=(\d+)', path).group(1)


@dataclass
class TurkevsDatasetParameters():
    label: str
    id: str
    transformation: str
    dataset_index: int



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
    params_ellipsoid_ripstype = EllipsoidParameters(nbhd_size=5, axes_ratios=np.array([3,1]), complex_subtype=ComplexSubtype.RIPS, save_simplex_tree=False, save_ellipsoid_list=False)
    params_ellipsoid_alpha = EllipsoidParameters(nbhd_size=5, axes_ratios=np.array([3,1]), complex_subtype=ComplexSubtype.ALPHA, save_simplex_tree=False, save_ellipsoid_list=False)
    params_rips = Parameters(complex_subtype=ComplexSubtype.RIPS, save_simplex_tree=False)
    params_alpha = Parameters(complex_subtype=ComplexSubtype.ALPHA, save_simplex_tree=False)

    with open(path, 'r') as f:
        data_pc_trnsfs = json.load(f)

    id = get_turkevs_dataset_id(path)

    save_location = os.path.join('data', f'test_20250314_id={id}')

    # subfolder = os.path.join(save_location, f'id={id}')
    # if not os.path.isdir(subfolder):
    #     os.makedirs(subfolder)
    #     print('Created folder ' + subfolder)

    labels = data_pc_trnsfs['labels']

    transformations = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

    for transformation in transformations:
        N = len(data_pc_trnsfs[transformation])
        for i in np.arange(N):

            print(transformation + ': dataset ' + str(i) + ' of ' + str(N))
            points = np.asarray(data_pc_trnsfs[transformation][i])
            data_type = f'turkevs-{id}-{transformation}-{str(i).zfill(3)}'
            dataset = Dataset(points=points, data_type=data_type)

            ambient_dim = dataset.ambient_dim()

            if ambient_dim == 2:
                axes_ratios = np.array([longest_axis,1])
                for ellipsoid_params in [params_ellipsoid_ripstype, params_ellipsoid_alpha]:
                    ellipsoid_params.axes_ratios = axes_ratios
            elif ambient_dim == 3:
                axes_ratios = np.array([longest_axis,longest_axis,1])
                for ellipsoid_params in [params_ellipsoid_ripstype, params_ellipsoid_alpha]:
                    ellipsoid_params.axes_ratios = axes_ratios

            turkevs_params = TurkevsParameters(label = str(labels[i]), id=id, transformation=transformation, dataset_index=i)

            experiments = []

            for params in ([params_rips, params_alpha, params_ellipsoid_ripstype, params_ellipsoid_alpha]):
                experiment = Experiment(dataset, params)
                experiment.turkevs_parameters = turkevs_params
                experiments.append(experiment)

            for experiment in experiments:
                experiment.run()
                experiment.save_to_json()
                experiment.print_execution_time()

            # plot_experiments(experiments)



if __name__ == '__main__':

    folder = 'datasets/turkevs'
    # filename = 'pc_test_trnsfs_N=100_n=20_id=0008.json'
    # path = os.path.join(folder,filename)

    # path = str(sys.argv[1])
    path = 'datasets/turkevs/pc_test_trnsfs_N=20_n=20_id=0012.json'

    calculate_turkevs(path)

