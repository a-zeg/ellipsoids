import numpy as np
import time
import sys
import os

sys.path.append(os.path.abspath('.'))

from ellipsoids.turkevs.ph import calculate_pds
from ellipsoids.topological_computations import calculate_ellipsoid_barcode
from ellipsoids.topological_computations import calculate_rips_barcode
from ellipsoids.data_handling import get_timestamp
from ellipsoids.turkevs.data_construction import sample_point_cloud


def turkevs_runtime_analysis():

    ########################################
    output_folder = 'data/runtime_analysis'
    shape = 'circle'
    n_turkevs = 1000
    n_ellipsoids = 301
    ########################################



    FIL_COMPLEX = "alpha"
    FIL_FUN = "dtm"
    DTM_M = 0.03
    DTM_P = 1
    point_cloud_turkevs = sample_point_cloud(1000, shape) # should be 1000 points
    t_start_turkevs = time.time()
    turkevs_pd0, turkevs_pd1 = calculate_pds(point_cloud_turkevs, fil_complex = FIL_COMPLEX, fil_fun = FIL_FUN, m = DTM_M, p = DTM_P)
    t_end_turkevs = time.time()
    t_turkevs = t_end_turkevs - t_start_turkevs


    point_cloud_ellipsoids = sample_point_cloud(301, shape) # should be 301 points
    nbhd_size = 9
    axes_ratios = np.array([3,1,1])
    t_start_ellipsoids = time.time()
    ellipsoid_barcode = calculate_ellipsoid_barcode(point_cloud_ellipsoids, nbhd_size, axes_ratios)
    t_end_ellipsoids = time.time()
    t_ellipsoids = t_end_ellipsoids - t_start_ellipsoids


    point_cloud_rips = sample_point_cloud(1000, shape) # should be 1000 points
    t_start_rips = time.time()
    rips_barcode = calculate_rips_barcode(point_cloud_rips, expansion_dim=2)
    t_end_rips = time.time()
    t_rips = t_end_rips - t_start_rips



    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    output_filename = 'runtime_analysis_%s' %get_timestamp()
    output_path = os.path.join(output_folder, output_filename)

    f = open(output_path, "w")
    f.write("Runtime analysis %s\n\n" %get_timestamp())

    f.write("turkevs\n")
    f.write("n_pts: %s\n" %len(point_cloud_turkevs))
    f.write("t_turkevs: %s\n\n" %t_turkevs)

    f.write("ellipsoids\n")
    f.write("n_pts: %s\n" %len(point_cloud_ellipsoids))
    f.write("t_ellipsoids: %s\n" %t_ellipsoids)

    f.write("rips\n")
    f.write("n_pts: %s\n" %len(point_cloud_rips))
    f.write("t_rips: %s\n" %t_rips)
    f.close()



if __name__ == '__main__':
    turkevs_runtime_analysis()



