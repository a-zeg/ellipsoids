import os
import sys

sys.path.append(os.path.abspath('.'))

from scipy.io import savemat

from ellipsoids.data_handling import savemat
from ellipsoids.data_handling import sample_from_circle
from ellipsoids.data_handling import sample_from_sphere
from ellipsoids.data_handling import figure_eight
from ellipsoids.data_handling import sample_from_annulus



n_pts = 1000
folder = os.path.join('datasets','shapes')

# circle:
points = sample_from_circle(n_pts=n_pts, variation=0)
vars = { 'points': points }
filename = 'circle' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# sphere:
points = sample_from_sphere(n_pts=n_pts, ambient_dim = 3, variation=0)
vars = { 'points': points }
filename = 'sphere' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# figure_eight:
points = figure_eight(n_pts, 1, 0.1)
vars = { 'points': points }
filename = 'figure_eight' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# annulus:
points = sample_from_annulus(n=n_pts, r=0.4, R=1, seed=0)
vars = { 'points': points }
filename = 'annulus' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)


