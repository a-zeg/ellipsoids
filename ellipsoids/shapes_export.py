import data_handling
from data_handling import savemat
from scipy.io import savemat
from data_handling import figure_eight
import os
import shapes

n_pts = 1000
folder = os.path.join('datasets','shapes')

# circle:
points = data_handling.sample_from_circle(n_pts=n_pts, variation=0)
vars = { 'points': points }
filename = 'circle' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# sphere:
points = data_handling.sample_from_sphere(n_pts=n_pts, ambient_dim = 3, variation=0)
vars = { 'points': points }
filename = 'sphere' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# figure_eight:
points = data_handling.figure_eight(n_pts, 1, 0.1)
vars = { 'points': points }
filename = 'figure_eight' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# annulus:
points = shapes.sample_from_annulus(n=n_pts, r=0.4, R=1, seed=0)
vars = { 'points': points }
filename = 'annulus' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)


