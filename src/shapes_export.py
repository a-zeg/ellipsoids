from data_handling import createData
from data_handling import savemat
from scipy.io import savemat
from data_handling import figure_eight
import os
from shapes import sample_from_annulus

n_pts = 1000
folder = os.path.join('datasets','shapes')

# circle:
points = createData(nPoints=n_pts, type='circle', variation=0)
vars = { 'points': points }
filename = 'circle' + '_' + f'{n_pts=}' + '.mat'
savemat(os.path.join(folder,filename), vars)

# sphere:
points = createData(nPoints=n_pts, type='sphere', ambient_dim = 3, variation=0)
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


