import os
from os import listdir
from os.path import isfile, join
import sys
import json
import re
from datetime import datetime

import pickle
import gudhi as gd
import numpy as np
from numpy import genfromtxt
from scipy.io import loadmat
from scipy.io import savemat

sys.path.append(os.path.abspath('.'))

from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.topological_computations import calculate_ellipsoid_barcode
from ellipsoids.topological_computations import calculate_rips_barcode
from ellipsoids.topological_computations import expandTreeAndCalculateBarcode



def sample_from_circle(n_pts=100, variation=0.1, outlier=False):
    if outlier is True: 
        n_pts = n_pts - 1

    r = 1
    t = np.linspace(0, 2*np.pi * (n_pts-1)/n_pts, n_pts)
    x = r*np.cos(t) + variation * np.random.rand(n_pts)
    y = r*np.sin(t) + variation * np.random.rand(n_pts)
    output = np.vstack((x,y)).transpose()

    if outlier is True:
        output = np.append(output,[[0,0]],axis=0)
    
    return output



def sample_from_ellipse(n_pts=100, a=2, b=1, variation=0.1):
    t = np.linspace(0, 2*np.pi * (n_pts-1)/n_pts, n_pts)
    x = a * np.cos(t) + variation * np.random.rand(n_pts)
    y = b * np.sin(t) + variation * np.random.rand(n_pts)
    return np.vstack((x,y)).transpose()



def sample_from_cassini_oval(n_pts=100, variation=0.1):
    t = np.linspace(-1, 1, int(n_pts/2))
    #t = np.sign(t)*np.abs(t)**(1/4)
    x = np.concatenate((t,t)) + variation * np.random.rand(n_pts)
    yh = (t**2 + 0.5) * np.sqrt(1 - t**2)
    y = np.concatenate((-yh, yh)) + variation * np.random.rand(n_pts)
    
    return np.vstack((x,y)).transpose()



def sample_from_sphere(n_pts=100, ambient_dim=3, r=1):
    vec = np.random.randn(ambient_dim, n_pts)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.transpose()



def sample_from_torus(n_pts=100, R=2, r=1):
    nPtsSampled = 0
    theta = np.zeros([n_pts])
    phi = np.zeros([n_pts])

    # rejection sampling
    while nPtsSampled < n_pts:
        thetaSample = 2 * np.pi * np.random.rand()
        phiSample = 2 * np.pi * np.random.rand()
        W = 2 * np.pi * np.random.rand()

        if W <= (R + r * np.cos(thetaSample))/(R+r):
            theta[nPtsSampled] = thetaSample
            phi[nPtsSampled] = phiSample
            nPtsSampled += 1


        x = (R + r * np.cos(theta)) * np.cos(phi)
        y = (R + r * np.cos(theta)) * np.sin(phi)
        z = r * np.sin(theta)

    return np.vstack((x,y,z)).transpose() 



def figure_eight(n, a, b, variation=0):
    # adapted from Bastian Rieck
    """Sample a set of points from a figure eight curve.

    Parameters
    ----------
    n : int
        Number of points to sample

    a : float
        Controls extents of the curve. A larger `a` parameter will
        result in larger scaling.

    b : float
        Controls neck size of the curve. A larger `b` parameter will
        result in an increased neck size.

    Returns
    -------
    np.array
        Array of shape (n, 2). Will contain the sampled points.
    """
    T = np.linspace(-3.14, 3.14, num=n)

    X = a * np.sin(T)
    Y = a * np.sin(T)**2 * np.cos(T) + b * np.cos(T)

    X = np.column_stack((X, Y))
    X += np.random.default_rng().uniform(0, variation, size=(n, 2))
    return X



def sample_from_annulus(n, r, R, seed=None):
    # taken from Bastian Rieck
    """Sample points from a 2D annulus.

    This function samples `N` points from an annulus with inner radius `r`
    and outer radius `R`.

    Parameters
    ----------
    n : int
        Number of points to sample

    r : float
        Inner radius of annulus

    R : float
        Outer radius of annulus

    seed : int, instance of `np.random.Generator`, or `None`
        Seed for the random number generator, or an instance of such
        a generator. If set to `None`, the default random number
        generator will be used.

    Returns
    -------
    torch.tensor of shape `(n, 2)`
        Tensor containing sampled coordinates.
    """
    if r >= R:
        raise RuntimeError(
            'Inner radius must be less than or equal to outer radius'
        )

    rng = np.random.default_rng(seed)
    thetas = rng.uniform(0, 2 * np.pi, n)

    # Need to sample based on squared radii to account for density
    # differences.
    radii = np.sqrt(rng.uniform(r ** 2, R ** 2, n))

    X = np.column_stack((radii * np.cos(thetas), radii * np.sin(thetas)))
    return X




def import_maxmin_mat(path: str):
    return np.asarray(loadmat(path)['maxmin_points']) # size: 6040 x 24 



def saveMAT(variable, filename, varname='points'):
    pointsDict = {varname: variable}
    savemat(filename, pointsDict)



def save_pentagon_points():
    filenameSaved = 'pentagonsamplesSmall2.mat'

    points = genfromtxt('pentagonsamplesSmall2.txt', delimiter=',')
    if np.isnan(points).any():
        print('The dataset contains NaN values, which will be removed.')
        points = points[~np.isnan(points).any(axis=1)]
    
    saveMAT(points, filenameSaved, varname='pentagonsamples')



def get_timestamp():
    return datetime.now().strftime("_%Y%m%d_%H%M%S")



class CustomEncoder(json.JSONEncoder):
    def default(self, obj):

        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
        elif isinstance(obj, Ellipsoid):
            obj_data = {
                "center": obj.center,
                "axes": obj.axes.tolist(),
                "axesLengths": obj.axesLengths.tolist()
            }
            return obj_data
        
        elif isinstance(obj,gd.SimplexTree):
            return list(obj.get_filtration())
        
        elif isinstance(obj, np.integer):
            return int(obj)
        
        return json.JSONEncoder.default(self, obj)



def save_variables(
        dictOfVars,
        filename=datetime.now().strftime("data/test.json"), 
        timestamp=True):
    
    print('Saving data to file...')

    if timestamp:
        filename = filename + '_' + get_timestamp()

    if not filename.endswith('.json'):
        filename = filename + '.json'

    json_string = json.dumps(dictOfVars, cls=CustomEncoder, indent=4)
    with open(filename, 'w') as outfile:
        outfile.write(json_string)
    print("Data saved to file " + filename + '.')

    return filename



def continuously_save_variables(
        dictOfVars,
        filename=datetime.now().strftime("data/test.json"), 
        addToStart=False):
    
    print('Saving data to file...')
    json_string = json.dumps(dictOfVars, cls=CustomEncoder, indent=4)
    if os.path.exists(filename): # if the file already exists, append to it
        if addToStart:
            index = 0
            with open(filename, "r") as f:
                contents = f.readlines()
            contents.insert(index, dictOfVars)
            with open(filename, "w") as f:
                contents = "".join(str(c) for c in contents)
                f.write(contents)
        else:
            with open(filename, "a") as f:
                f.write(json_string)
    else:
        with open(filename, 'w') as outfile:
            outfile.write(json_string)
    print("Data saved to file " + filename + '.')



def read_variable_from_dict(var_name, dictionary):
   if var_name in dictionary:
        vars[var_name] = dictionary[var_name] 



def read_variables(filename):
    with open(filename, "r") as f:
        json_vars = json.load(f)
    
    return json_process_variables(json_vars)



def _json_process_ellipsoid_list(ellipsoid_list_raw):

    ellipsoidList = []
    for ellipsoid in ellipsoid_list_raw:
        ellipsoidList.append(Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                    np.asarray(ellipsoid['axesLengths'])))
    
    return ellipsoidList



def _json_process_simplex_tree(simplex_tree_raw):

    simplex_tree = gd.SimplexTree()
    for simplex_tree_entry in simplex_tree_raw:
        simplex_tree.insert(simplex_tree_entry[0],simplex_tree_entry[1])
   
    return simplex_tree



def json_process_variables(json_vars: dict):

    vars = json_vars

    if 'rValues' in json_vars:
        vars['rValues'] = np.asarray(json_vars['rValues'])

    for name in ['ellipse_list', 'ellipsoid_list', 'ellipseList']: # TODO: for backwards compatibility, maybe can remove it
        if name in json_vars:
            vars['ellipsoid_list'] = _json_process_ellipsoid_list(json_vars[name])
            break
        # TODO delete the wrong name here.

    for name in ['simplex_tree_ellipsoids', 'simplexTreeEllipsoids']: # TODO: for backwards compatibility, maybe can remove it
        if name in json_vars:
            vars['simplex_tree_ellipsoids'] = _json_process_simplex_tree(json_vars[name])
            break

    for name in ['simplex_tree_rips', 'simplexTreeRips']: # TODO: for backwards compatibility, maybe can remove it
        if name in json_vars:
            vars['simplex_tree_rips'] = _json_process_simplex_tree(json_vars[name])
            break

    return vars
     

     

def readOFF(filename):
    file = open(filename, 'r')
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    nVerts, nFaces, nEdges = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')] for i_vert in range(nVerts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for i_face in range(nFaces)]
    return np.asarray(verts)



def printListOfSimplices(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    for splx in simplexList:
        print(splx)



def get_paths_of_files_in_a_folder(folder: str, extension='.mat'):
    filenames = [f for f in listdir(folder) if isfile(join(folder, f)) if f.endswith(extension)]
    paths = [os.path.join(folder, f) for f in filenames] 
    return paths



def get_paths_with_seed(paths, seed: int):

    output_paths = []

    for path in paths:
        match = re.search(r'seed=(\d+)', path)
        if match:
            path_seed = match.group(1)
            if str(seed) == path_seed:
                output_paths.append(path)
        else:
            ValueError('No match found.')

    return output_paths



def remove_dim_from_barcode(barcode):

    result_barcode = []
    for bar in barcode:
            result_barcode.append(bar[1])
    return result_barcode



def generate_filename(filename_parameters: dict, folder='data', timestamp = ''):
    '''
    Generates filename by creating a string from the variables in 
    filename_parameters and appends the timestamp.
    '''
    filename = 'ellipsoids'
    for key, value in filename_parameters.items():
        filename = "%s_%s=%s" % (filename, key, value)

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



def calculate_and_save_ellipsoids_and_rips_data(points, nbhd_size, axes_ratios, expansion_dim, filename, additional_vars_dict={}):
    
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
    barcode_ellipsoids, simplex_tree_ellipsoids, ellipsoid_list, t_ellipsoids \
        = calculate_ellipsoid_barcode(points, nbhd_size, axes_ratios, expansion_dim=expansion_dim)
    barcode_rips, simplex_tree_rips, t_rips \
        = calculate_rips_barcode(points, expansion_dim=expansion_dim)

    # Get the execution time
    t_ellipsoids_over_t_rips = t_ellipsoids / t_rips
    t_total = t_ellipsoids + t_rips
    print('\nThe total execution time is ' + str(t_total) + '\n')

    # Save variables to file
    params_dict = filter_dictionary(vars_to_save, locals())

    if additional_vars_dict: # if not empty
        for key, value in additional_vars_dict.items():
            params_dict[key] = value
    
    save_variables(params_dict, filename=filename)



def recalculateBarcodesFromFile(filename, expansionDim=2, collapseEdges=False):
    print('Reading in the variables... ', end='', flush=True)
    vars = read_variables(filename)
    if 'expansionDim' in vars and vars['expansionDim'] == expansionDim:
        print('The original barcode is already expanded to the specified dimension.')
        return None
    if 'simplexTreeEllipsoids' in vars:
        simplexTreeEllipsoids = vars['simplexTreeEllipsoids']
    else: simplexTreeEllipsoids = gd.SimplexTree()
    if 'simplexTreeRips' in vars:
        simplexTreeRips = vars['simplexTreeRips']
    else: simplexTreeRips = gd.SimplexTree()
    print('Done.')
    
    barcodeEllipsoids = expandTreeAndCalculateBarcode(simplexTreeEllipsoids, expansionDim, collapseEdges=collapseEdges)
    barcodeRips = expandTreeAndCalculateBarcode(simplexTreeRips, expansionDim, collapseEdges=collapseEdges)
    
    dictOfVars = {
        'originalFile': filename, 
        'expansionDim': expansionDim,
        'barcodeEllipsoids': barcodeEllipsoids,
        'barcodeRips': barcodeRips
    }

    filename =  filename[:filename.rfind('.')] + '-barcodes_expansionDim=' + f'{expansionDim}' + datetime.now().strftime("_%Y%m%d_%H%M%S") + '.json'
    save_variables(dictOfVars, filename=filename)




############################################
######   turkevs-specific functions   ######
############################################
# The rest of this file consists of functions specific to 
# handling the turkevs data.

# from Turkevs:
# Transform list of PDs with different number of cycles into an array of PDs with the same number of cycles.
def extend_pds_to_length(pds, length):
    pds_ext = np.zeros((len(pds), length, 2))
    for s, pd in enumerate(pds):
        for i in range(length):
            if i < len(pd):
                pds_ext[s][i] = pds[s][i]
            else:
                pds_ext[s][i] = np.asarray([0, 0]) 
    return pds_ext


def import_turkevs_data(picklepath):
    objects = []
    with (open(picklepath, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    
    transformations = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

    data_dict = objects[0]
    points_all = []
    labels_all = []

    for transformation in transformations:
        i = 0

        for data in data_dict[transformation]:
            points_all.append(data)
            labels_all.append('Turkevs-' + str(i).zfill(3) + transformation)
            i = i + 1

    return points_all, labels_all
 


def parse_turkevs_filename(path):
    '''
    Parses a filename of a JSON file generated by 'calculate_turkevs' and returns the 
    name of the transformation and the index corresponding to the given data.
    '''
    transformation = ''
    index = 0

    match = re.search(r'Turkevs-(\w+)-\d+', path)
    if match:
        transformation = match.group(1)
    else:
        ValueError('No match found.')

    match = re.search(r'Turkevs-\w+-(\d+)', path)
    if match:
        index = match.group(1)
        index = int(index)
    else:
        ValueError('No match found.')    

    return transformation, index



# def extract_barcodes_in_dim0_and_dim1(filename):

#     pd0 = []
#     pd1 = []
#     print('Reading in the variables... ', end='', flush=True)
#     vars = read_variables(filename)
#     print('Done.')
#     barcodeEllipsoids = vars['barcodeEllipsoids']

#     for bar in barcodeEllipsoids:
#         if bar[0] == 0:
#             pd0.append(bar[1])
#         elif bar[0] == 1:
#             pd1.append(bar[1])

#     return pd0, pd1


def extract_barcodes_in_dim0_and_dim1(barcode):

    pd0 = []
    pd1 = []

    for bar in barcode:
        if bar[0] == 0:
            pd0.append(bar[1])
        elif bar[0] == 1:
            pd1.append(bar[1])

    return pd0, pd1



def read_pd0_and_pd1(path):

    with open(path, "r") as f:
            vars = json.load(f)
    
    pdE0 = []
    pdE1 = []
    pdR0 = []
    pdR1 = []
    label = None

    if 'barcodeEllipsoids' in vars:
        barcode_ellipsoids = vars['barcodeEllipsoids']
    elif 'barcode_ellipsoids' in vars:
        barcode_ellipsoids = vars['barcode_ellipsoids']
    else:
        exit('The file ' + path + ' does not contain the ellipsoids barcode.')

    if 'barcodeRips' in vars:
        barcode_rips = vars['barcodeRips']
    elif 'barcode_rips' in vars:
        barcode_rips = vars['barcode_rips']
    else:
        exit('The file ' + path + ' does not contain the Rips barcode.')

    if 'points' in vars:
        points = vars['points']
    if 'label' in vars:
        label = int(vars['label'])

    for bar in barcode_ellipsoids:
        if bar[0] == 0:
            pdE0.append(bar[1])
        elif bar[0] == 1:
            pdE1.append(bar[1])

    for bar in barcode_rips:
        if bar[0] == 0:
            pdR0.append(bar[1])
        elif bar[0] == 1:
            pdR1.append(bar[1])

    return pdE0, pdE1, pdR0, pdR1, points, label



def import_turkevs_transformed(folder):
    '''
    Returns three dictionaries: pd0, pd1, and points.
    The keys of each dictionary are the transformations and the values are lists 
    ordered by the mesh numbers.
    
    For example, pd0 is a dictionary, where pd0['std'][0] is a list of 0-dimensional
    bars corresponding to the 0th standard mesh.
    Similarly, points is a dictionary, where points['stretch'][24] is a list of points
    corresponding to the 24th mesh after the stretch transformation has been applied to it.
    '''

    extension = '.json'
    paths = get_paths_of_files_in_a_folder(folder, extension=extension)
    print(str(len(paths)) + ' ' + extension + ' files found in the folder ' + folder + '.')
    if len(paths) == 0: 
        exit('Folder ' + folder + ' is empty.')
    transformations = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

    # ---------- initialising pds0 and pds1 -----------
    # (necessary because the data might be read in a random order and 
    # it needs to be placed at the correct index)

    number_of_files = len(paths)
    number_of_transformations = len(transformations)

    if number_of_files % number_of_transformations != 0:
        exit('Error: the number of files in ' + folder + ' is not consistent with \
              the number of transformations.')

    files_per_transformation = int(number_of_files / number_of_transformations)

    pdsE0 = {}
    pdsE1 = {}
    pdsR0 = {}
    pdsR1 = {}
    points = {}
    labels = [None] * files_per_transformation

    for transformation in transformations:
        pdsE0[transformation] = [None]*files_per_transformation
        pdsE1[transformation] = [None]*files_per_transformation
        pdsR0[transformation] = [None]*files_per_transformation
        pdsR1[transformation] = [None]*files_per_transformation
        points[transformation] = [None]*files_per_transformation
    # --------------------------------------------------

    for path in paths:
        transformation, index = parse_turkevs_filename(path)

        pdE0, pdE1, pdR0, pdR1, points_, label = read_pd0_and_pd1(path)
        pdsE0[transformation][index] = pdE0
        pdsE1[transformation][index] = pdE1
        pdsR0[transformation][index] = pdR0
        pdsR1[transformation][index] = pdR1
        points[transformation][index] = points_
        labels[index] = label


    # -------- Pad to max length ----------
    max_pdE0_length = []
    max_pdE1_length = []

    for transformation in transformations:

        # Transform list of 0-dim PDs with different number of cycles into an array of PDs with the same number of cycles. 
        pdsE0_length = [len(pd) for pd in pdsE0[transformation]]
        max_pdE0_length.append(max(pdsE0_length))

        # Transform list of 1-dim PDs with different number of cycles into an array of PDs with the same number of cycles.    
        pdsE1_length = [len(pd) for pd in pdsE1[transformation]]
        max_pdE1_length.append(max(pdsE1_length)) 
        
    for transformation in transformations:
        pdsE0[transformation] = extend_pds_to_length(pdsE0[transformation], max(max_pdE0_length))
        pdsE1[transformation] = extend_pds_to_length(pdsE1[transformation], max(max_pdE1_length))
    # ------------------------------------
    # -------- Pad to max length ----------
    max_pdR0_length = []
    max_pdR1_length = []

    for transformation in transformations:

        # Transform list of 0-dim PDs with different number of cycles into an array of PDs with the same number of cycles. 
        pdsR0_length = [len(pd) for pd in pdsR0[transformation]]
        max_pdR0_length.append(max(pdsR0_length))

        # Transform list of 1-dim PDs with different number of cycles into an array of PDs with the same number of cycles.    
        pdsR1_length = [len(pd) for pd in pdsR1[transformation]]
        max_pdR1_length.append(max(pdsR1_length)) 
        
    for transformation in transformations:
        pdsR0[transformation] = extend_pds_to_length(pdsR0[transformation], max(max_pdR0_length))
        pdsR1[transformation] = extend_pds_to_length(pdsR1[transformation], max(max_pdR1_length))
    # ------------------------------------

    labels = np.asarray(labels)

    return pdsE0, pdsE1, pdsR0, pdsR1, points, labels



def find_subfolder_with_given_id(parentfolder, id):

    subfolders = os.listdir(parentfolder)

    for subfolder in subfolders:
        print(subfolder)
        if id in subfolder:
            jsondatafolder = os.path.join(parentfolder, subfolder)

    try:
        jsondatafolder
    except NameError:
        print('No folder with ' + id + ' found.')
        exit()
    
    return jsondatafolder


