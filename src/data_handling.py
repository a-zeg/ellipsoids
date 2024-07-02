import numpy as np
from numpy import genfromtxt
# import figure_eight
from scipy.io import loadmat
from scipy.io import savemat
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
import random
from os import listdir
from os.path import isfile, join
import os
# import open3d
import json
import gudhi as gd
from datetime import datetime
import pickle
import re

import shapes
import topological_computations

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

def createData(nPoints, type, variation = 0.1, ambient_dim=2, outliers=False):
    ''' Generate point cloud data of type 'type' and consiting of 'nPoints' points.
    :nPoints: Number of points to generate
    :type: Type of data (e.g. 'circle', 'ellipse', 'Cassini_oval')
    :variation: Amount of variation from the :type:
    :dim: Not used yet
    :return: nPoints x dim numpy array
    '''
    np.random.seed(0)

    if type == 'circle':
        if outliers is True: nPoints = nPoints - 1
        r = 1
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        x = r*np.cos(t) + variation * np.random.rand(nPoints)
        y = r*np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

        if outliers is True:
            output = np.append(output,[[0,0]],axis=0)

    elif type == 'ellipse':
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        a = 2
        b = 1
        x = a * np.cos(t) + variation * np.random.rand(nPoints)
        y = b * np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

    elif type == 'Cassini_oval':
        t = np.linspace(-1, 1, int(nPoints/2))
        #t = np.sign(t)*np.abs(t)**(1/4)
        x = np.concatenate((t,t)) + variation * np.random.rand(nPoints)
        yh = (t**2 + 0.5) * np.sqrt(1 - t**2)
        y = np.concatenate((-yh, yh)) + variation * np.random.rand(nPoints)
        
        output = np.vstack((x,y)).transpose()

    elif type == 'sphere':
        r = 1

        vec = np.random.randn(ambient_dim, nPoints)
        vec /= np.linalg.norm(vec, axis=0)
        output = vec.transpose()

    elif type == 'torus':
        r = 1
        R = 2

        nPtsSampled = 0
        theta = np.zeros([nPoints])
        phi = np.zeros([nPoints])

        # rejection sampling
        while nPtsSampled < nPoints:
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

        output = np.vstack((x,y,z)).transpose()

    else:
        raise Exception(type + " is an invalid type of data.")
    
    return output

def figure_eight(n, a, b):
    # Taken from Bastian Rieck
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
    T = np.linspace(-5, 5, num=n)

    X = a * np.sin(T)
    Y = a * np.sin(T)**2 * np.cos(T) + b * np.cos(T)

    X = np.column_stack((X, Y))
    X += np.random.default_rng().uniform(0.05, 0.10, size=(n, 2))
    return X

def discretePDispersion(nPts, points):
    # Find a convex hull in O(N log N)

    # Returned 420 points in testing
    hull = ConvexHull(points)

    # Extract the points forming the hull
    hullpoints = points[hull.vertices,:]

    # Naive way of finding the best pair in O(H^2) time if H is number of points on
    # hull
    hdist = cdist(hullpoints, hullpoints, metric='euclidean')

    # Get the farthest apart points
    bestpair = np.unravel_index(hdist.argmax(), hdist.shape)

    P = np.array([hullpoints[bestpair[0]],hullpoints[bestpair[1]]])

    # Now we have a problem
    print("Finding optimal set...")
    while len(P)<nPts:
        print("P size = {0}".format(len(P)))
        distance_to_P        = cdist(points, P)
        minimum_to_each_of_P = np.min(distance_to_P, axis=1)
        best_new_point_idx   = np.argmax(minimum_to_each_of_P)
        best_new_point = np.expand_dims(points[best_new_point_idx,:],0)
        P = np.append(P,best_new_point,axis=0)

    return P

def imageToCoordinates(image, scaling = 10):
    [height,width] = image.shape
    positiveValues = sum(sum(image>0))
    coordinates = np.zeros([positiveValues,2], dtype=float)

    k = 0
    for i in np.arange(width):
        for j in np.arange(height):
            if image[i,j] > 0:
                coordinates[k,:] = [scaling*j,scaling*i]
                k = k+1

    return np.flip(coordinates,1)
    # return coordinates
    
def importPoints(**kwargs):
    #     nPts = 0, 
    #     dataType='circle', 
    #     dim=2, 
    #     variation=0.1, 
    #     outliers=False, 
    #     # meshFile='../data/61.off', 
    #     filename='',
    #     maxminSubsampling = False,
    #     randomSubsampling = True,
    #     seed=0
    # ):

    data_type_params = {}
    points = None

    if 'path' in kwargs:
        path = kwargs['path']

        if path.endswith('.mat'):
            try:
                points = np.asarray(loadmat(path)['maxmin_points']) # size: 6040 x 24 
            except IOError:
                print('Error: file does not exist.')
        
        elif path.endswith('.off'):
            try: 
                points = readOFF(filename)
            except IOError:
                print('Error: file does not exist.')

        # elif path.endswith('.ply'):
        #     try:
        #         pcl = open3d.io.read_point_cloud(path)
        #     except IOError:
        #         print('Error: file does not exist.')
        #     points = np.asarray(pcl.points) 

        elif path.endswith('.txt'):
            try:
                points = genfromtxt(path, delimeter=',')
            except IOError:
                print('Error: file does not exist.')

            if np.isnan(points).any():
                print('Warning: the dataset contains NaN values, which will be removed.')
                points = points[~np.isnan(points).any(axis=1)]


        filename = os.path.basename(path)
        data_type_params['file'] = filename 

    else:
        if 'data_type' in kwargs:
            dataType = kwargs['data_type']

            match dataType:
                case 'circle':
                    try: 
                        n_pts = kwargs['n_pts']
                    except NameError:
                        print('Please provide the number of points.')

                    if 'variation' in kwargs:
                        variation = kwargs['variation']
                    else:
                        variation = False

                    if 'outliers' in kwargs:
                        outliers = kwargs['outliers']
                    else:
                        outliers = False

                    points = createData(n_pts, 'circle', variation=variation, outliers=outliers)
                    data_type_params = f'{dataType=}' + '_' + f'{n_pts=}'

                case 'sphere':
                    try:
                        ambient_dim = kwargs['ambient_dim']
                    except NameError:
                        print('Please provide the ambient dimension.')

                    try:
                        n_pts = kwargs['n_pts']
                    except NameError:
                        print('Please provide the number of points.')

                    points = createData(n_pts, 'sphere', dim=ambient_dim)

                    # points= shapes.sample_from_sphere(n=nPts, d=(dim-1), seed=0)

                case 'torus':
                    try:
                        n_pts = kwargs['n_pts']
                    except NameError:
                        print('Please provide the ambient dimension.')

                    points = createData(n_pts, 'torus')

                case 'annulus':
                    try:
                        n_pts = kwargs['n_pts']
                    except NameError:
                        print('Please provide the ambient dimension.')

                    if 'inner_radius' in kwargs:
                        inner_radius = kwargs['inner_radius']
                    else:
                        inner_radius = 0.4
                    
                    if 'outer_radius' in kwargs:
                        outer_radius = kwargs['outer_radius']
                    else:
                        outer_radius = 1

                    if 'seed' in kwargs:
                        seed = kwargs['seed']
                    else:
                        seed = 0

                    points= shapes.sample_from_annulus(n=n_pts, r=inner_radius, R=outer_radius, seed=seed)
                
                case 'figure_eight':
                    try:
                        n_pts = kwargs['n_pts']
                    except NameError:
                        print('Please provide the ambient dimension.')

                    if 'figure_eight_scale' in kwargs:
                        figure_eight_scale = kwargs['figure_eight_scale']
                    else:
                        figure_eight_scale = 1

                    if 'figure_eight_neck_size' in kwargs:
                        figure_eight_neck_size = kwargs['figure_eight_neck_size']
                    else:
                        figure_eight_neck_size = 0.1

                    points = figure_eight(n_pts, figure_eight_scale, figure_eight_neck_size)

                case 'mnist':
                    from keras.datasets import mnist

                    (train_X, train_y), (test_X, test_y) = mnist.load_data()
                    coordinates = imageToCoordinates(train_X[1])
                    from matplotlib import pyplot as plt
                    # plt.scatter(coordinates[:,0], coordinates[:,1])
                    # ax = plt.gca()
                    # ax.set_aspect('equal')
                    # plt.show()

                    points = coordinates
                
                case 'turkevs':
                    try:
                        points_all = kwargs['points_all']
                    except NameError:
                        print('Error: please provide the list of all points.')
                    
                    try:
                        labels_all = kwargs['labels_all']
                    except NameError:
                        print('Error: please provide the list of labels.')

                    try:
                        index = kwargs['index']
                    except NameError:
                        print('Error: plese provide the index of the dataset you wish to import.')

                    return points[index], labels_all[index]

                case _:
                    print("Invalid data type specified.")
                    exit()

    return points, data_type_params
        




def get_timestamp():
    return datetime.now().strftime("_%Y%m%d_%H%M%S")


class CustomEncoder(json.JSONEncoder):
    def default(self, obj):

        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
        elif isinstance(obj, topological_computations.Ellipsoid):
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

def saveVarsToFile(dictOfVars,
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

def saveVarsToFile2(dictOfVars,
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

def loadVarsFromFile(filename):
    with open(filename, "r") as f:
        jsonVars = json.load(f)
    
    vars = {}

    if 'dim' in jsonVars:
        vars['dim'] = jsonVars['dim']
    elif 'ambient_dim' in jsonVars:
        vars['ambient_dim'] = jsonVars['ambient_dim']

    if 'rStart' in jsonVars:
        vars['rStart'] = jsonVars['rStart']

    if 'rEnd' in jsonVars:
        vars['rEnd'] = jsonVars['rEnd']

    if 'rStep' in jsonVars:
        vars['rStep'] = jsonVars['rStep']

    if 'rValues' in jsonVars:
        vars['rValues'] = np.asarray(jsonVars['rValues'])

    if 'nbhdSize' in jsonVars:
        vars['nbhdSize'] = jsonVars['nbhdSize']
    elif 'nbhd_size' in jsonVars:
        vars['nbhd_size'] = jsonVars['nbhd_size']

    if 'nPts' in jsonVars:
        vars['nPts'] = jsonVars['nPts']

    if 'points' in jsonVars:
        # vars['points'] = np.asarray(jsonVars['points'])
        vars['points'] = jsonVars['points']

    if 'ellipseList' in jsonVars:
        ellipsoidListRaw = jsonVars['ellipseList']
        ellipsoidList = []
        for ellipsoid in ellipsoidListRaw:
            ellipsoidList.append(topological_computations.Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                     np.asarray(ellipsoid['axesLengths'])))
        vars['ellipsoidList'] = ellipsoidList

    # if 'ellipsoidList' in jsonVars:
    #     ellipsoidListRaw = jsonVars['ellipsoidList']
    if 'ellipsoid_list' in jsonVars:
        ellipsoidListRaw = jsonVars['ellipsoid_list']
        ellipsoidList = []
        for ellipsoid in ellipsoidListRaw:
            ellipsoidList.append(topological_computations.Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                     np.asarray(ellipsoid['axesLengths'])))
        vars['ellipsoid_list'] = ellipsoidList


    if 'simplex_tree_ellipsoids' in jsonVars:
        simplexTreeEllipsoidsRaw = jsonVars['simplex_tree_ellipsoids']
        simplexTreeEllipsoids = gd.SimplexTree()
        for simplexTreeEntry in simplexTreeEllipsoidsRaw:
            simplexTreeEllipsoids.insert(simplexTreeEntry[0],simplexTreeEntry[1])
        vars['simplex_tree_ellipsoids'] = simplexTreeEllipsoids

    if 'simplexTreeRips' in jsonVars:        
        simplexTreeRipsRaw = jsonVars['simplexTreeRips']
        simplexTreeRips = gd.SimplexTree()
        for simplexTreeEntry in simplexTreeRipsRaw:
            simplexTreeRips.insert(simplexTreeEntry[0],simplexTreeEntry[1])
        vars['simplexTreeRips'] = simplexTreeRips

    if 'barcodeEllipsoids' in jsonVars:
        vars['barcode_ellipsoids'] = jsonVars['barcodeEllipsoids']
    elif 'barcode_ellipsoids' in jsonVars:
        vars['barcode_ellipsoids'] = jsonVars['barcode_ellipsoids']

    if 'barcode_rips' in jsonVars:
        vars['barcode_rips'] = jsonVars['barcode_rips']
    elif 'barcodeRips' in jsonVars:
        vars['barcode_rips'] = jsonVars['barcodeRips']

    if 'accs' in jsonVars:
        vars['accs'] = jsonVars['accs']

    if 'accs_file_paths' in jsonVars:
        vars['accs_file_paths'] = jsonVars['accs_file_paths']

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

def readTXT(filename):
    my_data = genfromtxt(filename, delimiter=',')

# from visualisation import visualisationFromFile
from datetime import datetime
from os import listdir
from os.path import isfile, join

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



def readBarcodesInDims0and1(folder, transformation = 'std'):
    filenames = [f for f in listdir(folder) if isfile(join(folder, f)) if f.endswith(".json")]

    pds0 = [None]*1000
    pds1 = [None]*1000

    for filename in filenames:

        if transformation in filename:

            name_to_match = 'Turkevs-' + transformation + '-'

            meshNumber = re.search(rf'{name_to_match}(\d+)', filename).group(1)
            meshNumber = int(meshNumber)
            pd0 = []
            pd1 = []
            print('Reading in the variables... ', end='', flush=True)
            vars = loadVarsFromFile(folder+'/'+filename)
            barcodeEllipsoids = vars['barcodeEllipsoids']
            #barcodeRips = vars['barcodeRips']
            print('Done.')

            for bar in barcodeEllipsoids:
                if bar[0] == 0:
                    pd0.append(bar[1])
                elif bar[0] == 1:
                    pd1.append(bar[1])

            pds0[meshNumber] = pd0
            pds1[meshNumber] = pd1


    # Transform list of 0-dim PDs with different number of cycles into an array of PDs with the same number of cycles.     
    pds0_length = [len(pd) for pd in pds0]
    max_pd_length  = max(pds0_length)    
    pds0 = extend_pds_to_length(pds0, max_pd_length)
    
    # Transform list of 1-dim PDs with different number of cycles into an array of PDs with the same number of cycles.    
    pds1_length = [len(pd) for pd in pds1]
    max_pd_length  = max(pds1_length)   
    pds1 = extend_pds_to_length(pds1, max_pd_length)

    return pds0, pds1
    
def saveMAT(variable, filename, varname='points'):

    pointsDict = {varname: variable}
    savemat(filename, pointsDict)

def savePentagonPoints():
    filenameSaved = 'pentagonsamplesSmall2.mat'

    points = genfromtxt('pentagonsamplesSmall2.txt', delimiter=',')
    if np.isnan(points).any():
        print('The dataset contains NaN values, which will be removed.')
        points = points[~np.isnan(points).any(axis=1)]
    
    saveMAT(points, filenameSaved, varname='pentagonsamples')

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

def get_paths_of_files_in_a_folder(folder, extension='.mat'):
    filenames = [f for f in listdir(folder) if isfile(join(folder, f)) if f.endswith(extension)]
    paths = [os.path.join(folder, f) for f in filenames] 
    return paths
    
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

def extract_barcodes_in_dim0_and_dim1(filename):

    pd0 = []
    pd1 = []
    print('Reading in the variables... ', end='', flush=True)
    vars = loadVarsFromFile(filename)
    print('Done.')
    barcodeEllipsoids = vars['barcodeEllipsoids']

    for bar in barcodeEllipsoids:
        if bar[0] == 0:
            pd0.append(bar[1])
        elif bar[0] == 1:
            pd1.append(bar[1])

    return pd0, pd1

def initialize_pds0_and_pds1(folder):
    paths = get_paths_of_files_in_a_folder(folder)
    pass

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

def remove_dim_from_barcode(barcode):

    result_barcode = []
    for bar in barcode:
            result_barcode.append(bar[1])
    return result_barcode



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

    # print(labels)

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
        = topological_computations.calculate_ellipsoid_barcode(points, nbhd_size, axes_ratios, expansion_dim=expansion_dim)
    barcode_rips, simplex_tree_rips, t_rips \
        = topological_computations.calculate_rips_barcode(points, expansion_dim=expansion_dim)

    # Get the execution time
    t_ellipsoids_over_t_rips = t_ellipsoids / t_rips
    t_total = t_ellipsoids + t_rips
    print('\nThe total execution time is ' + str(t_total) + '\n')

    # Save variables to file
    params_dict = filter_dictionary(vars_to_save, locals())

    if additional_vars_dict: # if not empty
        for key, value in additional_vars_dict.items():
            params_dict[key] = value
    
    saveVarsToFile(params_dict, filename=filename)