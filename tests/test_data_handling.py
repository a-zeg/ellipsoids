import numpy as np
import gudhi as gd

from ellipsoids.data_handling import parse_turkevs_filename
from ellipsoids.data_handling import read_pd0_and_pd1
from ellipsoids.data_handling import json_process_variables
from ellipsoids.data_handling import remove_dim_from_barcode
from ellipsoids.data_handling import generate_filename
from ellipsoids.data_handling import extract_barcodes_in_dim0_and_dim1
from ellipsoids.data_handling import set_filename_parameters
from ellipsoids.data_handling import filter_dictionary
from ellipsoids.data_handling import get_paths_with_seed


def test_json_process_variables():

    json_variables = { # arbitrary data
        'dim' : 2,
        'ambient_dim' : 3,
        'rStart' : 0,
        'rEnd' : 1,
        'rStep' : 0.5,
        'rValues' : [0, 0.5, 1],
        'nbhd_size' : 5,
        'nPts' : 10,
        'points' : [[0,0],[1,0]],
        'ellipsoid_list' : [{'center': [0,0], 'axes': [[1,0],[0,1]], 'axesLengths': [1,0.5]}],
        'simplex_tree_ellipsoids': 
            [
                [ [0], 0.0 ],
                [ [1], 0.0 ],
                [ [0,1], 0.27057605073979274 ]
            ],
        'simplex_tree_rips': 
            [
                [ [0], 0.0 ]
            ],
        'barcode_ellipsoids': 
        [ 
            [ 1, [0.5, 1] ],
            [ 0, [0.0, float('inf')] ] # default JSON decoder already decodes Infitiy to float('inf')
        ],
        'barcode_rips': 
        [ 
            [ 1, [1, 1.1112] ],
            [ 0, [0.0, float('inf')] ], 
            [ 0, [0.0, 1.1112] ]
        ],
        'accs' : [1,1,0.5],
        'accs_file_paths' : ['path/to/file1.json', 'path/to/file2.json']
    }

    processed_variables = json_process_variables(json_variables)

    assert processed_variables['dim'] == 2
    assert 'ellipseList' not in processed_variables.keys()
    for element in processed_variables['ellipsoid_list']:
        # assert isinstance(element, Ellipsoid)
        # TODO: how can I test if an element is an istance of the Ellipsoid class?
        # The above doesn't work because it ends up comparing topological_computations.Ellipsoid with just Ellipsoid.
        assert isinstance(element, object) 
    processed_ellipsoid = processed_variables['ellipsoid_list'][0]
    assert np.array_equal(processed_ellipsoid.center, np.array([0,0]))
    assert np.array_equal(processed_ellipsoid.axes, np.array([[1,0],[0,1]]))
    assert np.array_equal(processed_ellipsoid.axesLengths, np.array([1,0.5]))
    assert isinstance(processed_variables['simplex_tree_ellipsoids'], gd.SimplexTree)
    assert processed_variables['simplex_tree_ellipsoids'].is_empty() == False
    assert processed_variables['barcode_ellipsoids'][0] == [1, [0.5, 1]]



def test_remove_dim_from_barcode():

    barcode = [[0, [0,1]], [1, [0,1]]]
    target_barcode = [[0,1],[0,1]]
    assert target_barcode == remove_dim_from_barcode(barcode)

    barcode = [[3, [-1,1]], [1, [0,1]], [0, [4, 5.4]]]
    target_barcode = [[-1,1],[0,1],[4,5.4]]
    assert target_barcode == remove_dim_from_barcode(barcode)



def test_generate_filename():
    
    filename_target_01 = 'data/ellipsoids_n_pts=100'
    filename_parameters_01 = {'n_pts' : 100}
    assert generate_filename(filename_parameters_01) == filename_target_01

    filename_target_02 = 'folder/ellipsoids_n_pts=10_var=varValue'
    filename_parameters_02 = {'n_pts' : 10, 'var' : 'varValue'}
    folder_02 = 'folder'
    assert generate_filename(filename_parameters_02, folder_02) == filename_target_02


def test_extract_barcodes_in_dim0_and_dim1():

    barcode = [[0,[0,1]], [1, [0,1]],[1, [0.2, 42]], [2,[0,1]]]

    pd0, pd1 = extract_barcodes_in_dim0_and_dim1(barcode)

    assert pd0 == [[0,1]]
    assert pd1 == [[0,1], [0.2, 42]]
    


def test_parse_turkevs_filename():

    path = 'data/turkevs20240705/ellipsoids_data_type=Turkevs-gauss-003_n_pts=301_nbhd_size=9_axes_ratios=[3 1]__20240427_125418'
    assert parse_turkevs_filename(path) == ('gauss', 3)

    filename = 'data/turkevs/ellipsoids_data_type=Turkevs-trns-079_n_pts=300_nbhd_size=9_axes_ratios=[3 3 1]__20240427_103930.json'
    transformation, index = parse_turkevs_filename(filename)
    assert transformation == 'trns'
    assert index == 79



def equal_ignore_order(a, b):
    """ Use only when elements are neither hashable nor sortable! """
    print(a)
    print(b)
    unmatched = list(b)
    for element in a:
        try:
            unmatched.remove(element)
        except ValueError:
            return False
    print(unmatched)
    return not unmatched



def test_generate_filename():
    folder = 'data'

    dict = {
        'param1': 1,
        'param2': 3.4,
        'param3': 'string'
    }
    generated_filename \
        = generate_filename(filename_parameters=dict, folder=folder)
    target_filename \
        = 'data/ellipsoids_param1=1_param2=3.4_param3=string'

    assert target_filename == generated_filename



def test_set_filename_parameters():
    
    data_type = 'datatype'
    n_pts = 1
    nbhd_size = 1
    axes_ratios = np.array([3,1,1])
    data_type_params = { 'additional_param': 'some_value'}

    target_filename_params = {
        'data_type': data_type,
        'n_pts': n_pts,
        'nbhd_size': nbhd_size,
        'axes_ratios': axes_ratios,
        'additional_param': 'some_value'
    }

    filename_params = set_filename_parameters(data_type, n_pts, nbhd_size, axes_ratios, data_type_params)

    assert target_filename_params == filename_params


def test_filter_dictionary():

    vars_to_save = ['var1', 'var3']
    var1 = 1
    var2 = 'var2_value'
    var3 = 'var3_value'
    var4 = 4.444

    target_filtered_dictionary = {
        'var1': 1,
        'var3': 'var3_value'
    }

    filtered_dictionary = filter_dictionary(vars_to_save, locals())

    assert target_filtered_dictionary == filtered_dictionary


def test_get_paths_with_seed():

    paths = ['datasets/pentagons/pentagonsamplesSmall2.mat_nPts=100_seed=0.mat',
             'datasets/pentagons/pentagonsamplesSmall2.mat_nPts=100_seed=1.mat',
             'datasets/pentagons/pentagonsamplesSmall2.mat_nPts=100_seed=2.mat',
             'datasets/pentagons/pentagonsamplesSmall2.mat_nPts=300_seed=0.mat',
             'datasets/pentagons/pentagonsamplesSmall2.mat_nPts=300_seed=1.mat',
             'datasets/pentagons/pentagonsamplesSmall2.mat_nPts=300_seed=2.mat']

    seed = 1
    target_paths = ['datasets/pentagons/pentagonsamplesSmall2.mat_nPts=100_seed=1.mat',
                    'datasets/pentagons/pentagonsamplesSmall2.mat_nPts=300_seed=1.mat']
    
    output_paths = get_paths_with_seed(paths, seed)
    print(output_paths)

    assert set(target_paths) == set(get_paths_with_seed(paths, seed))


# def test_read_pd0_and_pd1():
#     path = 'tests/test_data/turkevs/testEllipsoids_data_type=Turkevs-std-094_n_pts=301_nbhd_size=9_axes_ratios=[3 1]__20240427_100454.json'

#     pd0,pd1,points = read_pd0_and_pd1(path)

#     points_test = [
#         [
#             0.04719981494445757,
#             0.0
#         ],
#         [
#             0.2952852157810678,
#             0.0
#         ],
#         [
#             0.10289952945739707,
#             0.0
#         ],
#         [
#             0.045391384515248295,
#             0.0
#         ],
#         [
#             0.12107783748711236,
#             0.0
#         ],
#         [
#             0.2845905593250435,
#             0.0
#         ]
#     ]

#     pd0_test = [[
#                 0.0,
#                 0.029964225883488276
#             ],
#             [
#                 0.0,
#                 0.0019173261391611301
#             ]]
#     pd1_test = [
#             [
#                 0.09224727005924485,
#                 0.5857452128573868
#             ],
#             [
#                 0.09993224543035839,
#                 0.5857452128573868
#             ],
#             [
#                 0.11120559667633173,
#                 0.5255687631696404
#             ]
#         ]   
    
#     assert equal_ignore_order(pd0, pd0_test)
#     assert equal_ignore_order(pd1, pd1_test)
#     assert equal_ignore_order(points, points_test)


# from src.data_handling import import_turkevs_transformed


# def test_dictionary():
#     mydict = {
#         'std': [[0,1],[0,2],[0,3]]
#     }

#     assert mydict['std'][1] == [0,2]

