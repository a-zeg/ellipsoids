from ellipsoids.data_handling import json_process_variables
from ellipsoids.data_handling import remove_dim_from_barcode
from ellipsoids.data_handling import generate_filename
from ellipsoids.data_handling import extract_barcodes_in_dim0_and_dim1
from ellipsoids.data_handling import parse_turkevs_filename
from ellipsoids.topological_computations import Ellipsoid
import gudhi as gd
import numpy as np
import inspect

def arbitrary_json_variable():

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

    return json_variables

def test_json_process_variables():

    json_variables = arbitrary_json_variable()
    processed_variables = json_process_variables(json_variables)

    assert processed_variables['dim'] == 2
    assert 'ellipseList' not in processed_variables.keys()
    for element in processed_variables['ellipsoid_list']:
        # assert isinstance(element, Ellipsoid)
        # TODO_Bastian: how can I test if an element is an istance of the Ellipsoid class?
        # The above doesn't work because it ends up comparing topological_computations.Ellipsoid with just Ellipsoid.
        assert isinstance(element, object) 
    assert isinstance(processed_variables['simplex_tree_ellipsoids'], gd.SimplexTree)
    assert processed_variables['simplex_tree_ellipsoids'].is_empty() == False
    assert processed_variables['barcode_ellipsoids'][0] == [1, [0.5, 1]]



def test_remove_dim_from_barcode():

    barcode = [[0, [0,1]], [1, [0,1]]]
    assert remove_dim_from_barcode(barcode) == [[0,1],[0,1]]



def test_generate_filename():
    
    filename_01 = 'data/ellipsoids_n_pts=100'
    filename_parameters_01 = {'n_pts' : 100}
    assert generate_filename(filename_parameters_01) == filename_01

    filename_02 = 'folder/ellipsoids_n_pts=10_var=varValue'
    filename_parameters_02 = {'n_pts' : 10, 'var' : 'varValue'}
    folder_02 = 'folder'
    assert generate_filename(filename_parameters_02, folder_02) == filename_02


def test_extract_barcodes_in_dim0_and_dim1():

    barcode = [[0,[0,1]], [1, [0,1]],[1, [0.2, 42]], [2,[0,1]]]

    pd0, pd1 = extract_barcodes_in_dim0_and_dim1(barcode)

    assert pd0 == [[0,1]]
    assert pd1 == [[0,1], [0.2, 42]]
    


def test_parse_turkevs_filename():

    path = 'data/turkevs20240705/ellipsoids_data_type=Turkevs-gauss-003_n_pts=301_nbhd_size=9_axes_ratios=[3 1]__20240427_125418'
    assert parse_turkevs_filename(path) == ('gauss', 3)



    