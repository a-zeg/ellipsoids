from data_handling import json_process_variables
from topological_computations import Ellipsoid

def test_json_process_variables():
    json_variables = {
        'dim' : 2,
        'ambient_dim' : 3,
        'rStart' : 0,
        'rEnd' : 1,
        'rStep' : 0.5,
        'rValues' : [0, 0.5, 1],
        'nbhd_size' : 5, # could also be nbhdSize in older versions
        'nPts' : 10,
        'points' : [[0,0],[1,0]],
        'ellipsoid_list' : [{'center': [0,0], 'axes': [[1,0],[0,1]], 'axesLengths': [1,0.5]}],
        'simplex_tree_ellipsoids': None,
        'simplex_tree_rips': None,
        'barcode_ellipsoids': None,
        'barcode_rips': None,
        'accs' : [1,1,0.5],
        'accs_file_paths' : ['path/to/file.json']
    }

    processed_variables = json_process_variables(json_variables)

    assert processed_variables['dim'] == 2
    assert processed_variables['ellipse_list'] is None
    assert type(processed_variables['ellipsoid_list']) is Ellipsoid
    assert type(processed_variables['ellipsoid_list']) is int
