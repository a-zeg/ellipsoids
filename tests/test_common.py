from ellipsoids.common import Ellipsoid
from ellipsoids.common import Dataset
from ellipsoids.common import ComplexType
from ellipsoids.common import ComplexSubtype
from ellipsoids.common import EllipsoidParameters
from ellipsoids.common import Results
from ellipsoids.common import EllipsoidResults
from ellipsoids.common import Parameters


import pytest
import numpy as np
import json

import gudhi as gd

@pytest.fixture
def ellipsoid_data():
    center = np.array([0, 0, 0])
    axes = np.array([1, 0, 0])
    axes_lengths = np.array([2, 3, 4])
    ellipsoid1 = Ellipsoid(center, axes, axes_lengths)
    ellipsoid2 = Ellipsoid(center, axes, axes_lengths)
    ellipsoid3 = Ellipsoid(np.array([1, 1, 1]), axes, axes_lengths)
    return ellipsoid1, ellipsoid2, ellipsoid3, center, axes, axes_lengths

def test_Ellipsoid_constructor(ellipsoid_data):
    ellipsoid1, _, _, center, axes, axes_lengths = ellipsoid_data

    assert np.array_equal(ellipsoid1.center, center)
    assert np.array_equal(ellipsoid1.axes, axes)
    assert np.array_equal(ellipsoid1.axes_lengths, axes_lengths)

def test_Ellipsoid_equality(ellipsoid_data):
    ellipsoid1, ellipsoid2, ellipsoid3, _, _, _ = ellipsoid_data

    assert ellipsoid1 == ellipsoid2
    assert ellipsoid1 != ellipsoid3

def test_Ellipsoid_to_dict(ellipsoid_data):
    ellipsoid1, _, _, center, axes, axes_lengths = ellipsoid_data
    expected_dict = {
        "center": center,
        "axes": axes,
        "axes_lengths": axes_lengths
    }

    assert ellipsoid1.to_dict() == expected_dict




@pytest.fixture
def dataset_data():
    points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    data_type = "example_type"
    dataset1 = Dataset(points, data_type)
    return dataset1, points, data_type

def test_Dataset_constructor(dataset_data):
    dataset1, points, data_type = dataset_data

    assert np.array_equal(dataset1.points, points)
    assert dataset1.data_type == data_type

def test_Dataset_n_points(dataset_data):
    dataset1, _, _ = dataset_data

    assert dataset1.n_points() == len(dataset1.points)

def test_Dataset_ambient_dim(dataset_data):
    dataset1, points, _ = dataset_data

    assert dataset1.ambient_dim() == len(points[0])

def test_Dataset_to_dict(dataset_data):
    dataset1, points, data_type = dataset_data

    expected_dict = {
        "points": points.tolist(),
        "data_type": data_type
    }
    assert dataset1.to_dict() == expected_dict


# Test to_dict for Dataset
def test_dataset_to_dict():
    """ Test that to_dict serializes the Dataset object correctly """
    points = np.array([[1, 2], [3, 4], [5, 6]])
    dataset = Dataset(points=points, data_type="example_type")

    serialized = dataset.to_dict()

    # Assertions to check if the dictionary contains the expected values
    assert serialized["data_type"] == "example_type"
    assert serialized["points"] == [[1, 2], [3, 4], [5, 6]]  # np.ndarray converted to list


# Test from_dict for Dataset
def test_dataset_from_dict():
    """ Test that from_dict deserializes the dictionary correctly into Dataset """
    data = {
        "data_type": "example_type",
        "points": [[1, 2], [3, 4], [5, 6]]  # List to be converted to np.ndarray
    }

    dataset = Dataset.from_dict(data)

    # Assertions to check if the deserialized Dataset object has the expected values
    assert dataset.data_type == "example_type"
    assert np.array_equal(dataset.points, np.array([[1, 2], [3, 4], [5, 6]]))  # Ensure the points are the same


# Test invalid dictionary for from_dict (e.g., missing required fields)
def test_dataset_from_dict_invalid():
    """ Test that from_dict raises an error when required fields are missing """
    invalid_data = {
        "data_type": "example_type",
        # Missing "points" here
    }

    with pytest.raises(KeyError):
        Dataset.from_dict(invalid_data)







def test_EllipsoidComplexType_enum_values():
    assert ComplexSubtype.RIPS.value == "RIPS"
    assert ComplexSubtype.ALPHA.value == "ALPHA"

def test_EllipsoidComplexType_enum_names():
    # Test that the Enum has the correct names
    assert ComplexSubtype.RIPS.name == "RIPS"
    assert ComplexSubtype.ALPHA.name == "ALPHA"

def test_EllipsoidComplexType_enum_members():
    # Test that the Enum contains the correct members
    assert ComplexSubtype("RIPS") == ComplexSubtype.RIPS
    assert ComplexSubtype("ALPHA") == ComplexSubtype.ALPHA





# Parameters

# Test to_dict
def test_parameters_to_dict():
    """ Test that to_dict serializes the Parameters object correctly """
    params = Parameters(expansion_dim=3, complex_subtype=ComplexSubtype.ALPHA)

    serialized = params.to_dict()

    # Assertions to check if the dictionary contains the expected values
    assert serialized["expansion_dim"] == 3
    assert serialized["collapse_edges"] is True
    assert serialized["complex_type"] == 'BALL'  # Serialized as the enum value (string)
    assert serialized["complex_subtype"] == 'ALPHA'    # Default value of the enum
    assert serialized["save_simplex_tree"] is False


# Test from_dict
def test_parameters_from_dict():
    """ Test that from_dict deserializes the dictionary correctly """
    data = {
        "expansion_dim": 4,
        "collapse_edges": False,
        "complex_type": 'BALL',   # This should map to ComplexType.BALL
        "complex_subtype": 'RIPS',  # This should map to ComplexSubtype.ALPHA
        "save_simplex_tree": True
    }

    params = Parameters.from_dict(data)

    # Assertions to check if the deserialized Parameters object has the expected values
    assert params.expansion_dim == 4
    assert params.collapse_edges is False
    assert params.complex_type == ComplexType.BALL  # Deserialized back to the Enum member
    assert params.complex_subtype == ComplexSubtype.RIPS  # Deserialized back to the Enum member
    assert params.save_simplex_tree is True


# Test invalid dictionary for from_dict (e.g., missing required fields)
def test_parameters_from_dict_invalid():
    """ Test that from_dict raises an error when required fields are missing """
    invalid_data = {
        "expansion_dim": 4,
        "collapse_edges": False,
        # Missing "complex_type"
        "complex_subtype": 'ALPHA',
        "save_simplex_tree": True
    }

    with pytest.raises(KeyError):
        Parameters.from_dict(invalid_data)






# EllipsoidParameters

# Test to_dict for EllipsoidParameters
def test_ellipsoid_parameters_to_dict():
    """ Test that to_dict serializes the EllipsoidParameters object correctly """
    params = EllipsoidParameters(
        nbhd_size=5,
        axes_ratios=np.array([3, 1]),
        r_spherisize=1.5,
        save_ellipsoid_list=False
    )

    serialized = params.to_dict()

    # Assertions to check if the dictionary contains the expected values
    assert serialized["expansion_dim"] == 2  # Inherited from Parameters
    assert serialized["collapse_edges"] is True  # Inherited from Parameters
    assert serialized["complex_type"] == 'ELLIPSOID'  # Inherited from Parameters
    assert serialized["complex_subtype"] == 'RIPS'  # Inherited from Parameters
    assert serialized["save_simplex_tree"] is False  # Inherited from Parameters
    assert serialized["nbhd_size"] == 5
    assert all(x == y for x, y in zip(serialized["axes_ratios"], [3, 1]))  # Convert np.ndarray to list
    assert serialized["r_spherisize"] == 1.5
    assert serialized["save_ellipsoid_list"] is False


# Test from_dict for EllipsoidParameters
def test_ellipsoid_parameters_from_dict():
    """ Test that from_dict deserializes the dictionary correctly into EllipsoidParameters """
    data = {
        "expansion_dim": 2,
        "collapse_edges": True,
        "complex_type": 'ELLIPSOID',  # This should map to ComplexType.ELLIPSOID
        "complex_subtype": 'RIPS',    # This should map to ComplexSubtype.RIPS
        "save_simplex_tree": False,
        "nbhd_size": 6,
        "axes_ratios": [4, 3],  # This should map to np.ndarray([4, 3])
        "r_spherisize": 2.5,
        "save_ellipsoid_list": True
    }

    params = EllipsoidParameters.from_dict(data)

    # Assertions to check if the deserialized Parameters object has the expected values
    assert params.expansion_dim == 2
    assert params.collapse_edges is True
    assert params.complex_type == ComplexType.ELLIPSOID
    assert params.complex_subtype == ComplexSubtype.RIPS
    assert params.save_simplex_tree is False
    assert params.nbhd_size == 6
    assert np.array_equal(params.axes_ratios, np.array([4, 3]))
    assert params.r_spherisize == 2.5
    assert params.save_ellipsoid_list is True


# Test invalid dictionary for from_dict (e.g., missing required fields)
def test_ellipsoid_parameters_from_dict_invalid():
    """ Test that from_dict raises an error when required fields are missing """
    invalid_data = {
        "expansion_dim": 2,
        "collapse_edges": True,
        "complex_type": 'ELLIPSOID',  # This should map to ComplexType.ELLIPSOID
        "complex_subtype": 'RIPS',    # This should map to ComplexSubtype.RIPS
        # Missing "nbhd_size" here
        "axes_ratios": [4, 3],
        "r_spherisize": 2.5,
        "save_ellipsoid_list": True
    }

    with pytest.raises(KeyError):
        EllipsoidParameters.from_dict(invalid_data)










# import pytest
# from unittest.mock import patch, MagicMock
# from datetime import datetime
# from ellipsoids.common import Experiment, Dataset, Parameters, Results, ComplexType, ComplexSubtype
# from ellipsoids.data_handling import CustomEncoder

# @pytest.fixture
# def mock_objects():

#     mock_dataset = MagicMock(Dataset, autospec = True)
#     mock_parameters = MagicMock(Parameters)
#     mock_results = MagicMock(Results)

#     mock_dataset.data_type = 'example_data'
#     mock_dataset.points = np.asarray([[0,0],[0,1]])
#     # mock_dataset.n_points.return_value = len(mock_dataset.points)  # Mock n_points() to return 2
#     # mock_dataset.n_points = 2
#     mock_parameters.complex_type = ComplexType.BALL
#     mock_parameters.complex_subtype = ComplexSubtype.RIPS
#     mock_results.execution_time = 123.45

#     mock_dataset.to_dict.return_value = {'data_type': 'example_data', 'points': [[0,0],[0,1]]}
#     mock_parameters.to_dict.return_value = {'complex_type': 'BALL', 'complex_subtype': 'Rips'}
#     mock_results.to_dict.return_value = {'execution_time': 123.45}

#     return mock_dataset, mock_parameters, mock_results


# @pytest.fixture
# def mock_current_time():
#     with patch('ellipsoids.common.datetime') as mock_datetime:
#         mock_datetime.now.return_value = datetime(2025, 2, 25, 12, 0, 0)
#         yield mock_datetime


# def test_save_to_json(mock_objects, mock_current_time):
#     mock_dataset, mock_parameters, mock_results = mock_objects

#     experiment = Experiment(mock_dataset, mock_parameters)
#     experiment.results = mock_results

#     # Patch open and json.dumps to mock file I/O
#     with patch("builtins.open", new_callable=MagicMock) as mock_open, patch("json.dumps") as mock_json_dumps:
#         expected_filename = 'example_data-2_BALL-Rips__20250225_120000.json'

#         # Mock the return value of json.dumps to prevent actual file writing
#         mock_json_dumps.return_value \
#             = '{"dataset": {"data_type": "example_data", "points": [[0,0],[0,1]]}, \
#                 "parameters": {"complex_type": "BALL", "complex_subtype": "Rips"}, \
#                 "results": {"execution_time": 123.45}}'

#         # Print debug output for the filename generation
#         print(f"Generated filename: {expected_filename}")

#         experiment.save_to_json()

#         print(f"open call arguments: {mock_open.call_args}")
#         mock_open.assert_called_once_with(expected_filename, 'w')

#         expected_data = {
#             'dataset': {'data_type': 'example_data', 'points': [[0,0],[0,1]]},
#             'parameters': {'complex_type': 'BALL', 'complex_subtype': 'Rips'},
#             'results': {'execution_time': 123.45},
#         }
#         print(f"json.dumps call arguments: {mock_json_dumps.call_args}")
#         mock_json_dumps.assert_called_once_with(expected_data, cls=CustomEncoder, indent=4)



import pytest
from datetime import datetime
from unittest.mock import patch, MagicMock
import numpy as np
from ellipsoids.common import Experiment, Dataset, Parameters, Results, ComplexType, ComplexSubtype
from ellipsoids.data_handling import CustomEncoder
import json

@pytest.fixture
def mock_current_time():
    with patch('ellipsoids.common.datetime') as mock_datetime:
        mock_datetime.now.return_value = datetime(2025, 2, 25, 12, 0, 0)
        yield mock_datetime

# Integration test function for save_to_json
def test_save_to_json_integration(mock_current_time):
    dataset = Dataset(data_type="example_data", points=np.array([[0, 0], [0, 1]]))
    parameters = Parameters(complex_type=ComplexType.BALL, complex_subtype=ComplexSubtype.RIPS)
    results = Results()
    results.execution_time = 123.45
    experiment = Experiment(dataset, parameters)
    experiment.results = results

    with patch("builtins.open", new_callable=MagicMock) as mock_open, patch("json.dumps") as mock_json_dumps:
        # Define the expected filename and mock return value for json.dumps
        expected_filepath = 'data/example_data-2_BALL-Rips__20250225_120000.json'
        # Call the save_to_json method (this is where we want to test the real behavior)
        experiment.save_to_json()

        # Assert that the correct file was created
        mock_open.assert_called_once_with(expected_filepath, 'w')

        # Assert that json.dumps was called with the correct data
        expected_data = {
            'dataset': {'data_type': 'example_data', 'points': [[0, 0], [0, 1]]},
            "parameters": {"expansion_dim": 2, "collapse_edges": True,
                            "complex_type": "BALL", "complex_subtype": "RIPS",
                            "save_simplex_tree": False},
            'results': {"barcode": [], "simplex_tree": gd.SimplexTree(), 'execution_time': 123.45}, # gd.SimplexTree() should just create an empty object
        }
        mock_json_dumps.assert_called_once_with(expected_data, cls=CustomEncoder, indent=4)
