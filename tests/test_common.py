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
    assert ComplexSubtype.RIPS.value == "rips"
    assert ComplexSubtype.ALPHA.value == "alpha"

def test_EllipsoidComplexType_enum_names():
    # Test that the Enum has the correct names
    assert ComplexSubtype.RIPS.name == "RIPS"
    assert ComplexSubtype.ALPHA.name == "ALPHA"

def test_EllipsoidComplexType_enum_members():
    # Test that the Enum contains the correct members
    assert ComplexSubtype("rips") == ComplexSubtype.RIPS
    assert ComplexSubtype("alpha") == ComplexSubtype.ALPHA





# Parameters

# Test to_dict
def test_parameters_to_dict():
    """ Test that to_dict serializes the Parameters object correctly """
    params = Parameters(expansion_dim=3, complex_subtype=ComplexSubtype.ALPHA)

    serialized = params.to_dict()

    # Assertions to check if the dictionary contains the expected values
    assert serialized["expansion_dim"] == 3
    assert serialized["collapse_edges"] is True
    assert serialized["complex_type"] == 'ball'  # Serialized as the enum value (string)
    assert serialized["complex_subtype"] == 'alpha'    # Default value of the enum
    assert serialized["save_simplex_tree"] is False


# Test from_dict
def test_parameters_from_dict():
    """ Test that from_dict deserializes the dictionary correctly """
    data = {
        "expansion_dim": 4,
        "collapse_edges": False,
        "complex_type": 'ball',   # This should map to ComplexType.BALL
        "complex_subtype": 'rips',  # This should map to ComplexSubtype.ALPHA
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
        "complex_subtype": 'alpha',
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
    assert serialized["complex_type"] == 'ellipsoid'  # Inherited from Parameters
    assert serialized["complex_subtype"] == 'rips'  # Inherited from Parameters
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
        "complex_type": 'ellipsoid',  # This should map to ComplexType.ELLIPSOID
        "complex_subtype": 'rips',    # This should map to ComplexSubtype.RIPS
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
        "complex_type": 'ellipsoid',  # This should map to ComplexType.ELLIPSOID
        "complex_subtype": 'rips',    # This should map to ComplexSubtype.RIPS
        # Missing "nbhd_size" here
        "axes_ratios": [4, 3],
        "r_spherisize": 2.5,
        "save_ellipsoid_list": True
    }

    with pytest.raises(KeyError):
        EllipsoidParameters.from_dict(invalid_data)
