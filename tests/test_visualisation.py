import numpy as np
import gudhi as gd
import pytest
import collections
from unittest.mock import MagicMock

from ellipsoids.common import Experiment
from ellipsoids.common import Results
from ellipsoids.common import PlotParameters
from ellipsoids.visualisation import quick_dim_to_bars
from ellipsoids.visualisation import reduce_barcode_descending
from ellipsoids.visualisation import find_max_end
from ellipsoids.visualisation import axis_end_experiments
# from ellipsoids.visualisation import BarsInDim



def test_quick_dim_to_bars():

    assert quick_dim_to_bars([0,1,2]) == {0:0, 1:1, 2:2}
    assert quick_dim_to_bars([0,0,1]) == {0:0, 1:0, 2:1}
    assert quick_dim_to_bars([10]) == {0:10}
    assert quick_dim_to_bars([10]) != {0:9}



def test_reduce_barcode_to_longest():

    barcode = [
        (0, (0,1)),
        (0, (-0.5, 1)),
        (1, (-0.3, 0.3)),
        (42, (-10, 353))
    ]

    target_barcode_1 = [
        (0, (0,1)),
        (0, (-0.5, 1))
    ]
    reduced_barcode_1 = reduce_barcode_descending(barcode, {0:5})
    assert collections.Counter(target_barcode_1) == collections.Counter(reduced_barcode_1) # compares lists if order doesn't matte

    target_barcode_2 = [
        (0, (-0.5,1))
    ]
    reduced_barcode_2 = reduce_barcode_descending(barcode, {0:1, 1:0})
    assert collections.Counter(target_barcode_2) == collections.Counter(reduced_barcode_2)

    target_barcode_3 = [
        (42, (-10, 353))
    ]
    reduced_barcode_3 = reduce_barcode_descending(barcode, {0:0, 42:2})
    assert collections.Counter(target_barcode_3) == collections.Counter(reduced_barcode_3)

    non_target_barcode_4 = [
        (41, (-10, 353))
    ]
    reduced_barcode_4 = reduce_barcode_descending(barcode, {0:0, 42:2})
    assert collections.Counter(non_target_barcode_4) != collections.Counter(reduced_barcode_4)


def test_find_max_end():

    barcode = [
        (0, (0,1)),
        (0, (-0.5, 1000)),
        (1, (-0.3, 0.3)),
        (42, (-10, 353))
    ]
    target_max_end = 1000
    assert target_max_end == find_max_end(barcode)

    barcode = [
        (0, (0,1)),
        (0, (-1000, -900)),
        (1, (-0.3, 0.3)),
        (42, (-10, 353))
    ]
    target_max_end = 353
    assert target_max_end == find_max_end(barcode)

    barcode = [
        (0, (0,1)),
        (0, (-1000, -900)),
        (1, (-0.3, 0.3)),
        (42, (-10, 353))
    ]
    non_target_max_end = 0
    assert non_target_max_end != find_max_end(barcode)



# Assuming axis_end_experiments, Experiment, Results, and other necessary classes are imported

@pytest.fixture
def mock_experiments():
    # Create mocked experiment instances
    experiment1 = MagicMock(spec=Experiment)
    experiment2 = MagicMock(spec=Experiment)

    # because results is not in the spec (it is a dynamic attribute), we need:
    experiment1.results = MagicMock(spec=Results)
    experiment2.results = MagicMock(spec=Results)
    experiment1.plot_parameters = MagicMock(spec=PlotParameters)
    experiment2.plot_parameters = MagicMock(spec=PlotParameters)

    # Mock the results and plot_parameters for experiment1
    experiment1.results.barcode = [
        (0, [0.0, 2.0]),  # Dimension 0, barcode [(0.0, 2.0)]
        (1, [1.0, 3.0]),  # Dimension 1, barcode [(1.0, 3.0)]
    ]
    experiment1.plot_parameters.n_bars = {0: 1, 1: 1}  # Use 1 bar per dimension

    # Mock the results and plot_parameters for experiment2
    experiment2.results.barcode = [
        (0, [0.5, 1.5]),  # Dimension 0, barcode [(0.5, 1.5)]
        (1, [0.0, 4.0]),  # Dimension 1, barcode [(0.0, 4.0)]
    ]
    experiment2.plot_parameters.n_bars = {0: 1, 1: 1}  # Use 1 bar per dimension

    return [experiment1, experiment2]


def test_axis_end_experiments(mock_experiments):
    # Call axis_end_experiments with the mocked experiments
    axis_end_value = axis_end_experiments(mock_experiments)

    # Assert that the result is what we expect
    # For this case, the max end value should be 4.0 (from experiment2, dimension 1)
    assert pytest.approx(axis_end_value, 0.001) == 1.1 * 4.0
