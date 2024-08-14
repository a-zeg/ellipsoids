import numpy as np
import gudhi as gd
import copy

from ellipsoids.topological_computations import fitEllipsoid
from ellipsoids.topological_computations import generateEllipsoidSimplexTree4
from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.topological_computations import findIntersectionRadius
from ellipsoids.topological_computations import get_max_axes_ratio
from ellipsoids.topological_computations import ellipsoidIntersection
from ellipsoids.topological_computations import reduceBarcode
from ellipsoids.topological_computations import padAxesRatios
# from src.topological_computations import get_axes_ratios

def test_fit_ellipsoid():

    center = [0,0]
    neighbourhood = np.array([[-1,0],[0,0],[1,0]])
    axes_ratios = np.array([2,1])

    fitted_ellipsoid = fitEllipsoid(center, neighbourhood, axesRatios = axes_ratios)

    assert fitted_ellipsoid.center == center
    assert np.allclose(fitted_ellipsoid.axes[0], np.array([1,0]))
    assert np.allclose(fitted_ellipsoid.axes[1], np.array([0,1]))
    assert np.allclose(fitted_ellipsoid.axesLengths, [1, 0.5])


def test_ellipsoid_intersection():

    # nearby ellipsoids
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([0.1,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True

    # far apart ellipsoids
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([5,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == False

    # ellipsoids that touch at one point (along the y-axis)
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    ellipsoid2 = Ellipsoid(np.array([0,1]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True

    # ellipsoids that touch at one point (along the x-axis)
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    ellipsoid2 = Ellipsoid(np.array([2,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True

    # ellipsoids that are the same
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))    
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True


def test_find_intersection_radius():

    center1 = np.asarray([0,0])
    axes1 = np.asarray([[1,0],[0,1]])
    axesLengths1 = np.asarray([2,1])

    center2 = np.asarray([1,0])
    axes2 = np.asarray([[1,0],[0,1]])
    axesLengths2 = np.asarray([2,1])


    ellipsoid1 = Ellipsoid(center1, axes1, axesLengths1)
    ellipsoid2 = Ellipsoid(center2, axes2, axesLengths2)

    intersection_radius = findIntersectionRadius(ellipsoid1, ellipsoid2)
    target_intersection_radius = 0.5

    assert np.isclose(intersection_radius, target_intersection_radius, atol=0.01)


def test_get_max_axes_ratio():

    center1 = np.array([0,0])
    axes1 = np.array([[1,0],[0,1]])
    axes_lengths1 = [3,2]
    ellipsoid1 = Ellipsoid(center1, axes1, axes_lengths1)
    target_axes_ratio1 = 3/2

    assert get_max_axes_ratio(ellipsoid1) == target_axes_ratio1 

    center1 = np.array([0,0,0])
    axes1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
    axes_lengths2 = [3,2,1]
    ellipsoid2 = Ellipsoid(center1, axes1, axes_lengths2)
    target_axes_ratio2 = 3

    assert get_max_axes_ratio(ellipsoid2) == target_axes_ratio2 


def _simplex_tree_to_list(simplex_tree: gd.SimplexTree):

    generator = simplex_tree.get_filtration()
    simplex_list = list(generator)

    return simplex_list


def _simplex_trees_equal(st1: gd.SimplexTree, st2: gd.SimplexTree, atol=0.01):
    ''' Compares if two simplex trees are equal.
    i
    The reason a simple list comparison doesn't work is because
    the filtrations might differ slightly due to numerical errors.
    
    This function instead first converts the simplex trees to lists
    of simplices and then checks if all the filtrations of the matching
    elements agree up to a given absolute tolerance atol.'''

    sl1 = _simplex_tree_to_list(st1)
    sl2 = _simplex_tree_to_list(st2)

    sl1_iter = copy.deepcopy(sl1)
    sl2_iter = copy.deepcopy(sl2)

    for splx1 in sl1_iter:
        match_found = False

        for splx2 in sl2_iter:

            if splx1[0] == splx2[0] and np.isclose(splx1[1], splx2[1], atol=atol):
                sl1.remove(splx1)
                sl2.remove(splx2)
                match_found = True
                break
            else: 
                continue

        if match_found == False:
            return False
        
    if (not sl1) and (not sl2):
        return True
    
    return False


def test_simplex_trees_equal():
    
    st1a = gd.SimplexTree()
    st1a.insert([0], filtration=0.0)
    st1b = gd.SimplexTree()
    st1b.insert([0], filtration=0.0)

    assert _simplex_trees_equal(st1a,st1b)

    st2a = gd.SimplexTree()
    st2a.insert([0], filtration=0.0)
    st2a.insert([1], filtration=1.1)
    st2a.insert([0,1], filtration=1.1)
    st2b = gd.SimplexTree()
    st2b.insert([0], filtration=0.0)
    st2b.insert([1], filtration=1.1)
    st2b.insert([0,1], filtration=1.1)

    assert _simplex_trees_equal(st2a, st2b)


def test_generate_ellipsoid_simplex_tree4():

    points = np.array([[-1,0], [0,0], [1,0]])
    nbhd_size = 3
    axes_ratios = np.array([2,1])

    simplex_list = [
        ( [0], 0.0 ),
        ( [1], 0.0 ),
        ( [2], 0.0 ),
        ( [0,1], 0.9985820312499998 ),
        ( [1,2], 0.9985820312499998 ),
        ( [0,2], 1.9995585937499998 )
    ]

    simplex_tree_target = gd.SimplexTree()
    for splx in simplex_list:
        simplex_tree_target.insert(splx[0], filtration=splx[1])

    simplex_tree = generateEllipsoidSimplexTree4(points, nbhd_size, axes_ratios)

    assert _simplex_trees_equal(simplex_tree[0], simplex_tree_target)

    simplex_list_approx = [
        ( [0], 0.0 ),
        ( [1], 0.0 ),
        ( [2], 0.0 ),
        ( [0,1], 1 ),
        ( [1,2], 1 ),
        ( [0,2], 2 )
    ]

    simplex_tree_target_approx = gd.SimplexTree()
    for splx in simplex_list_approx:
        simplex_tree_target_approx.insert(splx[0], filtration=splx[1])

    assert _simplex_trees_equal(simplex_tree[0], simplex_tree_target_approx)


def test_reduce_barcode():

    barcode = [
        [0, [0,1]],
        [0, [-0.5, 1]],
        [1, [-0.3, 0.3]],
        [42, [-10, 353]]
    ]

    target_barcode = [
        [0, [0,1]],
        [0, [-0.5, 1]]
    ]

    reduced_barcode, _ = reduceBarcode(barcode, nBarsDim0=2)
    assert target_barcode == reduced_barcode


    target_barcode = [
        [0, [0,1]]
    ]

    reduced_barcode, _ = reduceBarcode(barcode, nBarsDim0=1)
    assert target_barcode == reduced_barcode


def test_max_filtration():

    
    assert True
        

def pad_axes_ratios():

    axes_ratios1 = np.array([3,1])
    dim1 = 3
    target_axes_ratios1 = np.arary([3,1,1])   
    assert target_axes_ratios1 == padAxesRatios(axes_ratios1, dim1)


    axes_ratios2 = np.array([3,1,1,1,1])
    dim2 = 3
    target_axes_ratios2 = np.arary([3,1,1])   
    assert target_axes_ratios2 == padAxesRatios(axes_ratios2, dim2)