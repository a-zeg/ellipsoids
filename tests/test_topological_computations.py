import numpy as np
import gudhi as gd
import copy

from scipy.spatial import Delaunay, Voronoi

from ellipsoids.topological_computations import fit_ellipsoid, spherisize_axes
from ellipsoids.topological_computations import generate_ellipsoid_simplex_tree
from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.topological_computations import find_intersection_radius
from ellipsoids.topological_computations import get_max_axes_ratio
from ellipsoids.topological_computations import ellipsoid_intersection
from ellipsoids.topological_computations import reduce_barcode
from ellipsoids.topological_computations import pad_axes_ratios
from ellipsoids.topological_computations import spherisize
from ellipsoids.topological_computations import scale_to_01
from ellipsoids.data_handling import printListOfSimplices
from ellipsoids.topological_computations import adjacent_delaunay_vertices
# from src.topological_computations import get_axes_ratios

def test_fit_ellipsoid():

    center = [0,0]
    neighbourhood = np.array([[-1,0],[0,0],[1,0]])
    axes_ratios = np.array([2,1])

    fitted_ellipsoid = fit_ellipsoid(center, neighbourhood, axes_ratios = axes_ratios)

    assert fitted_ellipsoid.center == center
    assert np.allclose(fitted_ellipsoid.axes[0], np.array([1,0]))
    assert np.allclose(fitted_ellipsoid.axes[1], np.array([0,1]))
    assert np.allclose(fitted_ellipsoid.axes_lengths, [1, 0.5])


def test_ellipsoid_intersection():

    # nearby ellipsoids
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([0.1,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    assert ellipsoid_intersection(ellipsoid1, ellipsoid2, 1) == True

    # far apart ellipsoids
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([5,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    assert ellipsoid_intersection(ellipsoid1, ellipsoid2, 1) == False

    # ellipsoids that touch at one point (along the y-axis)
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    ellipsoid2 = Ellipsoid(np.array([0,1]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    assert ellipsoid_intersection(ellipsoid1, ellipsoid2, 1) == True

    # ellipsoids that touch at one point (along the x-axis)
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    ellipsoid2 = Ellipsoid(np.array([2,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    assert ellipsoid_intersection(ellipsoid1, ellipsoid2, 1) == True

    # ellipsoids that are the same
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))    
    assert ellipsoid_intersection(ellipsoid1, ellipsoid2, 1) == True


def test_find_intersection_radius():

    center1 = np.asarray([0,0])
    axes1 = np.asarray([[1,0],[0,1]])
    axesLengths1 = np.asarray([2,1])

    center2 = np.asarray([1,0])
    axes2 = np.asarray([[1,0],[0,1]])
    axesLengths2 = np.asarray([2,1])


    ellipsoid1 = Ellipsoid(center1, axes1, axesLengths1)
    ellipsoid2 = Ellipsoid(center2, axes2, axesLengths2)

    intersection_radius = find_intersection_radius(ellipsoid1, ellipsoid2)
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


# TODO: put these two into test/utils.py and have a separate test_utils file where they're tested
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

    simplex_tree = generate_ellipsoid_simplex_tree(points, nbhd_size, axes_ratios)
    # printListOfSimplices(simplex_tree_target)
    # printListOfSimplices(simplex_tree[0])

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

    reduced_barcode, _ = reduce_barcode(barcode, nBarsDim0=2, nBarsDim1=0, nBarsDim2=0)
    assert target_barcode == reduced_barcode


    target_barcode = [
        [0, [0,1]]
    ]

    reduced_barcode, _ = reduce_barcode(barcode, nBarsDim0=1, nBarsDim1=0, nBarsDim2=0)
    assert target_barcode == reduced_barcode


def test_max_filtration():

    
    assert True
        

def test_pad_axes_ratios():

    axes_ratios_1 = np.array([3,1])
    dim_1 = 3
    target_axes_ratios_1 = np.array([3,1,1])
    assert np.array_equal(target_axes_ratios_1, pad_axes_ratios(axes_ratios_1, dim_1))


    axes_ratios_2 = np.array([3,1,1,1,1])
    dim_2 = 3
    target_axes_ratios_2 = np.array([3,1,1])
    assert np.array_equal(target_axes_ratios_2, pad_axes_ratios(axes_ratios_2, dim_2))



def test_spherisize():

    axes_lengths = np.array([3,2,1])

    target_axes_lengths_1 = axes_lengths.astype(float)
    s_1 = 0
    assert np.allclose(spherisize(axes_lengths, s_1), target_axes_lengths_1)

    target_axes_lengths_2 = np.array([3,3,3])
    s_2 = 1
    assert np.allclose(spherisize(axes_lengths, s_2), target_axes_lengths_2)

    target_axes_lengths_3 = np.array([3,2.5,2])
    s_3 = 0.5
    assert np.allclose(spherisize(axes_lengths, s_3), target_axes_lengths_3)

    axes_lengths_4 = np.array([1,0.5])
    target_axes_lengths_4 = np.array([1,0.8])
    s_4 = 0.6
    assert np.allclose(spherisize(axes_lengths_4, s_4), target_axes_lengths_4)

    axes_lengths_5 = np.array([5,2,2,1])
    target_axes_lengths_5 = np.array([5, 3.2, 3.2, 2.6])
    s_5 = 0.4
    assert np.allclose(spherisize(axes_lengths_5, s_5), target_axes_lengths_5)



def test_scale_to_01():

    x1 = 2
    min1=1
    max1=3
    target1 = 0.5

    assert np.isclose(scale_to_01(x1,min1,max1), target1)

    x2 = 200
    min2 = 0
    max2 = np.inf
    target2 = 0

    assert np.isclose(scale_to_01(x2,min2,max2), target2)



def test_spherisize_by_filtration():

    axes_lengths_1 = np.array([5,2,2,1])
    spherisize_filtration_1 = 10
    r_1_1 = 12
    target_axes_lengths_1_1 = np.array([5,5,5,5])

    assert np.allclose(
        spherisize_axes(axes_lengths_1, r_1_1, r_spherisize=spherisize_filtration_1),
        target_axes_lengths_1_1
        )

    r_1_2 = -1
    target_axes_lengths_1_2 = axes_lengths_1

    assert np.allclose(
        spherisize_axes(axes_lengths_1, r_1_2, r_spherisize=spherisize_filtration_1),
        target_axes_lengths_1_2
        )

    r_1_3 = 4
    target_axes_lengths_1_3 = np.array([5, 3.2, 3.2, 2.6])

    assert np.allclose(
        spherisize_axes(axes_lengths_1, r_1_3, r_spherisize=spherisize_filtration_1),
        target_axes_lengths_1_3
        )


def test_adjacent_delaunay_vertices():
    points = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]])
    vertex_index = 0
    delaunay = Delaunay(points)
    adjacent_vertices = adjacent_delaunay_vertices(delaunay, vertex_index)
    assert adjacent_vertices == [1,3,4]
