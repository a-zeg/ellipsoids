from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.topological_computations import ellipsoidIntersection
import numpy as np

def test_ellipsoid_intersection_close():
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([0.1,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True

# def test_ellipsoid_intersection_same():
#     ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
#     ellipsoid2 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    
#     assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True

def test_ellipsoid_intersection_far():
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    ellipsoid2 = Ellipsoid(np.array([5,0]), np.array([[1,0],[0,1]]), np.array([1,1]))
    
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == False

def test_ellipsoid_intersection_ellipsoid_touch():
    ellipsoid1 = Ellipsoid(np.array([0,0]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    ellipsoid2 = Ellipsoid(np.array([0,1]), np.array([[1,0],[0,1]]), np.array([1,0.5]))
    
    assert ellipsoidIntersection(ellipsoid1, ellipsoid2, 1) == True