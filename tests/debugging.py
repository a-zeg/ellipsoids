from .test_topological_computations import simplex_trees_equal
import gudhi as gd


st2a = gd.SimplexTree()
st2a.insert([0], filtration=0.0)
st2a.insert([1], filtration=1.1)
st2a.insert([0,1], filtration=1.1)
st2b = gd.SimplexTree()
st2b.insert([0], filtration=0.0)
st2b.insert([1], filtration=1.1)
st2b.insert([0,1], filtration=1.1)

simplex_trees_equal(st2a, st2b)