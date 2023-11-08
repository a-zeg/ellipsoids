import gudhi as gd

def generateRipsSimplexTree(points):
    print('Creating the Rips complex... ', end='', flush=True)
    ripsComplex = gd.RipsComplex(points=points)
    print('Done.')
    print('Creating the Rips simplex tree... ', end='', flush=True)
    simplexTreeRips = ripsComplex.create_simplex_tree(max_dimension=1)
    print('Done.')
    return simplexTreeRips
