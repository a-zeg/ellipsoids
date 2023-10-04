from ellipsoids import loadVarsFromFile
import gudhi as gd
import numpy as np
import matplotlib.pyplot as plt

def calculatePersistenceIntervalsFromFile(filename):
    vars = loadVarsFromFile(filename)
    dim = len(np.asarray(vars['points'])[0])
    points = vars['points']
    ellipsoidsSimplexTree = vars['simplexTreeEllipsoids']
    ellipsoidsSimplexTree.compute_persistence()
    ripsComplex = gd.RipsComplex(points = points)
    ripsSimplexTree = ripsComplex.create_simplex_tree()
    ripsSimplexTree.expansion(dim)
    ripsSimplexTree.compute_persistence()

    return [ellipsoidsSimplexTree, ripsSimplexTree, dim]

def sumBottleneckDistancesAllDimensions(simplexTree1, simplexTree2, maxDim):
    totalDistance = 0
    for d in np.arange(maxDim+1):
        totalDistance = totalDistance + gd.bottleneck_distance(
            simplexTree1.persistence_intervals_in_dimension(d),
            simplexTree2.persistence_intervals_in_dimension(d)
        )
    return totalDistance

def calculateBottleneckDistances(filename1,filename2):
    ellipsoidsSimplexTree1, ripsSimplexTree1, dim \
        = calculatePersistenceIntervalsFromFile(filename1)
    ellipsoidsSimplexTree2, ripsSimplexTree2, NULL \
        = calculatePersistenceIntervalsFromFile(filename2)
    
    distanceEllipsoids = sumBottleneckDistancesAllDimensions(ellipsoidsSimplexTree1,\
                                                         ellipsoidsSimplexTree2,dim)
    distanceRips = sumBottleneckDistancesAllDimensions(ripsSimplexTree1,\
                                                       ripsSimplexTree2,dim)
    
    return [distanceEllipsoids,distanceRips]

def plotDistances(distanceEllipsoids, distanceRips): 
    width = 0.25
    fig, ax = plt.subplots(layout='constrained')
    ax.bar(0, distanceEllipsoids, width)
    ax.bar(0.3, distanceRips, width)

    ax.set_ylabel('Sum of bottleneck distances in all dimensions')
    ax.set_title('Comparison of bottleneck distances between Ellipsoids and Rips complexes')
    ax.set_xticks([0,0.3])
    ax.set_xticklabels(['Ellipsoid complex', 'Rips complex'])

    plt.show() 

def main():
    filename1 = \
        'data/ellipsoids_nPts=50_rStep=0.1_nbhdSize=5_20230927_145605.json'
    filename2 = \
        'data/ellipsoids_nPts=50_rStep=0.1_nbhdSize=5_20230927_144941.json'
    [distanceEllipsoids, distanceRips] = calculateBottleneckDistances(filename1,filename2)
    plotDistances(distanceEllipsoids, distanceRips)


if __name__ == "__main__":
    main()
