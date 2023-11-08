import numpy as np
import shapes
import figure_eight
from readWriteData import readOFF
from scipy.io import loadmat

def createData(nPoints, type, variation = 0.1, dim=2, outliers=False):
    ''' Generate point cloud data of type 'type' and consiting of 'nPoints' points.
    :nPoints: Number of points to generate
    :type: Type of data (e.g. 'circle', 'ellipse', 'Cassini_oval')
    :variation: Amount of variation from the :type:
    :dim: Not used yet
    :return: nPoints x dim numpy array
    '''
    np.random.seed(0)

    if type == 'circle':
        if outliers is True: nPoints = nPoints - 1
        r = 1
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        x = r*np.cos(t) + variation * np.random.rand(nPoints)
        y = r*np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

        if outliers is True:
            output = np.append(output,[[0,0]],axis=0)

    elif type == 'ellipse':
        t = np.linspace(0, 2*np.pi * (nPoints-1)/nPoints, nPoints)
        a = 2
        b = 1
        x = a * np.cos(t) + variation * np.random.rand(nPoints)
        y = b * np.sin(t) + variation * np.random.rand(nPoints)
        output = np.vstack((x,y)).transpose()

    elif type == 'Cassini_oval':
        t = np.linspace(-1, 1, int(nPoints/2))
        #t = np.sign(t)*np.abs(t)**(1/4)
        x = np.concatenate((t,t)) + variation * np.random.rand(nPoints)
        yh = (t**2 + 0.5) * np.sqrt(1 - t**2)
        y = np.concatenate((-yh, yh)) + variation * np.random.rand(nPoints)
        
        output = np.vstack((x,y)).transpose()

    elif type == 'sphere':
        r = 1

    elif type == 'torus':
        r = 1
        R = 2

        nPtsSampled = 0
        theta = np.zeros([nPoints])
        phi = np.zeros([nPoints])

        # rejection sampling
        while nPtsSampled < nPoints:
            thetaSample = 2 * np.pi * np.random.rand()
            phiSample = 2 * np.pi * np.random.rand()
            W = 2 * np.pi * np.random.rand()

            if W <= (R + r * np.cos(thetaSample))/(R+r):
                theta[nPtsSampled] = thetaSample
                phi[nPtsSampled] = phiSample
                nPtsSampled += 1


            x = (R + r * np.cos(theta)) * np.cos(phi)
            y = (R + r * np.cos(theta)) * np.sin(phi)
            z = r * np.sin(theta)

        output = np.vstack((x,y,z)).transpose()
        

    else:
        raise Exception(type + " is an invalid type of data.")
    
    return output
    
def importPoints(nPts, dataType='circle', dim=2, variation=0.1, outliers=False, meshFile='../data/61.off'):
    match dataType:
        case 'circle':
            return createData(nPts, 'circle', variation=variation, outliers=outliers)
        case 'torus':
            return createData(nPts, 'torus')
        case 'annulus':
            return shapes.sample_from_sphere(n=nPts, d=(dim-1), seed=0)
        case 'figureEight':
            return figure_eight.figure_eight(nPts, 1, 0.1)
        case 'mesh':
            return readOFF('data' + meshFile)
        case 'cyclooctane':
            points = np.asarray(\
                loadmat('pointsCycloOctane.mat')['pointsCycloOctane']) # size: 6040 x 24
            return points[0:nPts,:]

        
    pass