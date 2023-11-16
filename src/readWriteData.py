import json
import gudhi as gd
from datetime import datetime
from ellipsoidSimplexTree import Ellipsoid
import numpy as np
import os


class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Ellipsoid):
            obj_data = {
                "center": obj.center,
                "axes": obj.axes.tolist(),
                "axesLengths": obj.axesLengths.tolist()
            }
            return obj_data
        elif isinstance(obj,gd.SimplexTree):
            return list(obj.get_filtration())
        return json.JSONEncoder.default(self, obj)

def saveVarsToFile(dictOfVars,
                   filename=datetime.now().strftime("data/test.json")):
    print('Saving data to file...')
    json_string = json.dumps(dictOfVars, cls=CustomEncoder, indent=4)
    with open(filename, 'w') as outfile:
        outfile.write(json_string)
    print("Data saved to file " + filename + '.')

def saveVarsToFile2(dictOfVars,
                   filename=datetime.now().strftime("data/test.json"), 
                   addToStart=False):
    print('Saving data to file...')
    json_string = json.dumps(dictOfVars, cls=CustomEncoder, indent=4)
    if os.path.exists(filename): # if the file already exists, append to it
        if addToStart:
            index = 0
            with open(filename, "r") as f:
                contents = f.readlines()
            contents.insert(index, dictOfVars)
            with open(filename, "w") as f:
                contents = "".join(str(c) for c in contents)
                f.write(contents)
        else:
            with open(filename, "a") as f:
                f.write(json_string)
    else:
        with open(filename, 'w') as outfile:
            outfile.write(json_string)
    print("Data saved to file " + filename + '.')

def loadVarsFromFile(filename):
    with open(filename, "r") as f:
        jsonVars = json.load(f)
    
    vars = {}

    if 'dim' in jsonVars:
        vars['dim'] = jsonVars['dim']

    if 'rStart' in jsonVars:
        vars['rStart'] = jsonVars['rStart']

    if 'rEnd' in jsonVars:
        vars['rEnd'] = jsonVars['rEnd']

    if 'rStep' in jsonVars:
        vars['rStep'] = jsonVars['rStep']

    if 'rValues' in jsonVars:
        vars['rValues'] = np.asarray(jsonVars['rValues'])

    if 'nbhdSize' in jsonVars:
        vars['nbhdSize'] = jsonVars['nbhdSize']

    if 'nPts' in jsonVars:
        vars['nPts'] = jsonVars['nPts']

    if 'points' in jsonVars:
        vars['points'] = np.asarray(jsonVars['points'])

    if 'ellipseList' in jsonVars:
        ellipsoidListRaw = jsonVars['ellipseList']
        ellipsoidList = []
        for ellipsoid in ellipsoidListRaw:
            ellipsoidList.append(Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                     np.asarray(ellipsoid['axesLengths'])))
        vars['ellipsoidList'] = ellipsoidList

    if 'ellipsoidList' in jsonVars:
        ellipsoidListRaw = jsonVars['ellipsoidList']
        ellipsoidList = []
        for ellipsoid in ellipsoidListRaw:
            ellipsoidList.append(Ellipsoid(ellipsoid['center'], np.asarray(ellipsoid['axes']), \
                                     np.asarray(ellipsoid['axesLengths'])))
        vars['ellipsoidList'] = ellipsoidList

    if 'simplexTreeEllipsoids' in jsonVars:
        simplexTreeEllipsoidsRaw = jsonVars['simplexTreeEllipsoids']
        simplexTreeEllipsoids = gd.SimplexTree()
        for simplexTreeEntry in simplexTreeEllipsoidsRaw:
            simplexTreeEllipsoids.insert(simplexTreeEntry[0],simplexTreeEntry[1])
        vars['simplexTreeEllipsoids'] = simplexTreeEllipsoids

    if 'simplexTreeRips' in jsonVars:        
        simplexTreeRipsRaw = jsonVars['simplexTreeRips']
        simplexTreeRips = gd.SimplexTree()
        for simplexTreeEntry in simplexTreeRipsRaw:
            simplexTreeRips.insert(simplexTreeEntry[0],simplexTreeEntry[1])
        vars['simplexTreeRips'] = simplexTreeRips

    if 'barcodeEllipsoids' in jsonVars:
        vars['barcodeEllipsoids'] = jsonVars['barcodeEllipsoids']

    if 'barcodeRips' in jsonVars:
        vars['barcodeRips'] = jsonVars['barcodeRips']

    return vars
     
def readOFF(filename):
    file = open(filename, 'r')
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    nVerts, nFaces, nEdges = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')] for i_vert in range(nVerts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for i_face in range(nFaces)]
    return np.asarray(verts)

def printListOfSimplices(simplexTree):
    generator = simplexTree.get_filtration()
    simplexList = list(generator)
    for splx in simplexList:
        print(splx)
