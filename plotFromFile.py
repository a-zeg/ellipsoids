from ellipsoids import visualisationFromFile

###### User input ######
filenameLoad = \
    'data/ellipsoids_nPts=60_rStep=0.49999999999999994_nbhdSize=5_20230815_114301.json'
rPlot = 2
########################

visualisationFromFile(filenameLoad, rPlot=rPlot)
print('Warning: the simplex tree plot may be inaccurate if the calculations were ' \
        +'not performed for the chosen value of rPlot.')
