from ellipsoids import visualisationFromFile

###### User input ######
filenameLoad = \
    'data/ellipsoids_nPts=100_rStep=0.1_nbhdSize=5_20230926_015938.json'
    #'data/ellipsoids_nPts=50_rStep=0.49999999999999994_nbhdSize=5_20230816_124230.json'
rPlot = 2
########################

visualisationFromFile(filenameLoad, rPlot=rPlot, plotEllipsoids=True)
print('Warning: the simplex tree plot may be inaccurate if the calculations were ' \
        +'not performed for the chosen value of rPlot.')
