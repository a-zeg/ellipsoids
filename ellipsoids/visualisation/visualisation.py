import gudhi as gd
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime

from ellipsoids.topological_computations import Ellipsoid
from ellipsoids.data_handling import read_variables
from ellipsoids.topological_computations import reduceBarcode
from ellipsoids.visualisation.barcodePlotting import plot_persistence_barcode, plot_persistence_density


def plotEllipse(ellipse: Ellipsoid, color='grey', r=1, axes=None):
    sampleRate = 100
    t = np.linspace(0, 2*np.pi, sampleRate)
    xTemp = r*ellipse.axesLengths[0]*np.cos(t)
    yTemp = r*ellipse.axesLengths[1]*np.sin(t)
    x = ellipse.center[0] + ellipse.axes[0,0]*xTemp + ellipse.axes[1,0]*yTemp
    y = ellipse.center[1] + ellipse.axes[0,1]*xTemp + ellipse.axes[1,1]*yTemp
    if axes is None:
        plt.plot(x,y,c=color)
    else:
        axes.plot(x,y,c=color)

def plotEllipsoid(ellipsoid: Ellipsoid, color='grey', r=1, axes=None):
    # see https://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
    sampleRate = 100

    rx = r * ellipsoid.axesLengths[0]
    ry = r * ellipsoid.axesLengths[1]
    rz = r * ellipsoid.axesLengths[2]
    
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, sampleRate)
    v = np.linspace(0, np.pi, sampleRate)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    all = np.concatenate((np.reshape(x, [-1,1]), np.reshape(y, [-1,1]), np.reshape(z, [-1,1])), axis=1)
    allTransformed = ellipsoid.axes @ np.transpose(all)
    
    x = np.reshape(allTransformed[0,:],(100,100)) + ellipsoid.center[0]
    y = np.reshape(allTransformed[1,:],(100,100)) + ellipsoid.center[1]
    z = np.reshape(allTransformed[2,:],(100,100)) + ellipsoid.center[2]

    if axes is None:
        plt.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)
    else:
        axes.plot_surface(x,y,z, rstride=4, cstride=4, color=color, alpha = 0.2)

def plotEllipses(ellipseList, r, axes=None):
    for ellipse in ellipseList:
        plotEllipse(ellipse, r=r, axes=axes)

def plotEllipsoids(ellipsoidList, r, axes=None):
    for ellipsoid in ellipsoidList:
        plotEllipsoid(ellipsoid, r = r, axes=axes)

def plotCircle(point, r=1, color='grey', axes=None):
    sample_rate = 100
    t = np.linspace(0, 2*np.pi, sample_rate)
    x = point[0] + r*np.cos(t)
    y = point[1] + r*np.sin(t)
    if axes is None:
        plt.plot(x,y,c=color)
    else:
        axes.plot(x,y,c=color)

def plotCircles(points, r=1, axes=None):
    for point in points:
        plotCircle(point, r=r, axes=axes)

def plotSimplexTree(points, simplexTree, r, axes):
    dim = len(points[0])
    if dim > 3:
        raise Exception('Error: Attempting to plot simplex tree in dimension higher than 3.')
    generator = simplexTree.get_filtration()
    simplexList = list(generator)

    if axes is None:
        for splx in simplexList:
            if splx[1] <= r:
                vertices = splx[0]
                match len(vertices):
                    case 1:
                        plt.scatter(*np.transpose(points[vertices]), c='k', zorder=100)
                        # points[vertices] gives us points forming the vertices of the simplex
                        # transposing them and taking the * operator returns x-, y-. and z-coords separately
                    case 2:
                        plt.plot(*np.transpose(points[vertices]), c='r')
                    case 3:
                        if dim == 2:
                            plt.fill(*np.transpose(points[vertices]), c='r', alpha=0.1)
    else:
        idx = 1
        for splx in simplexList:
            idx = idx + 1
            
            # if splx[1] <= r:
            if splx[1] <= 2*r: # alt20230927_2: 2r so that it's comparable to Rips
                vertices = splx[0]
                match len(vertices):
                    case 1:
                        axes.scatter(*np.transpose(points[vertices]), c='k', zorder=100)
                    case 2:
                        axes.plot(*np.transpose(points[vertices]), c='r')
                    case 3:
                        if dim == 2:
                            axes.fill(*np.transpose(points[vertices]), c='r', alpha=0.1)

def plotDataPoints(points, axes=None):
    if axes is None:
        plt.scatter(points[:,0],points[:,1])
    else:
        axes.scatter(points[:,0],points[:,1])

def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a

def visualisation(**kwargs):
    print('Generating plots...')

    xAxisEnd = kwargs['xAxisEnd']
    barcodeEllipsoids = kwargs['barcodeEllipsoids']
    barcodeRips = kwargs['barcodeRips']
    plotDensity = False
    if 'plotDensity' in kwargs:
        plotDensity = kwargs['plotDensity']
    

    figsize = (10,10)

    drawPoints = False
    drawEllipsoids = False
    drawEllipsoidsSimplexTree = False

    if 'drawPoints' in kwargs:
        drawPoints = kwargs['drawPoints']
    if 'drawEllipsoids' in kwargs:
        drawEllipsoids = kwargs['drawEllipsoids']
    if 'drawEllipsoidsSimplexTree' in kwargs:
        drawEllipsoidsSimplexTree = kwargs['drawEllipsoidsSimplexTree']
    if 'rPlot' in kwargs:
        rPlot = kwargs['rPlot']
    else: rPlot = 0.5
    if 'persistenceDim' in kwargs:
        persistenceDim = kwargs['persistenceDim']
    else:
        persistenceDim = 0


    # if set(listOfPlotPointsVars).issubset(kwargs): 
    if drawPoints or drawEllipsoids: 
        figsize = (14,7)
        if 'points' in kwargs:
            points = kwargs['points']
        else: 
            print("Cannot plot points; no list of points provided.")
            exit()

        dim = len(points[0])
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2,2)

        if dim == 2:
            axData = fig.add_subplot(gs[:, 0])
        else:
            axData = fig.add_subplot(gs[:,0], projection='3d')

        axBarE = fig.add_subplot(gs[0, 1])
        axBarR = fig.add_subplot(gs[1, 1])

        for point in points:
            axData.scatter(*point, c='k')
        axData.set_title('Data (%d points)' %len(points), fontsize=12)

        if drawEllipsoids:
            if 'expansionDim' in kwargs:
                expansionDim = kwargs['expansionDim']
            else: expansionDim = 1
            if ('ellipsoidList' in kwargs) or ('ellipseList' in kwargs):
                ellipsoidList = kwargs['ellipsoidList'] if 'ellipsoidList' in kwargs else kwargs['ellipseList']
            else:
                print('Cannot plot ellipsoids, no list of ellipsoids provided.')
        
            if dim == 2:
                plotEllipses(ellipsoidList, rPlot, axes=axData)
            if dim == 3:
                plotEllipsoids(ellipsoidList, rPlot, axes=axData)

            if drawEllipsoidsSimplexTree:
                if 'simplexTreeEllipsoids' in kwargs:
                    simplexTreeEllipsoids = kwargs['simplexTreeEllipsoids']
                    plotSimplexTree(points, simplexTreeEllipsoids, rPlot, axes=axData)
                    axData.set_title(f'Point cloud data for {len(points)} points')# and the ellipsoid simplex tree for r = %0.2f' %(rPlot), fontsize=12)
                else:
                    print('Cannot plot ellipsoid simplex tree; no ellipsoid simplex tree provided.')
                    axData.set_title(f'Point cloud data for {len(points)} points')
                    
        axData.set_aspect('equal', adjustable='box')

    else:
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2,1)
        axBarE = fig.add_subplot(gs[0, 0])
        axBarR = fig.add_subplot(gs[1, 0])


    if 'filename' in kwargs:
        filename = kwargs['filename']
    else: filename = 'data/plotTest.png'

    if plotDensity:
        highestColor = None

        birth_min = 0
        birth_max = 0
        death_min = 0
        death_max = 0
        cmap = 'Blues'
        nBarsDim0 = 0
        nBarsDim1 = 0
        nBarsDim2 = 0
        maxNBars = 10000
        if persistenceDim == 0:
            nBarsDim0 = maxNBars
        elif persistenceDim == 1:
            nBarsDim1 = maxNBars
        elif persistenceDim == 2:
            nBarsDim2 = maxNBars
        else:
            nBarsDim0 = maxNBars
            nBarsDim1 = maxNBars
            nBarsDim2 = maxNBars
        
        maxBarEndEllipsoids = 1
        maxBarEndRips = 1 
        barcodeEllipsoidsReduced, maxBarEndEllipsoids = reduceBarcode(barcodeEllipsoids, nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2)
        barcodeRipsReduced, maxBarEndRips = reduceBarcode(barcodeRips, nBarsDim0=nBarsDim0, nBarsDim1=nBarsDim1, nBarsDim2=nBarsDim2)
        #print(barcodeRipsReduced)
        for bar in barcodeRipsReduced:
            print(bar[1][1]-bar[1][0])
            # if bar[1][0] > bar[1][1]:
            #     print(bar)
        #exit()
        death_max = max(maxBarEndEllipsoids, maxBarEndRips)
        birth_max = death_max

    if plotDensity:
        print('Plotting ellipsoid density... ', end='', flush=True)
        #barcodeEllipsoids = totuple(barcodeEllipsoids)
        print(barcodeEllipsoids[0])
        #exit()
        plot_persistence_density(persistence=barcodeEllipsoidsReduced, axes=axBarE, fontsize=12,\
                                                    max_intervals=1000, #dimension=persistenceDim,
                                                    birth_min=birth_min, birth_max=birth_max,
                                                    death_min=death_min, death_max=death_max,
                                                    cmap=cmap,
                                                    highestColor=highestColor)
    else:
        print('Plotting ellipsoid barcode... ', end='', flush=True)
        plot_persistence_barcode(barcodeEllipsoids, inf_delta=0.5, axes=axBarE, fontsize=12,\
                                            axis_start = -0.1, infinity = xAxisEnd, max_intervals=100)
    print('Done.')
    axBarE.set_title('Ellipsoid barcode', fontsize=12)

    if plotDensity:
        print('Plotting Rips density... ', end='', flush=True)
        plot_persistence_density(persistence=barcodeRipsReduced, axes=axBarR, fontsize=12,\
                                                    max_intervals=1000, #dimension=persistenceDim,
                                                    birth_min=birth_min, birth_max=birth_max,
                                                    death_min=death_min, death_max=death_max,
                                                    cmap=cmap,
                                                    highestColor=highestColor)
    else:
        print('Plotting Rips barcode... ', end='', flush=True)
        plot_persistence_barcode(barcodeRips, inf_delta=0.5, axes=axBarR, fontsize=12,\
                                            axis_start = -0.1, infinity = xAxisEnd, max_intervals=100) #(0.1 + xAxisEnd))
    print('Done.')
    axBarR.set_title('Rips barcode', fontsize=12)
    
    if 'savePlot' in kwargs and kwargs['savePlot'] is True:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print('Plot saved to file.')

    if 'showPlot' in kwargs and kwargs['showPlot'] is True:
        plt.show()

def visualisationFromFile(\
        filename, \
        nBarsDim0=1, nBarsDim1=0, nBarsDim2=0, \
        rPlot=0.6, \
        drawEllipsoids=False, \
        drawEllipsoidsSimplexTree=False, \
        savePlot=False,
        showPlot=False,
        plotDensity=False,
        persistenceDim=0):

    print('Reading in the variables... ', end='', flush=True)
    vars = read_variables(filename)

    if 'barcodeEllipsoids' in vars:
        barcodeEllipsoids = vars['barcodeEllipsoids']
    elif 'barcode_ellipsoids' in vars:
        barcodeEllipsoids = vars['barcode_ellipsoids']
    else:
        print('Error: ellipsoids barcode not found.')
        exit()

    if 'barcodeRips' in vars:
        barcodeRips = vars['barcodeRips']
    elif 'barcode_rips' in vars:
        barcodeRips = vars['barcode_rips']
    else:
        print('Error: Rips barcode not found.')
        exit()

    print('Done.')

    print('Calculating the reduced barcodes... ', end='', flush=True)
    reducedBarcodeEllipsoids, maxBarEndEllipsoids = reduceBarcode( \
                                barcodeEllipsoids, \
                                nBarsDim0=nBarsDim0, \
                                nBarsDim1=nBarsDim1, \
                                nBarsDim2=nBarsDim2)
    reducedBarcodeRips, maxBarEndRips = reduceBarcode( \
                                barcodeRips, \
                                nBarsDim0=nBarsDim0, \
                                nBarsDim1=nBarsDim1, \
                                nBarsDim2=nBarsDim2)
    print('Done.')

    barcodeEllipsoids = reducedBarcodeEllipsoids
    barcodeRips = reducedBarcodeRips
    xAxisEnd = max(maxBarEndEllipsoids, maxBarEndRips) * 1.1 # TODO this won't work well if the filtration is negative

    print('Plotting...')

    filename = filename[:filename.rfind('.')] + '-barcodesDim=0-' + f'{nBarsDim0}' + '_1-' + f'{nBarsDim1}' + '_2-' + f'{nBarsDim2}' + datetime.now().strftime("_%Y%m%d_%H%M%S")
    filename = filename + '.png'
    
    if drawEllipsoids is True:
        simplexTreeEllipsoids = vars['simplexTreeEllipsoids']
        simplexTreeRips = vars['simplexTreeRips']
        points = vars['points']
        visualisation(points = points,\
                    ellipsoidList = vars['ellipsoidList'], \
                    rPlot = rPlot, \
                    simplexTreeEllipsoids = simplexTreeEllipsoids, \
                    barcodeEllipsoids = barcodeEllipsoids, \
                    barcodeRips = barcodeRips, \
                    xAxisEnd = xAxisEnd, \
                    showPlot = showPlot, \
                    savePlot = savePlot, \
                    filename = filename, \
                    drawEllipsoids = drawEllipsoids, 
                    drawEllipsoidsSimplexTree = drawEllipsoidsSimplexTree,
                    plotDensity = plotDensity,
                    persistenceDim = persistenceDim
                    )
    else:
        visualisation( \
                    xAxisEnd = xAxisEnd, \
                    barcodeEllipsoids = barcodeEllipsoids, \
                    barcodeRips = barcodeRips, \
                    showPlot = showPlot, \
                    savePlot = savePlot, \
                    filename = filename,
                    plotDensity = plotDensity,
                    persistenceDim = persistenceDim
                    )



# def plot_barcode(
#     persistence=[], # rename to 'barcode
#     alpha=0.6,
#     max_intervals=20000,
#     inf_delta=0.1,
#     legend=None,
#     colormap=None,
#     axes=None,
#     fontsize=16,
#     axis_start=None,
#     infinity=None,
#     bar_height=None,
# ):
 
#     try:
#         import matplotlib.pyplot as plt
#         import matplotlib.patches as mpatches
#         from matplotlib import rc

#         if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
#             plt.rc("text", usetex=True)
#             plt.rc("font", family="serif")
#         else:
#             plt.rc("text", usetex=False)
#             plt.rc("font", family="DejaVu Sans")

#         try:
#             persistence, nx2_array = _array_handler(persistence)
#             persistence = _limit_to_max_intervals(
#                 persistence, max_intervals, key=lambda life_time: life_time[1][1] - life_time[1][0]
#             )
#             (min_birth, max_death) = __min_birth_max_death(persistence)
#             persistence = sorted(persistence, key=lambda life_time: life_time[1][1] - life_time[1][0])
#             persistence = sorted(persistence, key=lambda birth: birth[1][0])
#         except IndexError:
#             min_birth, max_death = 0.0, 1.0
#             pass

#         delta = (max_death - min_birth) * inf_delta
#         # Replace infinity values with max_death + delta for bar code to be more
#         # readable
#         if infinity is None:
#             infinity = max_death + delta
#         if axis_start is None:
#             axis_start = min_birth - delta

#         if axes is None:
#             _, axes = plt.subplots(1, 1)
#         if colormap is None:
#             colormap = plt.cm.Set1.colors

#         x = [birth for (dim, (birth, death)) in persistence]
#         y = [(death - birth) if death != float("inf") else (infinity - birth) for (dim, (birth, death)) in persistence]
#         c = [colormap[dim] for (dim, (birth, death)) in persistence]

#         if bar_height is None:
#             bar_height = 0.6
     
#         axes.barh(range(len(x)), y, left=x, alpha=alpha, color=c, linewidth=0, height=bar_height)

#         if legend is None and not nx2_array:
#             # By default, if persistence is an array of (dimension, (birth, death)), display the legend
#             legend = True

#         if legend:
#             dimensions = {item[0] for item in persistence}
#             axes.legend(
#                 handles=[mpatches.Patch(color=colormap[dim], label=str(dim)) for dim in dimensions],
#                 loc="best",
#             )

#         axes.set_title("Persistence barcode", fontsize=fontsize)
#         axes.set_yticks([])

#         # -------------------------------- 
#         # the next part fixes the scaling 
#         margin = 0.4
#         n_bars = len(persistence)
#         padding = 1.5


#         fig = axes.get_figure()
#         axes.set_ylim([-padding, (n_bars-1) + padding])
#         fig_height = n_bars * (bar_height + margin)# + padding
#         axes.invert_yaxis() # temp changing this
        
#         fig.set_figheight(fig_height)

#         width_inches = 10
#         height_inches = n_bars * margin + padding
#         fig.set_size_inches(width_inches, height_inches)
#         # fig.subplots_adjust(top=0.85)
#         fig.subplots_adjust(top=1)
#         # --------------------------------          


#         # Ends plot on infinity value and starts a little bit before min_birth
#         if len(x) != 0:
#             axes.set_xlim((axis_start, infinity))
#         return axes

