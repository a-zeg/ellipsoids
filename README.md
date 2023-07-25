# Ellipsoids

To run the code, run ellipsoids.py. Below is a more detailed description of all the files in this project.

## ellipsoids.py
This is the main file. It:
- calculates the ellipsoid simplex tree of the data specified in the variable `points` in the function `main()` at the filtration levels specified by `rValues`;
- calculates the Rips simplex tree at the same filtration level;
- creates a plot containing the initial data, the ellipsoid simplex tree and the ellipses used to create it (both at filtration level `rPlot`), and the corresponding barcodes.

## barcodePlotting.py
This is a slightly modified version of https://gudhi.inria.fr/python/latest/_modules/gudhi/persistence_graphical_tools.html#plot_persistence_barcode, because it was not possible to specify the start and the end of x-axis in the original.

## figure_eight.py
Bastian's file.

## shapes.py
Slightly modified version of Bastian's code https://github.com/aidos-lab/pytorch-topological/blob/main/torch_topological/data/shapes.py. The original code outputs torch tensor, and I just needed a numpy array.