<div align="center">

  <h1 align="center">Ellipsoids project</h3>

  <p align="center">
    Calculate persistence homology of ellipsoid complexes.
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li><a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

# About the project
![Ellipsoid plots](images/ellipsoids_nPts=10_rStep=0.01_nbhdSize=3_20230731_115042.png)

The main goal of this project is to investigate properties of ellipsoid complexes. Ellipsoid complexes are Rips-like complexes created by considering ellipsoids centered at the data points (as opposed to balls in the case of Rips complexes).

This software generates ellipsoid complexes from the given point cloud data and calculates barcodes of their persistent homology.

# Getting started
## Prerequisites
- [GUDHI](https://gudhi.inria.fr/index.html) Python interface 3.8.0
- matplotlib 3.7.1
- numpy 1.25.0
- scikit_learn 1.3.0
- scipy 1.11.1

See also [requirements.txt](./requirements.txt).

## Installation
There is no need to install anything else, just download the project.


# Usage
To run the programme, run `ellipsoids.py`. 
Depending on the value of `boolPlotFromFile`, the programme will either:

1. Load data from a file specified in the variable `filenameLoad` and plot it; or
2. Calculate barcodes of ellipsoid and Rips complexes from data specified in `points`. 

## User input
All the parameters need to be specified within this file, in places marked as follows:
```python
###### User input ######
...
########################
```

In particular, the user is expected to specify the following:
| Variable name | Possible values | Description | 
| :----------- | :-------------- | :--------- | 
| `boolPlotFromFile` | `True` <br> `False` | Load data from file. <br> Perform calculations. |
| `filenameLoad` | String | Name of the JSON file from which to load data.|
| `rPlot` | Non-negative float | Filtration for which to plot the ellipsoid complex. If `rPlot = 0`, the ellipsoid complex will not be plotted.
| `boolSaveData` | Bool | If `True`, the data will be saved to a file. |
| `boolShowPlot` | Bool | If `True`, the plot will be shown. |
| `boolSavePlot` | Bool | If `True`, the plot will be saved as a PNG. |
| `dim`          | Integer | Dimension of the data. |
| `rValues` | Numpy array of non-negative floats | Filtration values. |
| `nbhdSize` | Positive integer | Number of points for doing PCA |
| `points` | Numpy array of numpy arrays | An array of data points |

## Output
Depending on the user input (see the <a href="#user-input">User input</a> section above), the programme will:
- Output a plot containing the data points, ellipsoid complex barcode, Rips barcode.
- Output a plot containing the data points with the ellipsoid complex, ellipsoid complex barcode, Rips barcode.
- Save the generated plot.
- Save the generated data to a JSON file.


# Roadmap
- [x] Minimal working example 
- [x] Save generated data to file.
- [x] Import saved data from a file and generate figures.
- [ ] Support for 3+ dimensional data.

# Acknowledgements
- Kalisnik, Lesnik - [Finding the Homology of Manifolds using Ellipsoids (2020)](https://arxiv.org/abs/2006.09194) is the inspiration for this project.
- [GUDHI Library](https://gudhi.inria.fr/index.html) is used throughout the project. The file `barcodePlotting.py` is a slight modification of [this file](https://gudhi.inria.fr/python/latest/_modules/gudhi/persistence_graphical_tools.html#plot_persistence_barcode). <The modification makes it possible to specify the start and the end of x-axis.>
- The file `figure_eight.py` is courtesy of Bastian Rieck.
- The file `shapes.py` comes from [pytorch-topological](https://github.com/aidos-lab/pytorch-topological/blob/main/torch_topological/data/shapes.py).



