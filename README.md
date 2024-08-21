<div align="center">

  <h1 align="center">Ellipsoids project</h3>

  <p align="center">
    Persistence homology of ellipsoid complexes
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
       <ul>
        <li><a href="#user-input">To be updated</a></li>
      </ul>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

# About the project
![Ellipsoid plots](images/ellipsoids_nPts=10_rStep=0.01_nbhdSize=3_20230731_115042.png)

The main goal of this project is to investigate properties of ellipsoid complexes. Ellipsoid complexes are Rips-like complexes created by considering ellipsoids centered at the data points (as opposed to balls in the case of Rips complexes).

This software generates ellipsoid complexes from the given point cloud data and calculates barcodes of their persistent homology.

# Getting started
## Prerequisites

TO BE UPDATED 

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
The programme will calculate barcodes of ellipsoid and Rips complexes from data specified in `points` and possibly save the results to a JSON file and plot the barcodes, data, and possibly the ellipsoids.

To plot results from a JSON file, run `plotFromFile.py`.


## TO BE UPDATED



# Acknowledgements

TO BE UPDATED

- Kalisnik, Lesnik - [Finding the Homology of Manifolds using Ellipsoids (2020)](https://arxiv.org/abs/2006.09194) is the inspiration for this project.
- [GUDHI Library](https://gudhi.inria.fr/index.html) is used throughout the project. The file `barcodePlotting.py` is a slight modification of [this file](https://gudhi.inria.fr/python/latest/_modules/gudhi/persistence_graphical_tools.html#plot_persistence_barcode). <The modification makes it possible to specify the start and the end of x-axis.>
- The file `figure_eight.py` is courtesy of Bastian Rieck.
- The file `shapes.py` comes from [pytorch-topological](https://github.com/aidos-lab/pytorch-topological/blob/main/torch_topological/data/shapes.py).



