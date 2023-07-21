"""Create figure eight data set.

The purpose of this script is to create a figure eight data set with
different 'neck widths' to illustrate that many data sets require an
analysis at multiple scales.
"""

import argparse
import sys

import numpy as np


def figure_eight(n, a, b):
    """Sample a set of points from a figure eight curve.

    Parameters
    ----------
    n : int
        Number of points to sample

    a : float
        Controls extents of the curve. A larger `a` parameter will
        result in larger scaling.

    b : float
        Controls neck size of the curve. A larger `b` parameter will
        result in an increased neck size.

    Returns
    -------
    np.array
        Array of shape (n, 2). Will contain the sampled points.
    """
    T = np.linspace(-5, 5, num=n)

    X = a * np.sin(T)
    Y = a * np.sin(T)**2 * np.cos(T) + b * np.cos(T)

    X = np.column_stack((X, Y))
    X += np.random.default_rng().uniform(0.05, 0.10, size=(n, 2))
    return X


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--num-points', '-n',
        default=100,
        type=int,
        help='Number of points'
    )
    parser.add_argument(
        '-a',
        default=1.0,
        type=float,
        help='Curve scale parameter'
    )
    parser.add_argument(
        '-b',
        default=0.0,
        type=float,
        help='Curve neck width parameter'
    )

    args = parser.parse_args()

    X = figure_eight(args.num_points, args.a, args.b)
    np.savetxt(sys.stdout, X, fmt='%.04f')
