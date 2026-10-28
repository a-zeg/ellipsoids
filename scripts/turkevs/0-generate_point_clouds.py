#!/usr/bin/env python3
"""Generate point clouds and save directly in the Turkevs dataset format.

Run from the project root, as with 0-generate_turkevs_datasets.py.
Requires NumPy, SciPy, and the project's ellipsoids package. No input NPZ or
other point-cloud-generation script is needed.
"""

import logging
import math
from pathlib import Path
import os
import sys

import numpy as np
from scipy.spatial.distance import pdist

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import Dataset, TurkevsDatasetInfo, TurkevsTransformation
from ellipsoids.turkevs.config import TURKEVS_DATA_DIR, get_datasets_path, get_setup_name
from ellipsoids.turkevs.turkevs_utils import add_dataset_id, generate_datasets_hash
from ellipsoids.turkevs.data_construction import calculate_point_clouds_under_trnsf
from ellipsoids.data_handling import save_to_json
from ellipsoids.logging_setup import setup_logging

N_SHAPE_VARIATIONS = 5
N_POINTS = 600
SEED = 0  # Controls original, bottleneck, torus, and transformation sampling.
TRANSFORMATIONS_MODE = "all"  # "all" or "standard"
# TRANSFORMATIONS_MODE = "standard"  # "all" or "standard"

_TAU = 2.0 * np.pi

def _original_parabola(n=600):
    x = np.random.uniform(-1, 1, size=n)
    y = x ** 2
    X = np.array([x, y]).T
    return X

def _original_circle(n=600, r=1):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    x = r * np.cos(t)
    y = r * np.sin(t)
    X = np.array([x, y]).T
    return X

def _original_closed_curve(n=600, num_vertices=3):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    b = 1
    k = num_vertices
    a = k * b
    x = (a - b) * np.cos(t) + b * np.cos(t * (a / b - 1))
    y = (a - b) * np.sin(t) - b * np.sin(t * (a / b - 1))
    X = np.array([x, y]).T
    return X

def _original_annulus(n=600):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    x = (1 + np.random.uniform(0, 1, size=n)) * np.cos(t)
    y = (1 + np.random.uniform(0, 1, size=n)) * np.sin(t)
    X = np.array([x, y]).T
    return X

def _original_two_circles(n=600):
    n1 = int(n / 2)
    n2 = n - n1
    circle1 = _original_circle(n1, r=1)
    circle2 = 2 + _original_circle(n2, r=1)
    X = np.concatenate([circle1, circle2])
    return X

def _original_lemniscate(n=600):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    x = np.cos(t)
    y = np.sin(2 * t)
    X = np.array([x, y]).T
    return X

def _original_concentric_circles(n=600):
    n1 = int(n / 2)
    n2 = n - n1
    circle1 = _original_circle(n1, r=1)
    circle2 = _original_circle(n2, r=0.5)
    X = np.concatenate([circle1, circle2])
    return X

def _original_rose_curve(n=600, a=2):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    x = np.cos(a * t) * np.cos(t)
    y = np.cos(a * t) * np.sin(t)
    X = np.array([x, y]).T
    return X

def _original_rose_nine_petals(n=600):
    n_rose_three_petals = int(n / 3)
    rose1 = _original_rose_curve(n_rose_three_petals, a=3)
    rose2 = _original_rotate(rose1, t=np.pi / 4.5)
    rose3 = _original_rotate(rose2, t=np.pi / 4.5)
    n_rmn = n - 3 * n_rose_three_petals
    X = np.concatenate([rose1, rose2, rose3])
    if n_rmn > 0:
        epsilon = 1e-06
        center_noise = np.random.normal(scale=epsilon, size=(n_rmn, 2))
        center = np.zeros((n_rmn, 2)) + center_noise
        X = np.concatenate([X, center])
    return X

def _original_eyes_1(n=600):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    a = 3
    b = 1
    c = 3
    d = 1
    x = np.cos(a * t) - np.cos(b * t) ** 3
    y = np.sin(c * t) - np.sin(d * t) ** 3
    X = np.array([x, y]).T
    return X

def _original_olympics(n=600):
    t = np.random.uniform(0, 2 * np.pi, size=int(n / 5))
    circle = np.array([np.cos(t), np.sin(t)])
    l = 1
    center1 = np.array([0, 0.75])
    center2 = np.array([-1 * l, 0])
    center3 = np.array([-2 * l, 0.75])
    center4 = np.array([1 * l, 0])
    center5 = np.array([2 * l, 0.75])
    X1 = np.transpose(np.transpose(circle) + center1)
    X2 = np.transpose(np.transpose(circle) + center2)
    X3 = np.transpose(np.transpose(circle) + center3)
    X4 = np.transpose(np.transpose(circle) + center4)
    X5 = np.transpose(np.transpose(circle) + center5)
    X = np.concatenate([X1, X2, X3, X4, X5], 1)
    X = X.T
    return X

def _original_square(n=600):
    n_line = int(n / 4)
    n_last_line = n - 3 * n_line
    x = np.random.uniform(0, 1, size=n_line)
    x1 = x
    x2 = x
    x3 = [0] * n_line
    x4 = [1] * n_last_line
    y1 = [0] * n_line
    y2 = [1] * n_line
    y3 = np.random.uniform(0, 1, size=n_line)
    y4 = np.random.uniform(0, 1, size=n_last_line)
    x = np.concatenate([x1, x2, x3, x4])
    y = np.concatenate([y1, y2, y3, y4])
    X = np.array([x, y]).T
    return X

def _original_nine_squares(n=600):
    n_square = int(n / 9)
    X1 = _original_square(n_square)
    X_shift = np.zeros(X1.shape)
    X_shift[:, 0] = 1.1
    Y_shift = np.zeros(X1.shape)
    Y_shift[:, 1] = 1.1
    X2 = X1 + X_shift
    X3 = X2 + X_shift
    X4 = X1 + Y_shift
    X5 = X4 + X_shift
    X6 = X5 + X_shift
    X7 = X4 + Y_shift
    X8 = X7 + X_shift
    X9 = X8 + X_shift
    n_first_square = n - 8 * n_square
    X1 = _original_square(n_first_square)
    X = np.concatenate([X1, X2, X3, X4, X5, X6, X7, X8, X9])
    return X

def _original_sphere(n=600, r=1):
    theta = np.random.uniform(0, 2 * np.pi, size=n)
    phi = np.arccos(np.random.uniform(-1, 1, size=n))
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    X = np.array([x, y, z]).T
    return X

def _original_fountain(n=600):
    theta = np.random.uniform(0, 2 * np.pi, size=n)
    phi = np.arccos(np.random.uniform(-1, 1, size=n))
    x = (1 + 0.4 * np.sin(7 * phi) * np.cos(theta)) * np.cos(theta) * np.sin(phi)
    y = (1 + 0.4 * np.sin(7 * phi) * np.cos(theta)) * np.sin(theta) * np.sin(phi)
    z = (1 + 0.4 * np.sin(7 * phi) * np.cos(theta)) * np.cos(phi)
    X = np.array([x, y, z]).T
    return X

def _original_two_spheres(n=600):
    n1 = int(n / 2)
    n2 = n - n1
    sphere1 = _original_sphere(n1, r=1)
    sphere2 = 2 + _original_sphere(n2, r=1)
    X = np.concatenate([sphere1, sphere2])
    return X

def _original_3d_curve(n=600, gamma=1):
    t = np.random.uniform(0, 2 * np.pi, size=n)
    Xdata = np.cos(t)
    Ydata = np.sin(t)
    XXnormal = gamma * np.cos(t) ** 2
    YYnormal = gamma * np.sin(t) ** 2
    XYnormal = gamma * np.cos(t) * np.sin(t)
    X = np.vstack((Xdata, Ydata, XXnormal, XYnormal, XYnormal, YYnormal)).T
    X = np.array([X[:, 0], X[:, 1], X[:, 2]]).T
    return X

def _original_torus(n=600, r1=1, r2=0.5):
    theta = np.random.uniform(0, 2 * np.pi, size=n)
    phi = np.random.uniform(0, 2 * np.pi, size=n)
    x = (r1 + r2 * np.cos(theta)) * np.cos(phi)
    y = (r1 + r2 * np.cos(theta)) * np.sin(phi)
    z = r2 * np.sin(theta)
    X = np.array([x, y, z]).T
    return X

def _original_double_torus(n=600, r1=1, r2=0.5):
    n1 = int(n / 2)
    n2 = n - n1
    torus1 = _original_torus(n1, r1, r2)
    torus2 = 1 + 0.2 + _original_torus(n2, r1, r2)
    X = np.concatenate([torus1, torus2])
    return X

def _original_necklace(n=600):
    X1 = _original_sphere(int(n / 4)) + [2, 0, 0]
    X2 = _original_sphere(int(n / 4)) + [-1, 2 * 0.866, 0]
    X3 = _original_sphere(int(n / 4)) + [-1, -2 * 0.866, 0]
    X4 = 2 * _original_circle(int(n / 4))
    X4 = np.array((X4[:, 0], X4[:, 1], np.zeros(int(n / 4)))).T
    X = np.concatenate((X1, X2, X3, X4))
    return X

def _original_handlebody(n=600, r=1, r1=0.5, r2=0.25):
    n_torus = int(n / 4)
    torus = _original_torus(n_torus, r1, r2)
    torus1 = np.asarray([r, 0, 0]) + torus
    handle1 = []
    for p in range(n_torus):
        x = torus1[p, 0]
        y = torus1[p, 1]
        z = torus1[p, 2]
        if x ** 2 + y ** 2 + z ** 2 > r ** 2:
            handle1.append([x, y, z])
    handle1 = np.asarray(handle1).reshape(-1, 3)
    n_handle1 = handle1.shape[0]
    torus2 = np.asarray([-r, 0, 0]) + torus
    handle2 = []
    for p in range(n_torus):
        x = torus2[p, 0]
        y = torus2[p, 1]
        z = torus2[p, 2]
        if x ** 2 + y ** 2 + z ** 2 > r ** 2:
            handle2.append([x, y, z])
    handle2 = np.asarray(handle2).reshape(-1, 3)
    n_handle2 = handle2.shape[0]
    torus3 = np.asarray([0, r, 0]) + torus
    handle3 = []
    for p in range(n_torus):
        x = torus3[p, 0]
        y = torus3[p, 1]
        z = torus3[p, 2]
        if x ** 2 + y ** 2 + z ** 2 > r ** 2:
            handle3.append([x, y, z])
    handle3 = np.asarray(handle3).reshape(-1, 3)
    n_handle3 = handle3.shape[0]
    torus4 = np.asarray([0, -r, 0]) + torus
    handle4 = []
    for p in range(n_torus):
        x = torus4[p, 0]
        y = torus4[p, 1]
        z = torus4[p, 2]
        if x ** 2 + y ** 2 + z ** 2 > r ** 2:
            handle4.append([x, y, z])
    handle4 = np.asarray(handle4).reshape(-1, 3)
    n_handle4 = handle4.shape[0]
    n_rmn = n - n_handle1 - n_handle2 - n_handle3 - n_handle4
    n_sphere = int(0.9 * n_rmn)
    n_circle = n_rmn - n_sphere
    sphere = _original_sphere(n_sphere, r)
    circle_2d = _original_circle(n_circle, r)
    circle = np.zeros((n_circle, 3))
    for p in range(n_circle):
        circle[p] = np.asarray([circle_2d[p][0], circle_2d[p][1], 2 * r])
    X = np.concatenate([circle, sphere, handle1, handle2, handle3, handle4])
    return X

def _original_normalize(X):
    dim = X.shape[1]
    for i in range(dim):
        X[:, i] = (X[:, i] - np.min(X[:, i])) / (np.max(X[:, i]) - np.min(X[:, i]))
    return X

def _original_transform(X, trnsf_matrix):
    X_trnsf = np.zeros(X.shape)
    for i in range(X.shape[0]):
        X_trnsf[i] = np.matmul(trnsf_matrix, X[i])
    return X_trnsf

def _original_rotate(X, t=None):
    if t is None:
        t = np.random.uniform(-np.pi / 9, np.pi / 9)
    dim = X.shape[1]
    if dim == 2:
        rotation_matrix = np.asarray([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    if dim == 3:
        rotation_matrix = np.asarray([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
    X = _original_transform(X, rotation_matrix)
    return X

def _orthogonal_union_polygon(rectangles):
    """Return the CCW vertices of a simple, hole-free rectangle union.

    ``rectangles`` contains (xmin, ymin, xmax, ymax). Each rectangle must have
    strictly positive width and height. Boundary edges are traced on the grid
    of rectangle coordinates, then collinear intermediate vertices are removed.
    """
    rects = np.asarray(rectangles, dtype=float)
    if rects.ndim != 2 or rects.shape[1] != 4 or (not len(rects)):
        raise ValueError('rectangles must be a nonempty (m,4) array')
    if not np.isfinite(rects).all() or np.any(rects[:, 2:] <= rects[:, :2]):
        raise ValueError('each rectangle must have finite, increasing bounds')
    xs = np.unique(rects[:, [0, 2]])
    ys = np.unique(rects[:, [1, 3]])
    xc = (xs[:-1] + xs[1:]) / 2
    yc = (ys[:-1] + ys[1:]) / 2
    occupied = np.zeros((len(xc), len(yc)), dtype=bool)
    for x0, y0, x1, y1 in rects:
        occupied |= (xc[:, None] > x0) & (xc[:, None] < x1) & (yc[None, :] > y0) & (yc[None, :] < y1)
    edges = []
    nx, ny = occupied.shape
    for i, j in np.argwhere(occupied):
        i, j = (int(i), int(j))
        if j == 0 or not occupied[i, j - 1]:
            edges.append(((i, j), (i + 1, j)))
        if i == nx - 1 or not occupied[i + 1, j]:
            edges.append(((i + 1, j), (i + 1, j + 1)))
        if j == ny - 1 or not occupied[i, j + 1]:
            edges.append(((i + 1, j + 1), (i, j + 1)))
        if i == 0 or not occupied[i - 1, j]:
            edges.append(((i, j + 1), (i, j)))
    outgoing, incoming = ({}, {})
    for a, b in edges:
        if a in outgoing or b in incoming:
            raise ValueError('union has a nonsimple vertex-only contact')
        outgoing[a] = b
        incoming[b] = a
    start = min(outgoing)
    chain, vertex = ([], start)
    while True:
        chain.append(vertex)
        vertex = outgoing[vertex]
        if vertex == start:
            break
        if len(chain) > len(edges):
            raise ValueError('failed to close rectangle-union boundary')
    if len(chain) != len(edges):
        raise ValueError('union is disconnected or has a hole')
    polygon = np.array([[xs[i], ys[j]] for i, j in chain])
    d0 = polygon - np.roll(polygon, 1, axis=0)
    d1 = np.roll(polygon, -1, axis=0) - polygon
    cross = d0[:, 0] * d1[:, 1] - d0[:, 1] * d1[:, 0]
    polygon = polygon[cross != 0]
    area2 = np.sum(polygon[:, 0] * np.roll(polygon[:, 1], -1) - polygon[:, 1] * np.roll(polygon[:, 0], -1))
    if area2 <= 0:
        raise ValueError('expected a counterclockwise outer boundary')
    return polygon

def _rounded_union(rectangles, radius):
    """Replace every corner of a rectangle union by an exact circular fillet.

    ``radius`` is either a positive scalar or a dictionary with positive
    ``convex`` and ``concave`` values. Adjacent fillets may touch but may not
    overlap along an original edge. This local check cannot exclude nonlocal
    collisions in arbitrary unions: the intended inputs are chamber-and-corridor
    trees with separated chambers. Validate those geometric clearances upstream.
    """
    vertices = _orthogonal_union_polygon(rectangles)
    n = len(vertices)
    incoming = vertices - np.roll(vertices, 1, axis=0)
    outgoing = np.roll(vertices, -1, axis=0) - vertices
    edge_lengths = np.linalg.norm(outgoing, axis=1)
    incoming /= np.linalg.norm(incoming, axis=1)[:, None]
    outgoing /= edge_lengths[:, None]
    turns = np.sign(incoming[:, 0] * outgoing[:, 1] - incoming[:, 1] * outgoing[:, 0])
    if isinstance(radius, dict):
        radii = np.where(turns > 0, radius['convex'], radius['concave']).astype(float)
    else:
        radii = np.full(n, float(radius))
    if not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError('fillet radii must be finite and positive')
    tol = 1e-12 * max(1.0, np.max(np.abs(vertices)))
    if np.any(radii + np.roll(radii, -1) > edge_lengths + tol):
        raise ValueError('adjacent corner fillets overlap along an edge')
    starts = vertices - radii[:, None] * incoming
    ends = vertices + radii[:, None] * outgoing
    left_normals = np.column_stack((-incoming[:, 1], incoming[:, 0]))
    centers = starts + (turns * radii)[:, None] * left_normals
    primitives = []
    for i in range(n):
        delta = starts[i] - centers[i]
        primitives.append({'type': 'arc', 'center': centers[i].tolist(), 'radius': float(radii[i]), 'start_angle': float(np.arctan2(delta[1], delta[0])), 'sweep': float(turns[i] * math.pi / 2), 'start': starts[i].tolist(), 'end': ends[i].tolist()})
        j = (i + 1) % n
        if np.linalg.norm(starts[j] - ends[i]) > tol:
            primitives.append({'type': 'line', 'start': ends[i].tolist(), 'end': starts[j].tolist()})
    return primitives

def _primitive_length(primitive):
    if primitive['type'] == 'line':
        return float(np.linalg.norm(np.asarray(primitive['end']) - primitive['start']))
    if primitive['type'] == 'arc':
        return float(primitive['radius'] * abs(primitive['sweep']))
    raise ValueError('unknown primitive type')

def _boundary_length(primitives):
    return float(sum((_primitive_length(p) for p in primitives)))

def _boundary_at_arclength(primitives, s, return_tangent=False, return_indices=False):
    """Evaluate periodic arclength coordinates, optionally with unit tangents.

    Output points have shape ``np.shape(s)+(2,)``. If requested, tangents have
    the same shape and primitive indices have shape ``np.shape(s)``. Returns
    points alone or a tuple (points, tangents?, indices?) in that order.
    """
    s = np.asarray(s, dtype=float)
    lengths = np.array([_primitive_length(p) for p in primitives])
    ends = np.cumsum(lengths)
    if len(ends) == 0 or ends[-1] <= 0:
        raise ValueError('empty or zero-length boundary')
    flat = np.mod(s.ravel(), ends[-1])
    indices = np.minimum(np.searchsorted(ends, flat, side='right'), len(ends) - 1)
    offsets = np.r_[0.0, ends[:-1]]
    local = flat - offsets[indices]
    points, tangents = (np.empty((len(flat), 2)), np.empty((len(flat), 2)))
    for i, p in enumerate(primitives):
        mask = indices == i
        if not mask.any():
            continue
        if p['type'] == 'line':
            direction = (np.asarray(p['end']) - p['start']) / lengths[i]
            points[mask] = p['start'] + local[mask, None] * direction
            tangents[mask] = direction
        else:
            direction = np.sign(p['sweep'])
            angle = p['start_angle'] + direction * local[mask] / p['radius']
            radial = np.column_stack((np.cos(angle), np.sin(angle)))
            points[mask] = p['center'] + p['radius'] * radial
            tangents[mask] = direction * np.column_stack((-radial[:, 1], radial[:, 0]))
    output = [points.reshape(s.shape + (2,))]
    if return_tangent:
        output.append(tangents.reshape(s.shape + (2,)))
    if return_indices:
        output.append(indices.reshape(s.shape))
    return output[0] if len(output) == 1 else tuple(output)

def _sample_boundary(primitives, n, rng, return_tangent=False, return_indices=False):
    """Draw iid samples uniformly with respect to exact boundary arclength."""
    s = rng.uniform(0, _boundary_length(primitives), size=n)
    return _boundary_at_arclength(primitives, s, return_tangent, return_indices)

def _distance_to_primitive(points, primitive):
    """Exact Euclidean distances to a closed straight segment or circular arc."""
    points = np.asarray(points, dtype=float)
    if points.shape[-1] != 2:
        raise ValueError('points must have final dimension 2')
    if primitive['type'] == 'line':
        a, b = (np.asarray(primitive['start']), np.asarray(primitive['end']))
        direction = b - a
        t = np.clip(np.sum((points - a) * direction, axis=-1) / np.dot(direction, direction), 0, 1)
        return np.linalg.norm(points - (a + t[..., None] * direction), axis=-1)
    c, r = (np.asarray(primitive['center']), primitive['radius'])
    delta = points - c
    theta = np.arctan2(delta[..., 1], delta[..., 0])
    a0, sweep = (primitive['start_angle'], primitive['sweep'])
    if sweep > 0:
        inside_angle = np.mod(theta - a0, _TAU) <= sweep + 1e-13
    else:
        inside_angle = np.mod(a0 - theta, _TAU) <= -sweep + 1e-13
    start = c + r * np.array([math.cos(a0), math.sin(a0)])
    end = c + r * np.array([math.cos(a0 + sweep), math.sin(a0 + sweep)])
    endpoints = np.minimum(np.linalg.norm(points - start, axis=-1), np.linalg.norm(points - end, axis=-1))
    radial = np.abs(np.linalg.norm(delta, axis=-1) - r)
    return np.where(inside_angle, radial, endpoints)

def _distances_to_boundary(points, primitives):
    """Distances to each primitive, with final axis indexing the primitives."""
    return np.stack([_distance_to_primitive(points, p) for p in primitives], axis=-1)

def _boundary_distance(points, primitives):
    return np.min(_distances_to_boundary(points, primitives), axis=-1)

def _contains_points(points, primitives, include_boundary=True, tolerance=1e-12):
    """Analytic even-odd ray containment for orthogonal filleted boundaries.

    Every quarter-circle arc is monotone in y. A horizontal ray therefore has
    at most one intersection with each primitive. The half-open y convention
    avoids double-counting their shared endpoints. Boundary points are included
    by default (to the supplied absolute floating-point tolerance).
    """
    points = np.asarray(points, dtype=float)
    shape = points.shape[:-1]
    points = points.reshape(-1, 2)
    x, y = (points[:, 0], points[:, 1])
    inside = np.zeros(len(points), dtype=bool)
    for p in primitives:
        a, b = (np.asarray(p['start']), np.asarray(p['end']))
        ymin, ymax = (min(a[1], b[1]), max(a[1], b[1]))
        active = (y >= ymin) & (y < ymax)
        if not active.any():
            continue
        if p['type'] == 'line':
            xx = a[0] + (y[active] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        else:
            cx, cy = p['center']
            sign = 1.0 if math.cos(p['start_angle'] + p['sweep'] / 2) >= 0 else -1.0
            xx = cx + sign * np.sqrt(np.maximum(0.0, p['radius'] ** 2 - (y[active] - cy) ** 2))
        inside[active] ^= xx > x[active]
    if include_boundary and len(points):
        for start in range(0, len(points), 4096):
            batch = points[start:start + 4096]
            inside[start:start + 4096] |= _boundary_distance(batch, primitives) <= tolerance
    return inside.reshape(shape)

def _chamber_layout(chambers, theme, spacing=3.2):
    """Place chambers along a straight or branched tree of narrow corridors."""
    b = int(chambers)
    if b == 1:
        return np.array([[0., 0.]]), []
    if theme == "straight":
        nodes = np.column_stack((np.arange(b, dtype=float), np.zeros(b)))
        segments = [(nodes[j].copy(), nodes[j+1].copy()) for j in range(b-1)]
    elif theme == "branched":
        if b == 2:
            nodes = np.array([[-1., 0.], [1., 0.]])
            segments = [(nodes[0], np.zeros(2)), (np.zeros(2), nodes[1]),
                        (np.zeros(2), np.array([0., 0.65]))]
        elif b == 3:
            nodes = np.array([[-1., 0.], [1., 0.], [0., 1.]])
            segments = [(np.zeros(2), p.copy()) for p in nodes]
        else:
            spine = int(np.ceil(b / 3))
            base = [np.array([float(j), 0.]) for j in range(spine)]
            segments = [(base[j].copy(), base[j+1].copy())
                        for j in range(spine-1)]
            nodes = list(base)
            for j in range(spine):
                for side in (-1., 1.):
                    if len(nodes) < b:
                        tip = np.array([float(j), side])
                        nodes.append(tip)
                        segments.append((base[j].copy(), tip.copy()))
            nodes = np.asarray(nodes)
    else:
        raise ValueError(f"Unknown bottleneck theme {theme!r}")
    return (nodes * spacing,
            [(np.asarray(a)*spacing, np.asarray(b)*spacing) for a, b in segments])

def _rectangle(cx, cy, rx, ry):
    # Rounding prevents spurious 1e-15-wide cells at shared corridor edges.
    return np.round([cx-rx, cy-ry, cx+rx, cy+ry], 12).tolist()

def _bottleneck_component(chambers, theme, rng):
    """One rounded, contractible chamber-and-corridor region."""
    gap = 0.85
    nodes, segments = _chamber_layout(chambers, theme)
    half_sizes = rng.uniform(0.90, 1.10, size=(len(nodes), 2))
    rectangles = [_rectangle(*p, *r) for p, r in zip(nodes, half_sizes)]
    for a, b in segments:
        center = (a+b)/2
        radii = np.abs(b-a)/2 + gap/2
        rectangles.append(_rectangle(*center, *radii))
    primitives = _rounded_union(
        rectangles, radius={"convex": min(0.35, .45*gap), "concave": 0.10}
    )
    bounds = [min(r[0] for r in rectangles), min(r[1] for r in rectangles),
              max(r[2] for r in rectangles), max(r[3] for r in rectangles)]
    return dict(primitives=primitives, bounds=bounds,
                boundary_length=float(_boundary_length(primitives)))

def _allocate_component_points(lengths, n):
    """Allocate integer sample counts proportionally to boundary length."""
    expected = np.asarray(lengths) * n / np.sum(lengths)
    counts = np.floor(expected).astype(int)
    counts[np.argsort(-(expected-counts), kind="stable")[:n-int(counts.sum())]] += 1
    return counts

def _sample_filled_component(component, n, rng):
    """Uniform area samples by rejection from the bounding rectangle."""
    xmin, ymin, xmax, ymax = component["bounds"]
    accepted, have = [], 0
    while have < n:
        proposal = rng.uniform([xmin, ymin], [xmax, ymax], size=(max(2048, 3*n), 2))
        good = proposal[_contains_points(proposal, component["primitives"])]
        accepted.append(good)
        have += len(good)
    return np.concatenate(accepted)[:n]

def _sample_bottleneck(beta, theme, seed):
    """Generate one planar cloud: filled region if beta=0, boundaries otherwise.

    The positive-label clouds consist of beta disjoint simple closed curves.
    Chamber counts per component are 12; 6+6; or 6+2+2+2, respectively.
    """
    rng = np.random.default_rng(seed)
    partition = {0: [12], 1: [12], 2: [6, 6], 4: [6, 2, 2, 2]}[beta]
    components = [_bottleneck_component(b, theme, rng) for b in partition]
    counts = _allocate_component_points([c["boundary_length"] for c in components], N_POINTS)
    points, ids = [], []
    offset = 0.
    for j, (component, count) in enumerate(zip(components, counts)):
        xmin, ymin, xmax, ymax = component["bounds"]
        translation = np.array([offset-xmin, -(ymin+ymax)/2])
        offset += xmax-xmin + 2.0
        if beta == 0:
            xy = _sample_filled_component(component, int(count), rng)
        else:
            xy = _sample_boundary(component["primitives"], int(count), rng)
        points.append(xy + translation)
        ids.extend([j]*int(count))
    raw = np.concatenate(points)
    theta = rng.uniform(0., 2*np.pi)
    rotation = np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta), np.cos(theta)]])
    raw = raw @ rotation.T
    clean_center = raw.mean(axis=0)
    clean_scale = 2. / float(pdist(raw).max())
    scaled_clean = (raw-clean_center)*clean_scale
    observed = scaled_clean + rng.normal(scale=0.001, size=raw.shape)
    observed_center = observed.mean(axis=0)
    observed_scale = 2. / float(pdist(observed).max())
    observed = (observed-observed_center)*observed_scale
    clean = (scaled_clean-observed_center)*observed_scale
    padded = np.zeros((N_POINTS, 3)); padded[:, :2] = observed
    padded_clean = np.zeros((N_POINTS, 3)); padded_clean[:, :2] = clean
    return dict(points=padded, clean=padded_clean,
                component_ids=np.asarray(ids, dtype=np.int16))

def _oval_surface(R, r, stretch, u, v):
    """Return surface coordinates and the exact surface-area Jacobian."""
    u, v = np.broadcast_arrays(np.asarray(u, dtype=float), np.asarray(v, dtype=float))
    cu, su, cv, sv = np.cos(u), np.sin(u), np.cos(v), np.sin(v)
    rho = R + r * cv
    xyz = np.stack((stretch * rho * cu, rho * su, r * sv), axis=-1)
    # For A=diag(stretch,1,1), area scales by |det A| ||A^{-T} n||.
    normal_factor = np.sqrt((cv * cu / stretch) ** 2 + (cv * su) ** 2 + sv ** 2)
    jacobian = r * rho * stretch * normal_factor
    return xyz, jacobian

def _sample_oval_torus(seed):
    """Draw N_POINTS from one oval torus, then rotate and add small noise."""
    rng = np.random.default_rng(seed)
    R = 1.0
    r = float(rng.uniform(0.62, 0.70))
    stretch = float(rng.uniform(1.30, 1.45))
    bound = float(r * (R+r) * max(1.0, stretch))
    accepted, have = [], 0
    while have < N_POINTS:
        batch = max(2048, 3 * (N_POINTS-have))
        angles = rng.uniform(0.0, 2.0*np.pi, size=(batch, 2))
        points, jacobian = _oval_surface(R, r, stretch, angles[:, 0], angles[:, 1])
        if np.any(jacobian <= 0) or np.any(jacobian > bound * (1 + 1e-13)):
            raise RuntimeError("Invalid surface-area rejection envelope")
        keep = rng.uniform(0.0, bound, size=batch) < jacobian
        accepted.append(points[keep])
        have += int(np.sum(keep))
    raw = np.concatenate(accepted, axis=0)[:N_POINTS]
    rotation, triangular = np.linalg.qr(rng.normal(size=(3, 3)))
    signs = np.sign(np.diag(triangular))
    signs[signs == 0] = 1
    rotation = rotation @ np.diag(signs)
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    raw = raw @ rotation.T
    clean_center = raw.mean(axis=0)
    clean_scale = 2.0 / float(pdist(raw).max())
    scaled_clean = (raw-clean_center) * clean_scale
    observed = scaled_clean + rng.normal(scale=0.001, size=raw.shape)
    observed_center = observed.mean(axis=0)
    observed_scale = 2.0 / float(pdist(observed).max())
    return dict(points=(observed-observed_center)*observed_scale,
                clean=(scaled_clean-observed_center)*observed_scale,
                component_ids=np.zeros(N_POINTS, dtype=np.int16))

_ORIGINAL_FAMILIES = (
    "parabola", "sphere", "fountain", "two spheres",
    "circle", "pentagon", "annulus", "3d curve",
    "two circles", "lemniscate", "concentric circles", "torus",
    "rose curve 4", "face", "double torus", "necklace",
    "rose curve 9", "olympics", "nine squares", "handlebody",
)

_ORIGINAL_SAMPLERS = (
    _original_parabola, _original_sphere, _original_fountain, _original_two_spheres,
    _original_circle, lambda n: _original_closed_curve(n, num_vertices=5),
    _original_annulus, _original_3d_curve,
    _original_two_circles, _original_lemniscate, _original_concentric_circles, _original_torus,
    lambda n: _original_rose_curve(n, a=2), _original_eyes_1,
    _original_double_torus, _original_necklace,
    _original_rose_nine_petals, _original_olympics, _original_nine_squares,
    _original_handlebody,
)


_BOTTLENECK_FAMILIES = (
    (0, "straight"), (0, "branched"),
    (1, "branched"), (2, "branched"),
    (4, "straight"), (4, "branched"),
)

# Shape counts follow the registered family definitions, not a fixed cloud count.
N_SHAPES = len(_ORIGINAL_FAMILIES) + len(_BOTTLENECK_FAMILIES) + 1
N_CLOUDS = N_SHAPES * N_SHAPE_VARIATIONS

def _generate_arrays():
    """Generate every family and variation using the one user-configured seed."""
    data = dict(
        points=np.zeros((N_CLOUDS, N_POINTS, 3), dtype=np.float64),
        labels=np.zeros(N_CLOUDS, dtype=np.int64),
        sources=np.empty(N_CLOUDS, dtype="<U11"),
        families=np.empty(N_CLOUDS, dtype="<U20"),
        shape_ids=np.empty(N_CLOUDS, dtype="<U51"),
        ambient_dimensions=np.zeros(N_CLOUDS, dtype=np.int64),
    )

    # Original formulas use np.random's legacy RNG. Preserve its previous state
    # so callers get repeatable clouds without affecting their random stream.
    previous_state = np.random.get_state()
    try:
        np.random.seed(SEED)
        for family_index, (family, sampler) in enumerate(
            zip(_ORIGINAL_FAMILIES, _ORIGINAL_SAMPLERS)
        ):
            for replicate in range(N_SHAPE_VARIATIONS):
                i = family_index * N_SHAPE_VARIATIONS + replicate
                cloud = _original_normalize(sampler(N_POINTS))
                dimension = cloud.shape[1]
                center = cloud.mean(axis=0)
                scale = 2 / float(pdist(cloud).max())
                data["points"][i, :, :dimension] = (cloud - center) * scale
                data["labels"][i] = (0, 1, 2, 4, 9)[family_index // 4]
                data["sources"][i] = "original"
                data["families"][i] = family
                data["shape_ids"][i] = f"original_{i:03d}"
                data["ambient_dimensions"][i] = dimension
    finally:
        np.random.set_state(previous_state)

    i = len(_ORIGINAL_FAMILIES) * N_SHAPE_VARIATIONS
    for family_index, (beta, theme) in enumerate(_BOTTLENECK_FAMILIES):
        for replicate in range(N_SHAPE_VARIATIONS):
            # Separate deterministic streams for each family/variation.
            child_seed = np.random.SeedSequence(
                SEED, spawn_key=(1, family_index, replicate)
            )
            sample = _sample_bottleneck(beta, theme, child_seed)
            data["points"][i] = sample["points"]
            data["labels"][i] = beta
            data["sources"][i] = "targeted"
            data["families"][i] = f"{theme}_b{beta}"
            data["shape_ids"][i] = f"targeted_{theme}_b{beta}_r{replicate}"
            data["ambient_dimensions"][i] = 2
            i += 1

    for replicate in range(N_SHAPE_VARIATIONS):
        child_seed = np.random.SeedSequence(SEED, spawn_key=(2, replicate))
        sample = _sample_oval_torus(child_seed)
        data["points"][i] = sample["points"]
        data["labels"][i] = 2
        data["sources"][i] = "toric"
        data["families"][i] = "oval_torus"
        data["shape_ids"][i] = f"torus_dev_{replicate:03d}"
        data["ambient_dimensions"][i] = 3
        i += 1

    assert i == N_CLOUDS
    return data


def build_datasets():
    """Generate native 2D/3D clouds and create the standard Turkevs objects."""
    data = _generate_arrays()
    if data["points"].shape != (N_CLOUDS, N_POINTS, 3):
        raise ValueError("Unexpected point-cloud dimensions")
    if not np.isfinite(data["points"]).all():
        raise ValueError("Generated non-finite coordinates")
    if np.any(~np.isin(data["ambient_dimensions"], [2, 3])):
        raise ValueError("Unexpected ambient dimension")

    seen = {}
    records = []
    for i in range(N_CLOUDS):
        source = str(data["sources"][i])
        family = str(data["families"][i])
        key = (source, family)
        variation = seen.get(key, 0)
        seen[key] = variation + 1
        dimension = int(data["ambient_dimensions"][i])
        if dimension == 2 and np.any(data["points"][i, :, 2] != 0):
            raise ValueError(f"Nonzero padding in 2D cloud {i}")
        points = data["points"][i, :, :dimension].copy()
        records.append((points, source, family, variation, int(data["labels"][i])))

    if len(seen) != N_SHAPES or any(count != N_SHAPE_VARIATIONS for count in seen.values()):
        raise ValueError("Unexpected number of families or variations")

    transformations = (list(TurkevsTransformation) if TRANSFORMATIONS_MODE == "all"
                       else [TurkevsTransformation.STANDARD])
    # Turkevs transformations use the legacy NumPy random stream too.
    np.random.seed(SEED)
    datasets = []
    originals = [record[0] for record in records]
    for transformation in transformations:
        if transformation is TurkevsTransformation.STANDARD:
            transformed = originals
        else:
            transformed = calculate_point_clouds_under_trnsf(
                originals, transformation=transformation.fullname
            )
        if len(transformed) != N_CLOUDS:
            raise ValueError(f"{transformation.fullname} changed the number of clouds")
        for i, (cloud, (original, source, family, variation, label)) in enumerate(
            zip(transformed, records)
        ):
            cloud = np.asarray(cloud, dtype=np.float64)
            if cloud.shape != original.shape or not np.isfinite(cloud).all():
                raise ValueError(f"Invalid transformed cloud {i}: {transformation.fullname}")
            info = TurkevsDatasetInfo(
                dataset_id="", seed=SEED, point_cloud_index=i,
                transformation=transformation, label=label,
                shape_group="turkevs" if source == "original" else "additional",
                shape_name=family, shape_variation_index=variation,
                n_point_clouds=N_CLOUDS, n_shapes=N_SHAPES,
                n_shape_variations=N_SHAPE_VARIATIONS,
            )
            datasets.append(Dataset(points=cloud, data_type="turkevs", additional_info=info))

    add_dataset_id(datasets, generate_datasets_hash(datasets))
    return datasets


def main():
    if TRANSFORMATIONS_MODE not in {"all", "standard"}:
        raise ValueError("TRANSFORMATIONS_MODE must be 'all' or 'standard'")
    setup_logging()
    logger = logging.getLogger(__name__)
    setup_name = get_setup_name(
        n_shapes=N_SHAPES, n_shape_variations=N_SHAPE_VARIATIONS,
        n_points=N_POINTS, seed=SEED,
    )
    experiment_folder = Path(TURKEVS_DATA_DIR) / setup_name
    if experiment_folder.exists():
        raise FileExistsError(f"Folder already exists: {experiment_folder}")

    logger.info("Generating %d clouds with %d points each...", N_CLOUDS, N_POINTS)
    datasets = build_datasets()
    output_path = Path(get_datasets_path(setup_name))
    if output_path.parent.resolve() != experiment_folder.resolve():
        raise ValueError("Dataset output directory does not match experiment directory")
    experiment_folder.mkdir(parents=True, exist_ok=False)
    save_to_json([d.to_dict() for d in datasets], filename=str(output_path), add_timestamp=False)
    logger.info("Saved %d datasets to %s", len(datasets), output_path)

    return experiment_folder

if __name__ == "__main__":
    main()
