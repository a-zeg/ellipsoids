import os
import sys

sys.path.append(os.path.abspath('.'))

import numpy as np


from ellipsoids.data_handling import sample_from_figure_eight



NUM_POINTS = 300




def sample_spiral(n=NUM_POINTS, turns=3, inner_radius=0.1, outer_radius=1.0):
    """
    Open 1D spiral in R^2.
    No holes.
    """
    t = np.random.uniform(0, 1, size=n)

    theta = 2 * np.pi * turns * t
    r = inner_radius + (outer_radius - inner_radius) * t

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    X = np.array([x, y]).T

    return X


def sample_3d_bottleneck(n=NUM_POINTS, gap=0.12):
    """
    Closed non-planar 1D curve in R^3.

    In xy projection it resembles a figure eight, but the two
    branches near the crossing are separated in z.

    Smaller gap -> stronger bottleneck / near-self-intersection.

    Topologically this is a circle (1 hole).
    """
    t = np.random.uniform(0, 2 * np.pi, size=n)

    x = np.sin(t)
    y = np.sin(2 * t)
    z = gap * np.cos(t)

    X = np.array([x, y, z]).T

    return X

def _sample_curve_uniform_arclength(curve, n, resolution=5000):
    """
    Randomly sample a parametric curve uniformly with respect
    to arc length.

    curve(t) takes t in [0,1] and returns an (N,d) array.
    """

    t_dense = np.linspace(0, 1, resolution)
    X_dense = curve(t_dense)

    ds = np.linalg.norm(
        np.diff(X_dense, axis=0),
        axis=1
    )

    s_dense = np.concatenate([
        [0],
        np.cumsum(ds)
    ])

    total_length = s_dense[-1]

    # Uniform random positions along arc length
    s = np.random.uniform(
        0,
        total_length,
        size=n
    )

    t = np.interp(
        s,
        s_dense,
        t_dense
    )

    return curve(t)


def _sample_curve_family(curves, n, resolution=5000):
    """
    Uniform random sampling over a collection of curves.

    Curves are selected proportional to their arc length,
    so the resulting sampling density is uniform over the
    entire 1D shape.
    """

    lengths = []

    for curve in curves:

        t = np.linspace(
            0,
            1,
            resolution
        )

        X = curve(t)

        length = np.sum(
            np.linalg.norm(
                np.diff(X, axis=0),
                axis=1
            )
        )

        lengths.append(length)

    lengths = np.asarray(lengths)

    probabilities = (
        lengths / lengths.sum()
    )

    # Pick which branch each point belongs to
    branch_ids = np.random.choice(
        len(curves),
        size=n,
        p=probabilities
    )

    dimension = curves[0](
        np.array([0.5])
    ).shape[1]

    X = np.empty((n, dimension))

    for j, curve in enumerate(curves):

        indices = np.where(
            branch_ids == j
        )[0]

        if len(indices) == 0:
            continue

        X[indices] = (
            _sample_curve_uniform_arclength(
                curve,
                len(indices),
                resolution=resolution
            )
        )

    return X


# ============================================================
# 2 holes
# ============================================================

def sample_2_holes_bottleneck(
    n=NUM_POINTS,
    gap=0.20
):
    """
    1D theta graph in R^2.

    Three distinct paths connect the same two endpoints,
    giving beta_1 = 2.

    Decreasing gap makes the branches increasingly close,
    producing a harder problem for Rips.
    """

    amplitudes = [
        -gap,
        0.0,
        gap
    ]

    curves = []

    for amplitude in amplitudes:

        def curve(t, a=amplitude):

            x = 2 * t - 1
            y = a * np.sin(np.pi * t)

            return np.array([
                x,
                y
            ]).T

        curves.append(curve)

    return _sample_curve_family(
        curves,
        n
    )


# ============================================================
# 4 holes
# ============================================================

def sample_4_holes_bottleneck(
    n=NUM_POINTS,
    gap=0.15
):
    """
    Five nearby curved paths joining the same two endpoints.

    Number of independent cycles:
        5 paths - 1 = 4 holes.

    All branches are long, so sampling is much more balanced
    than with a thin rectangular grid.
    """

    amplitudes = gap * np.array([
        -2,
        -1,
         0,
         1,
         2
    ])

    curves = []

    for amplitude in amplitudes:

        def curve(t, a=amplitude):

            x = 2 * t - 1
            y = a * np.sin(np.pi * t)

            # Slight asymmetry makes it less grid-like
            y += (
                0.08
                * a
                * np.sin(2 * np.pi * t)
            )

            return np.array([
                x,
                y
            ]).T

        curves.append(curve)

    return _sample_curve_family(
        curves,
        n
    )


# ============================================================
# 9 holes
# ============================================================

def sample_9_holes_bottleneck(
    n=NUM_POINTS,
    radius=0.32,
    twist=0.7
):
    """
    3D cage / multi-theta graph.

    Ten curved strands connect the same two endpoints.
    Therefore

        beta_1 = 10 - 1 = 9.

    The strands are arranged around a narrow 3D tube.
    Small radius produces many Euclidean-close but
    topologically distinct branches.
    """

    number_of_paths = 10

    curves = []

    for j in range(number_of_paths):

        phi = (
            2 * np.pi
            * j
            / number_of_paths
        )

        def curve(t, phi=phi):

            x = 2 * t - 1

            # Radius vanishes at both endpoints,
            # so every branch meets there.
            r = (
                radius
                * np.sin(np.pi * t)
            )

            # A small twist makes the object genuinely 3D
            angle = (
                phi
                + twist
                * np.sin(np.pi * t)
            )

            y = (
                r
                * np.cos(angle)
            )

            z = (
                r
                * np.sin(angle)
            )

            return np.array([
                x,
                y,
                z
            ]).T

        curves.append(curve)

    return _sample_curve_family(
        curves,
        n
    )



















    


if __name__ == "__main__":

    # X0 = sample_spiral()
    # X1 = sample_3d_bottleneck(gap=0.08)

    # X2 = sample_2_holes_bottleneck(gap=0.12)
    # X4 = sample_4_holes_bottleneck(gap=0.12)
    # X9 = sample_9_holes_bottleneck(gap=0.12)


    import plotly.graph_objects as go
    from plotly.subplots import make_subplots


# Generate point clouds
    X0 = sample_spiral()
    X1 = sample_3d_bottleneck(gap=0.08)

    X2 = sample_2_holes_bottleneck(gap=0.12)
    X4 = sample_4_holes_bottleneck(gap=0.12)
    X9 = sample_9_holes_bottleneck(radius=0.22)


# ============================================================
        # Plot
        # ============================================================

    fig = make_subplots(
        rows=2,
        cols=3,
        specs=[
            [{"type": "xy"}, {"type": "scene"}, {"type": "xy"}],
            [{"type": "xy"}, {"type": "xy"}, None],
        ],
        subplot_titles=[
            "Spiral — 0 holes",
            "3D bottleneck — 1 hole",
            "2 holes",
            "4 holes",
            "9 holes",
        ],
        horizontal_spacing=0.06,
        vertical_spacing=0.12,
        )


    # ------------------------------------------------------------
        # Spiral
        # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter(
            x=X0[:, 0],
            y=X0[:, 1],
            mode="markers",
            marker=dict(size=4),
            name="Spiral",
        ),
        row=1,
        col=1,
        )


    # ------------------------------------------------------------
        # 3D bottleneck
        # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter3d(
            x=X1[:, 0],
            y=X1[:, 1],
            z=X1[:, 2],
            mode="markers",
            marker=dict(size=3),
            name="3D bottleneck",
        ),
        row=1,
        col=2,
        )


    # ------------------------------------------------------------
        # 2 holes
        # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter(
            x=X2[:, 0],
            y=X2[:, 1],
            mode="markers",
            marker=dict(size=4),
            name="2 holes",
        ),
        row=1,
        col=3,
        )


    # ------------------------------------------------------------
        # 4 holes
        # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter(
            x=X4[:, 0],
            y=X4[:, 1],
            mode="markers",
            marker=dict(size=4),
            name="4 holes",
        ),
        row=2,
        col=1,
        )


    # ------------------------------------------------------------
        # 9 holes
        # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter(
            x=X9[:, 0],
            y=X9[:, 1],
            mode="markers",
            marker=dict(size=4),
            name="9 holes",
        ),
        row=2,
        col=2,
        )


    # ============================================================
        # Formatting
        # ============================================================

    fig.update_layout(
        width=1100,
        height=750,
        showlegend=False,
        title="1D Point Cloud Test Shapes",
        )


    # Hide 2D axes
    for row, col in [
        (1, 1),
        (1, 3),
        (2, 1),
        (2, 2),
        ]:
        fig.update_xaxes(
            visible=False,
            row=row,
            col=col,
        )

        fig.update_yaxes(
            visible=False,
            row=row,
            col=col,
        )


    # 3D appearance
    fig.update_scenes(
        xaxis_visible=False,
        yaxis_visible=False,
        zaxis_visible=False,
        aspectmode="data",
        )


    fig.show()
