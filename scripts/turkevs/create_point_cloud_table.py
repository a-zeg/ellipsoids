#!/usr/bin/env python3
"""Render point clouds from 0-generate_point_clouds.py and write a LaTeX table.

Default behavior: save exactly one representative image per shape family, so the
output contains 27 images total (one per family) rather than all 135 replicate
variations. The LaTeX table uses these same representative images.

If you ever want every variation again, pass --all-variations.

Run: python create_point_cloud_table.py --source 0-generate_point_clouds.py
Requires: numpy, scipy, matplotlib. No ellipsoids installation is needed.
"""

import argparse
import ast
import csv
from collections import Counter, defaultdict
from pathlib import Path
import math
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist

CONSTANT_NAMES = {
    "N_SHAPE_VARIATIONS", "N_POINTS", "SEED", "_TAU",
    "_ORIGINAL_FAMILIES", "_ORIGINAL_SAMPLERS", "_BOTTLENECK_FAMILIES",
    "N_SHAPES", "N_CLOUDS",
}
ROW_END = "\\\\"


def load_clouds(source: Path):
    """Execute only the source's data-generation definitions, not its I/O/main."""
    parsed = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    selected = []
    for node in parsed.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("_") and node.name != "__main__":
                selected.append(node)
        elif isinstance(node, ast.Assign):
            names = {target.id for target in node.targets if isinstance(target, ast.Name)}
            if names & CONSTANT_NAMES:
                selected.append(node)
    namespace = {"np": np, "math": math, "pdist": pdist}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(source), "exec"), namespace)
    required = ("_generate_arrays", "_ORIGINAL_FAMILIES", "_BOTTLENECK_FAMILIES")
    if not all(k in namespace for k in required):
        raise ValueError("The supplied file lacks expected generator definitions")
    return namespace["_generate_arrays"](), namespace


def safe_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def compute_shared_limits(data, limit_padding: float):
    """Compute one common symmetric limit per ambient dimension."""
    maxima = {2: 0.0, 3: 0.0}
    for i in range(len(data["labels"])):
        dim = int(data["ambient_dimensions"][i])
        points = np.asarray(data["points"][i, :, :dim], dtype=float)
        maxima[dim] = max(maxima[dim], float(np.max(np.abs(points))))
    return {dim: maxima[dim] * (1.0 + limit_padding) for dim in maxima}


def render_cloud(points, dimension: int, destination: Path, dpi: int, limits):
    """Use identical square canvases, marker styles and shared spatial limits."""
    fig = plt.figure(figsize=(1.42, 1.42), dpi=dpi, facecolor="white")
    limit = float(limits[dimension])
    ticks = np.linspace(-1, 1, 5)
    if dimension == 2:
        ax = fig.add_axes((0.13, 0.13, 0.75, 0.75))
        ax.scatter(points[:, 0], points[:, 1], s=0.9, c="black",
                   alpha=0.82, linewidths=0, rasterized=True)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.grid(True, color="#d8d8d8", linewidth=0.35)
        ax.tick_params(length=1.6, width=0.35, labelbottom=False, labelleft=False)
        for spine in ax.spines.values():
            spine.set_linewidth(0.45)
            spine.set_color("#777777")
    elif dimension == 3:
        ax = fig.add_axes((0.03, 0.03, 0.94, 0.94), projection="3d")
        ax.scatter(points[:, 0], points[:, 1], points[:, 2],
                   s=0.9, c="black", alpha=0.8, depthshade=False,
                   linewidths=0, rasterized=True)
        for set_lim in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
            set_lim(-limit, limit)
        ax.set_box_aspect((1, 1, 1), zoom=0.90)
        ax.view_init(elev=23, azim=-55)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_ticks(ticks)
            axis.set_ticklabels([])
            axis._axinfo["grid"]["color"] = (0.80, 0.80, 0.80, 1)
            axis._axinfo["grid"]["linewidth"] = 0.4
            axis._axinfo["axisline"]["linewidth"] = 0.4
        ax.tick_params(length=0, pad=-3)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
    else:
        raise ValueError(f"Unsupported ambient dimension {dimension}")
    fig.savefig(destination, dpi=dpi, facecolor="white")
    plt.close(fig)


def write_table(groups, output: Path, image_width: float, column_sep_cm: float):
    """Write a simple two-column table with direct includegraphics commands.

    The second column contains the example images directly, with no nested
    tabular and no rendered horizontal space between images.
    """
    _ = column_sep_cm  # kept for CLI compatibility; intentionally unused
    tex = [
        "% Generated by create_point_cloud_table.py",
        "% Requires \\usepackage{graphicx,booktabs,array}",
        "% Compile from this directory, or adjust the figures/ path.",
        "\\begin{tabular}{@{}c@{\\hspace{0.2em}}c@{}}",
        "\\toprule",
        f"No. holes & Example point clouds {ROW_END}",
        "\\midrule",
    ]
    for label in sorted(groups):
        imgs = groups[label]
        commands = ''.join(
            f"\\includegraphics[width={image_width}cm]{{figures/{path.stem}}}"
            for path in imgs
        )
        tex.append(f"\\raisebox{{0.85em}}{{{label}}} & {commands} {ROW_END}")
    tex += ["\\bottomrule", "\\end{tabular}", ""]
    output.write_text("\n".join(tex), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path(__file__).with_name("0-generate_point_clouds.py"))
    parser.add_argument("--output", type=Path, default=Path("point_cloud_table"))
    parser.add_argument("--image-width", type=float, default=0.96,
                        help="Width of each LaTeX image in cm (default 0.96)")
    parser.add_argument("--column-sep", type=float, default=0.01,
                        help="Horizontal spacing between example images in cm (default 0.04)")
    parser.add_argument("--limit-padding", type=float, default=0.06,
                        help="Fractional safety margin added to the shared plot limits (default 0.06)")
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--all-variations", action="store_true",
                        help="Save all replicate variations (135 images) instead of one representative per family")
    args = parser.parse_args()

    if args.image_width <= 0 or args.column_sep < 0 or args.limit_padding < 0 or args.dpi < 50:
        parser.error("image width and dpi must be positive; column-sep and limit-padding must be nonnegative")

    data, config = load_clouds(args.source)
    limits = compute_shared_limits(data, args.limit_padding)
    images = args.output / "figures"
    images.mkdir(parents=True, exist_ok=True)

    n = len(data["labels"])
    groups = defaultdict(list)
    metadata = []
    families_seen = set()
    rendered = 0

    for i in range(n):
        label = int(data["labels"][i])
        dim = int(data["ambient_dimensions"][i])
        family = str(data["families"][i])
        source = str(data["sources"][i])
        key = (source, family)
        variation = len([m for m in metadata if (m["source"], m["family"]) == key])
        points = np.asarray(data["points"][i, :, :dim])

        if points.shape != (int(config["N_POINTS"]), dim) or dim not in (2, 3):
            raise ValueError(f"Unexpected shape/dimension at index {i}")
        if not np.isfinite(points).all():
            raise ValueError(f"Nonfinite coordinates at index {i}")
        if dim == 2 and np.any(data["points"][i, :, 2] != 0):
            raise ValueError(f"2D cloud {i} has nonzero third coordinates")

        # Default: save exactly one representative image per family.
        if (not args.all_variations) and key in families_seen:
            continue

        suffix = f"_v{variation:02d}" if args.all_variations else ""
        filename = f"holes_{label}_{safe_name(source)}_{safe_name(family)}{suffix}.png"
        output_path = images / filename
        render_cloud(points, dim, output_path, args.dpi, limits)

        metadata.append(dict(index=i, holes=label, dimension=dim, source=source,
                             family=family, variation=variation, image=str(output_path)))
        rendered += 1

        if key not in families_seen:
            groups[label].append(output_path)
            families_seen.add(key)

    if len(families_seen) != int(config["N_SHAPES"]):
        raise AssertionError("Missing/duplicate families")

    counts = Counter((m["source"], m["family"]) for m in metadata)
    if args.all_variations:
        expected = int(config["N_SHAPE_VARIATIONS"])
        if set(counts.values()) != {expected}:
            raise AssertionError("Unexpected variations per family")
    else:
        if len(counts) != int(config["N_SHAPES"]) or set(counts.values()) != {1}:
            raise AssertionError("Expected exactly one saved image per family")

    with (args.output / "point_cloud_index.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=metadata[0].keys())
        writer.writeheader()
        writer.writerows(metadata)

    write_table(groups, args.output / "point_cloud_table.tex", args.image_width, args.column_sep)

    print(f"Saved {rendered} PNG plots from {len(families_seen)} families in {images}")
    print("Family counts by hole class:", dict(sorted((k, len(v)) for k, v in groups.items())))
    print("Shared plot limits:", limits)
    print(f"Table: {args.output / 'point_cloud_table.tex'}")
    print(f"Index: {args.output / 'point_cloud_index.csv'}")


if __name__ == "__main__":
    main()
