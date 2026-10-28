#!/usr/bin/env python3

import os
import sys

sys.path.append(os.path.abspath('.'))

import numpy as np
from ellipsoids.common import Dataset
from ellipsoids.common import EllipsoidParameters, Parameters
from ellipsoids.common import Experiment
from ellipsoids.data_handling import sample_from_circle
from ellipsoids.visualisation import plot_experiment


import pandas as pd
import time
import matplotlib.pyplot as plt


OUTPUT_PATH = "data/benchmarking"

def benchmark(func, points, repeats=5):
    times = []

    for _ in range(repeats):
        start = time.perf_counter()
        func(points)
        end = time.perf_counter()

        times.append(end - start)

    return np.median(times)


def run_experiment(dataset, parameters):
    experiment = Experiment(dataset, parameters)
    experiment.run()

def main():

   
    # Point cloud sizes to test
    sizes = [50, 100, 200, 400, 600, 800, 1000, 1200]

    results = []

    rng = np.random.default_rng(0)

    for n in sizes:

        params_ellipsoid = EllipsoidParameters(
            nbhd_size=np.floor(n/10),
            save_simplex_tree=False,
            save_ellipsoid_list=False,
            axes_ratios=np.asarray([2,1]),
            collapse_edges=True,
            )
        params_ball = Parameters(
            collapse_edges=True,
            ) 


        print(f"benchmarking n={n}")
        points = rng.random((n, 2))

        dataset = Dataset(points, "random")

        for name, func in [
            ("BALL-RIPS", lambda dataset: run_experiment(dataset, params_ball) ),
            ("ELLIPSOIDS-RIPS", lambda dataset: run_experiment(dataset, params_ellipsoid)),
        ]:
            runtime = benchmark(func, dataset, 10)

            results.append({
                "algorithm": name,
                "n_points": n,
                "runtime": runtime,
            })

    df = pd.DataFrame(results)
    print(df)

    output_name = os.path.join(OUTPUT_PATH, "timing_benchmark_results")
    csv_name = f"{output_name}.csv"
    png_name = f"{output_name}.png"
    df.to_csv(csv_name, index=False)

    plot_benchmark(csv_name, png_name)




def plot_benchmark(
    csv_path="benchmark_results.csv",
    output_path="benchmark_plot.png",
    log_scale=False,
):
    df = pd.read_csv(csv_path)

    # Use column positions, not column names
    algorithm_col, size_col, runtime_col = df.columns[:3]

    df = (
        df.groupby([algorithm_col, size_col], as_index=False)
        [runtime_col].median()
    )

    algorithms = list(df[algorithm_col].unique())

    if len(algorithms) != 2:
        raise ValueError("Expected exactly two algorithms.")

    fig, ax1 = plt.subplots(figsize=(8, 5))

    # Left axis: runtimes
    for name, group in df.groupby(algorithm_col, sort=False):
        group = group.sort_values(size_col)

        ax1.plot(
            group[size_col],
            group[runtime_col],
            marker="o",
            linewidth=2,
            label=name,
        )

    ax1.set_xlabel("Number of points")
    ax1.set_ylabel("Runtime (seconds)")
    ax1.grid(True, alpha=0.3)

    # Right axis: runtime ratio
    pivot = df.pivot(
        index=size_col,
        columns=algorithm_col,
        values=runtime_col,
    )

    ratio = (
        pivot[algorithms[0]] / pivot[algorithms[1]]
    ).dropna()

    ax2 = ax1.twinx()

    ax2.plot(
        ratio.index,
        ratio.values,
        marker="s",
        linestyle="--",
        color="purple",
        linewidth=2,
        label=f"{algorithms[0]} / {algorithms[1]}",
    )

    ax2.axhline(1, color="purple", linestyle=":", alpha=0.5)
    ax2.set_ylabel("Runtime ratio", color="purple")
    ax2.tick_params(axis="y", colors="purple")
    ax2.spines["right"].set_color("purple")

    if log_scale:
        ax1.set_xscale("log")
        ax1.set_yscale("log")

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()

    ax1.legend(lines1 + lines2, labels1 + labels2)

    ax1.set_title("Persistent Homology Runtime Comparison")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)    

    

if __name__ == '__main__':
    main()
