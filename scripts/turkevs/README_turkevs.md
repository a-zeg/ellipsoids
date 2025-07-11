# Introduction

This README contains instructions for running experiments analogous to [[1]](#1), including classification based on ellipsoid and Rips complexes.

These scripts make heavy use of the code published [here](https://github.com/renata-turkes/turkevs2022on).

# Organisation

The following scripts are available. If you're running these experiments from scratch, the number in the filename indicates the order in which the scripts should be executed.

Core:
- `0-generate_turkevs_datasets.py` generates datasets according to the specified number of point clouds, number of points and the seed and saves them in a `.json` file.
- `1-calculate_ellipsoid_barcodes.py` calculates barcodes corresponding to the ellipsoid complexes with the specified parameters and saves them to a `.jsonl.gz` file. Each line of this file corresponds to a single experiment summary and, depending on the number of point clouds and parameters, there might be thousands of lines.
- `2-run_classification.py` reads in the `.jsonl.gz` containing the barcodes and runs (multiple) classification using those, as well as the barcodes and other signatures used in [[1]](#1).
- `3-aggregate_classification_results.py` calculates the mean classification accuracy per parameter set and creates plots.

Optional:
- `1.5-merge_results.py` merges summaries corresponding to different dataset chunks into a single summaries file. This is useful after running parallel jobs on a cluster.
- `4-compare_aggreated_results.py` is for checking which parameters lead to the best classification performance across all datasets. 
- `submit_euler_jobs.py` for submitting SLURM jobs on the ETH HPC.
- `submit_euler_jobs_by_dataset_chunks.py`for submitting SLURM jobs on the ETH HPC.


# References

<a id="1">[1]</a> 
R.Turkeš, G.F.Montúfar, and N.Otter. 
"On the Effectiveness of Persistent Homology”. 
In: Advances in Neural Information Processing Systems. Ed. by S. Koyejo et al. Vol. 35. Curran Associates, Inc., 2022, pp. 35432–35448.
