#!/usr/bin/env python3

import os
import sys
import gzip
import shutil
import glob
import argparse
from typing import cast
from collections import defaultdict

sys.path.append(os.path.abspath('.'))

from ellipsoids.common import TurkevsDatasetInfo
from ellipsoids.turkevs.turkevs_utils import flush_buffer, read_experiment_summaries_from_jsonl
from ellipsoids.turkevs.config import REL_SUMMARIES_DIR


def get_common_param_hashes(summaries):
    """Return the set of parameter hashes that appear for all point clouds."""
    from collections import defaultdict

    point_cloud_to_param_hashes = defaultdict(set)
    for summary in summaries:
        info = cast(TurkevsDatasetInfo, summary.dataset_summary.additional_info)
        pc_id = (info.transformation.shortname, info.point_cloud_index)
        param_hash = hash(summary.parameters)
        point_cloud_to_param_hashes[pc_id].add(param_hash)

    all_hash_sets = point_cloud_to_param_hashes.values()
    if not all_hash_sets:
        return set()
    return set.intersection(*all_hash_sets)



def write_filtered_summaries(summaries, allowed_hashes, output_path):
    filtered_summaries = [summary.to_dict() for summary in summaries if hash(summary.parameters) in allowed_hashes]
    flush_buffer(filtered_summaries, output_path, compress=True)
    # with gzip.open(output_path, "wt") as f_out:
    #     for summary in summaries:
    #         if hash(summary.parameters) in allowed_hashes:
    #             f_out.write(summary.serialize_json() + "\n")
    print(f"Wrote filtered summaries to {output_path}")



def merge_jsonl_gz(folder):
    summaries_folder = os.path.join(folder, REL_SUMMARIES_DIR)
    partial_summaries = [os.path.join(summaries_folder,path)
                         for path in os.listdir(summaries_folder)
                         if "summaries" in path and "part" in path and "merged" not in path]
    output_path = os.path.join(summaries_folder, f"summaries_{os.path.basename(folder)}_merged.jsonl.gz")

    with gzip.open(output_path, 'wb') as f_out:
        for part in partial_summaries:
            with gzip.open(part, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    print(f"Merged {len(partial_summaries)} parts into: {output_path}")
    return output_path



def validate_summaries(summaries):
    if not summaries:
        raise ValueError("No summaries to load.")

    dataset_ids = set()
    seeds = set()
    n_point_clouds_set = set()


    point_cloud_to_param_hashes = defaultdict(set)
    param_hash_set = set()
    point_cloud_metadata = {}

    for summary in summaries:

        # checking labels and transformations consistency
        dataset_info = cast(TurkevsDatasetInfo, summary.dataset_summary.additional_info)
        seeds.add(dataset_info.seed)
        dataset_ids.add(dataset_info.dataset_id)
        n_point_clouds_set.add(dataset_info.n_point_clouds)
        pc_identifier = dataset_info.point_cloud_index
        transformation = dataset_info.transformation
        pc_identifier = (transformation, pc_identifier)
        label = dataset_info.label
        if pc_identifier in point_cloud_metadata:
            if point_cloud_metadata[pc_identifier] != label:
                # raise ValueError(f"Inconsistent metadata for point cloud index {point_cloud_index}: {point_cloud_metadata[point_cloud_index]} vs {key}")
                print(f"Inconsistent metadata for point cloud index {pc_identifier}: {point_cloud_metadata[pc_identifier]} vs {label}")
        else:
            point_cloud_metadata[pc_identifier] = label

        # checking if all parameter sets were calculated for each point cloud index
        param_hash = hash(summary.parameters)
        point_cloud_to_param_hashes[pc_identifier].add(param_hash)
        param_hash_set.add(param_hash)

    def check_single_val(name, s):
        if len(s) != 1:
            raise ValueError(f"Inconsistent {name}: {s}")
        return next(iter(s))

    dataset_id = check_single_val("dataset_id", dataset_ids)
    seed = check_single_val("seed", seeds)
    n_point_clouds = check_single_val("n_point_clouds", n_point_clouds_set)

    expected_param_count = len(param_hash_set)
    incomplete = []

    for pc_identifier, hashes in point_cloud_to_param_hashes.items():
        if len(hashes) != expected_param_count:
            incomplete.append((pc_identifier, len(hashes), expected_param_count))

    # Reporting
    print(f"\nDataset ID: {dataset_id}")
    print(f"Seed: {seed}")
    print(f"Total point clouds: {len(point_cloud_to_param_hashes)}")
    print(f"Parameters per point cloud: {expected_param_count}")

    if incomplete:
        print(f"\n {len(incomplete)} point clouds are missing parameter results:")
        for idx, count, expected in incomplete:
            print(f"  - Index {idx}: {count}/{expected} parameter sets")
        raise RuntimeError("Validation failed: incomplete parameter coverage.")
    else:
        print("\nAll point clouds have full parameter coverage and consistent metadata.")




def extract_last_n_entries_from_jsonl_gz(input_path: str, output_path: str, n: int):
    """Extract the last n lines from a .jsonl.gz file and write them to a new .jsonl.gz file."""
    # Read all lines (decompressed)
    with gzip.open(input_path, "rt") as f_in:
        lines = f_in.readlines()

    if n > len(lines):
        print(f"Requested {n} lines, but file only contains {len(lines)} lines. Copying all.")
        n = len(lines)

    last_lines = lines[-n:]

    # Write selected lines to new file (compressed)
    with gzip.open(output_path, "wt") as f_out:
        for line in last_lines:
            f_out.write(line)

    print(f"Wrote last {n} entries to {output_path}")



def count_lines_in_jsonl_gz(path: str) -> int:
    """Returns the number of lines (entries) in a .jsonl.gz file."""
    with gzip.open(path, "rt") as f:
        return sum(1 for _ in f)




if __name__ == "__main__":
    folder = "data/turkevs/turkevs_npc=100_np=300_s=0"
    merged_filename = merge_jsonl_gz(folder)
    summaries = read_experiment_summaries_from_jsonl(merged_filename)
    validate_summaries(summaries)


    # common_params = get_common_param_hashes(summaries)
    # merged_filtered_basename = f"summaries_{os.path.basename(folder)}_merged_filtered.jsonl.gz"
    # merged_filtered_filename = os.path.join(folder, merged_filtered_basename)
    # # write_filtered_summaries(summaries, common_params, merged_filtered_filename)
    # mf_summaries = read_experiment_summaries_from_jsonl(merged_filtered_filename)
    # validate_summaries(mf_summaries)


    # filename = "data/turkevs/turkevs_npc=100_np=300_s=0/partial/summaries_turkevs_npc=100_np=300_s=0.jsonl.gz"
    # filename_out = "data/turkevs/turkevs_npc=100_np=300_s=0/partial/summaries_turkevs_npc=100_np=300_s=0_part=1-70.jsonl.gz"

    # filename_count = "data/turkevs/turkevs_npc=100_np=300_s=0/partial/summaries_turkevs_npc=100_np=300_s=0_part=70-140.jsonl.gz"
    # lines_n = count_lines_in_jsonl_gz(filename_count)

    # extract_last_n_entries_from_jsonl_gz(filename, filename_out, lines_n)
