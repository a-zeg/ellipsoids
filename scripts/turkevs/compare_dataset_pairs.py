#!/usr/bin/env python3

import json
import numpy as np

FILE_1 = "data/turkevs/turkevs_ns=27_nsv=5_np=600_s=0/datasets_turkevs_ns=27_nsv=5_np=600_s=0.json"
FILE_2 = "data/turkevs/turkevs_ns=27_nsv=5_np=600_s=0_/datasets_turkevs_ns=27_nsv=5_np=600_s=0.json"


with open(FILE_1, "r") as f:
    data1 = json.load(f)

with open(FILE_2, "r") as f:
    data2 = json.load(f)


print("Number of datasets:")
print(len(data1), len(data2))

if len(data1) != len(data2):
    raise ValueError("Different number of datasets.")


for i, (d1, d2) in enumerate(zip(data1, data2)):

    # Compare point clouds
    points1 = np.asarray(d1["points"])
    points2 = np.asarray(d2["points"])

    if not np.array_equal(points1, points2):
        raise ValueError(
            f"Points differ at dataset index {i}"
        )

    # Compare metadata
    if d1["additional_info"] != d2["additional_info"]:
        raise ValueError(
            f"Additional info differs at dataset index {i}:\n"
            f"{d1['additional_info']}\n"
            f"{d2['additional_info']}"
        )

    # Compare anything else
    keys_to_ignore = {"points", "additional_info"}

    rest1 = {
        k: v
        for k, v in d1.items()
        if k not in keys_to_ignore
    }

    rest2 = {
        k: v
        for k, v in d2.items()
        if k not in keys_to_ignore
    }

    if rest1 != rest2:
        raise ValueError(
            f"Other fields differ at dataset index {i}"
        )


print("✓ Files contain the same dataset values.")
