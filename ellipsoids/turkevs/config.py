import os

# folder structure
TURKEVS_DATA_DIR = os.path.join("data", "turkevs")
REL_DATASETS_DIR = ""
REL_EXPERIMENT_SUMMARIES_DIR = ""
REL_CLASSIFICATION_RESULTS_DIR = "classification_results"
REL_AGGREGATED_RESULTS_DIR = "aggregated_results"

# name templating
def get_setup_name(n_point_clouds: int, n_points: int, seed: int):
    return f"turkevs_npc={n_point_clouds}_np={n_points}_s={seed}"


def resolve_setup_path(setup_input: str):
    if os.path.isabs(setup_input) or os.path.sep in setup_input:
        setup_path = os.path.normpath(setup_input)
    else:
        setup_path = os.path.join(TURKEVS_DATA_DIR, setup_input)
    return setup_path


# path builders
def get_path(setup_input: str, rel_subdir: str, prefix: str = "", suffix: str = "") -> str:
    """
    General path builder for dataset-related files.

    Parameters:
    - setup_input: either a setup name (e.g., 'turkevs_npc=...') or full path
    - rel_subdir: relative subfolder name (e.g., 'datasets', 'classification_results')
    - prefix: file prefix (e.g. 'datasets_', 'summaries_')
    - suffix: file suffix (e.g. '.json', "part=1-20")

    Returns:
    - full path to the file
    """
    setup_path = resolve_setup_path(setup_input)
    setup_name = os.path.basename(setup_path)
    return os.path.join(setup_path, rel_subdir, f"{prefix}{setup_name}{suffix}")



def get_datasets_path(setup_input: str) -> str:
    return get_path(setup_input, REL_DATASETS_DIR, "datasets_", ".json")

def get_experiment_summaries_path(setup_input: str, suffix="") -> str:
    return get_path(setup_input, REL_EXPERIMENT_SUMMARIES_DIR, "summaries_", f"{suffix}.jsonl.gz")

def get_classification_results_path(setup_input: str) -> str:
    return get_path(setup_input, REL_CLASSIFICATION_RESULTS_DIR, "classification_")

def get_aggregated_results_basepath(setup_input: str) -> str:
    return get_path(setup_input, REL_AGGREGATED_RESULTS_DIR, "aggregated_results_")

def get_classification_results_folder(setup_input: str) -> str:
    if os.path.isabs(setup_input) or os.path.sep in setup_input:
        setup_path = os.path.normpath(setup_input)
    else:
        setup_path = os.path.join(TURKEVS_DATA_DIR, setup_input)

    return os.path.join(setup_path, REL_CLASSIFICATION_RESULTS_DIR)
