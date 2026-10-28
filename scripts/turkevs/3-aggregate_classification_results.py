import os
import sys
import logging
from statistics import mean

sys.path.append(os.path.abspath('.'))

import numpy as np
from ellipsoids.data_handling import save_to_json
from ellipsoids.data_handling import parse_args
from ellipsoids.turkevs.turkevs_utils import add_aggregate_results_labels, plot_aggregated_results
from ellipsoids.turkevs.turkevs_utils import calculate_aggregated_classification_results
from ellipsoids.turkevs.config import get_aggregated_results_basepath, get_classification_results_folder
from ellipsoids.logging_setup import setup_logging
from ellipsoids.turkevs.turkevs_utils import make_parameter_filter
from ellipsoids.common import ComplexType, ComplexSubtype, EllipsoidParameters


setup_logging()
logger = logging.getLogger(__name__)


def generate_latex_table(results_list):

    header = r"""\begin{tabular}{|l|l|l|c|c|}
\hline
\textbf{Results label} & \textbf{Complex type} & \textbf{Subtype} & \textbf{Neighbourhood size} & \textbf{Axes ratios} \\
\hline"""

    rows = []
    for result in results_list:
        label = result.results_label
        if result.parameters is None:
            complex_type = complex_subtype = nbhd_size = axes_ratios = ""
        else:
            complex_type = result.parameters.complex_type
            complex_subtype = result.parameters.complex_subtype
            if isinstance(result.parameters, EllipsoidParameters):
                nbhd_size = result.parameters.nbhd_size
                axes_ratios = result.parameters.axes_ratios if result.parameters.axes_ratios is not None else "PCA"
            else:
                nbhd_size = ""
                axes_ratios = ""

        row = f"{label} & {complex_type} & {complex_subtype} & {nbhd_size} & {axes_ratios} \\\\ \\hline"
        rows.append(row)

    table = "\n".join([header] + rows + [r"\end{tabular}"])
    return table



def main():

    setup_path = "data/turkevs/turkevs_ns=27_nsv=5_np=600_s=0"

    filter_fn = None
    filters_snipe = None

    # filter_fn = make_parameter_filter(
    #     complex_type=[
    #         ComplexType.ELLIPSOID
    #         ],
    #     # axes_ratios=None,
    #     r_spherisize=[np.inf]
    # )
    # filters_snipe = [
    #     # make_parameter_filter(
    #     #     complex_type = ComplexType.ELLIPSOID,
    #     #     complex_subtype = ComplexSubtype.RIPS,
    #     #     nbhd_size = 17,
    #     #     axes_ratios = np.array([2,1]),
    #     #     r_spherisize = np.inf),
    #     make_parameter_filter(
    #         complex_type = ComplexType.ELLIPSOID,
    #         complex_subtype = ComplexSubtype.RIPS,
    #         nbhd_size = 25,
    #         axes_ratios = None,
    #         r_spherisize = np.inf),
    #     # make_parameter_filter(
    #     #     complex_type = ComplexType.ELLIPSOID,
    #     #     complex_subtype = ComplexSubtype.RIPS,
    #     #     nbhd_size = 18,
    #     #     axes_ratios = None,
    #     #     r_spherisize = np.inf),
    #     # make_parameter_filter(
    #     #     complex_type = ComplexType.BALL,
    #     #     complex_subtype = ComplexSubtype.RIPS),
    # ]


    args = parse_args()
    setup_path = args.folder or setup_path
    logger.info(f"Calculating aggregated results for the setup {setup_path}...")

    classification_results_folder = get_classification_results_folder(setup_path)
    aggregated_results = calculate_aggregated_classification_results(classification_results_folder)

    filtered_results = []

    if filter_fn is None or filters_snipe is None:
        filtered_results = aggregated_results
    else:
        for results in aggregated_results:
            if results.parameters is None: filtered_results.append(results)
            if (filter_fn is not None and filter_fn(results.parameters)) \
            or (filters_snipe is not None and any(f(results.parameters) for f in filters_snipe)):
                filtered_results.append(results)
            else:
                continue

    # ad hoc fix for simplifying labels for paper
    for results in filtered_results:
        if results.parameters is not None:
            if results.parameters.complex_type == ComplexType.BALL:
                results.pipeline_name = "PHR"
            elif results.parameters.complex_type == ComplexType.ELLIPSOID:
                results.pipeline_name = "PHE"



    def compute_mean_accuracy(result):
        accs = getattr(result, 'mean_accuracy_per_trnsf', {})
        return sum(accs.values()) / len(accs) if accs else float('-inf')


    filtered_results.sort(key=compute_mean_accuracy, reverse=True)
    filtered_results = add_aggregate_results_labels(filtered_results)

    aggregated_results_path = get_aggregated_results_basepath(setup_path)
    save_to_json([a_res.to_dict() for a_res in filtered_results], aggregated_results_path, add_timestamp=False)
    plot_aggregated_results(filtered_results, aggregated_results_path, print_parameters=False)

    latex_table = generate_latex_table([res for res in filtered_results if res.pipeline_name == "PHE"])
    with open (f"{aggregated_results_path}_latex_table.txt", "w") as f:
        f.write(latex_table)



if __name__ == '__main__':
    main()
