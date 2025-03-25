from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Optional
from datetime import datetime
import json
import os

import numpy as np
import gudhi as gd




def restore_numpy_array(field_data: any) -> np.ndarray:
    """
    Restores a NumPy array from a list (if necessary) or returns the original data.

    Args:
        field_data: The field data which might be a list or a NumPy array.

    Returns:
        np.ndarray: The restored NumPy array.
    """
    if isinstance(field_data, list):
        return np.array(field_data)
    return field_data


class Ellipsoid:
    def __init__(self, center: np.ndarray, axes: np.ndarray, axes_lengths: np.ndarray):
        self.center = center
        self.axes = axes
        self.axes_lengths = axes_lengths

    def __eq__(self, other):
        if isinstance(other, Ellipsoid):
            return (np.array_equal(self.center, other.center) and
                    np.array_equal(self.axes, other.axes) and
                    np.array_equal(self.axes_lengths, other.axes_lengths))
        return False

    def to_dict(self):
        obj_data = {
                "center": self.center,
                "axes": self.axes,
                "axes_lengths": self.axes_lengths
            }
        return obj_data

    @classmethod
    def from_dict(cls, data: dict):
        center = restore_numpy_array(data.get("center", None))
        axes = restore_numpy_array(data.get("axes", None))
        axes_lengths = restore_numpy_array(data.get("axes_lengths", None))
        return cls(
            center=center,
            axes=axes,
            axes_lengths=axes_lengths
        )



@dataclass
class Dataset:
    points: np.ndarray
    data_type: str

    def n_points(self) -> int:
        return len(self.points)

    def ambient_dim(self) -> int:
        return len(self.points[0])

    def to_dict(self) -> dict:
        return {
            "data_type": self.data_type,
            "points": self.points.tolist(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Dataset":
        points = restore_numpy_array(data["points"])
        return cls(points=points, data_type=data["data_type"])



class ComplexType(Enum):
    BALL = "BALL"
    ELLIPSOID = "ELLIPSOID"

    def __str__(self):
        return self.name.upper()



class ComplexSubtype(Enum):
    RIPS = "RIPS"
    ALPHA = "ALPHA"

    def __str__(self):
        return self.name.title()



@dataclass
class Parameters:
    """ for storing the calculation parameters """
    expansion_dim: int = 2
    collapse_edges: bool = True
    complex_type: ComplexType = ComplexType.BALL
    complex_subtype: ComplexSubtype = ComplexSubtype.RIPS
    save_simplex_tree: bool = False

    def to_dict(self):
        return {
            "expansion_dim": self.expansion_dim,
            "collapse_edges": self.collapse_edges,
            "complex_type": self.complex_type.value,
            "complex_subtype": self.complex_subtype.value,
            "save_simplex_tree": self.save_simplex_tree,
        }

    @classmethod
    def from_dict(cls, data):
        return cls(
            expansion_dim=data["expansion_dim"],
            collapse_edges=data["collapse_edges"],
            complex_type=ComplexType(data["complex_type"]),
            complex_subtype=ComplexSubtype(data["complex_subtype"]),
            save_simplex_tree=data["save_simplex_tree"],
        )



@dataclass
class EllipsoidParameters(Parameters):
    """ for storing all the parameters for ellipsoid complexes """

    # WARNING: the first ComplexType in the next line is not just a type hint.
    # Without it, the default value won't be set correctly.
    complex_type: ComplexType = ComplexType.ELLIPSOID
    nbhd_size: int = 3
    axes_ratios: np.ndarray = field(default_factory=lambda: np.array([2,1]))
    r_spherisize: float = np.inf
    save_ellipsoid_list: bool = True

    def to_filename(self) -> str:
        return f"nbhd_size={self.nbhd_size}_\
                 axes_ratios={self.axes_ratios}_\
                 r_spherisize={self.r_spherisize}_\
                 expansion_dim={self.expansion_dim}_\
                 {self.complex_type}"

    def to_dict(self):
        data = super().to_dict()
        data.update({
            "nbhd_size": self.nbhd_size,
            "axes_ratios": self.axes_ratios,
            "r_spherisize": self.r_spherisize,
            "save_ellipsoid_list": self.save_ellipsoid_list,
        })
        return data

    @classmethod
    def from_dict(cls, data: dict):
        # Deserialize the parent class fields first
        parent_data = {key: data[key] for key in Parameters.__annotations__ if key in data}
        parameters = Parameters.from_dict(parent_data)

        axes_ratios = restore_numpy_array(data.get("axes_ratios", None))  # Default to None if missing

        return cls(
            **parameters.__dict__,
            nbhd_size=data["nbhd_size"],
            axes_ratios=axes_ratios,
            r_spherisize=data["r_spherisize"],
            save_ellipsoid_list=data["save_ellipsoid_list"]
        )



@dataclass
class Results:
    barcode: Optional[list[tuple]] = field(default_factory=list)
    simplex_tree: Optional[list] = field(default_factory=gd.SimplexTree)
    execution_time: Optional[float] = None

    def to_dict(self):
        # return {key: value for key, value in vars(self).items() if value not in (None, [], {})}
        return {
            "barcode": self.barcode,
            "simplex_tree": self.simplex_tree,
            "execution_time": self.execution_time,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Results":
        return cls(
            barcode = data.get("barcode"),
            simplex_tree = data.get("simplex_tree"),
            execution_time = data.get("execution_time"),
            )



@dataclass
class EllipsoidResults(Results):
    ellipsoid_list: Optional[list[Ellipsoid]] = field(default_factory=list)

    def to_dict(self):
        data = super().to_dict()
        data["ellipsoid_list"] = self.ellipsoid_list
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "EllipsoidResults":
        results = super().from_dict(data)
        ellipsoid_list = data.get("ellipsoid_list", [])

        return cls(
            **results.__dict__,
            ellipsoid_list = ellipsoid_list
        )




@dataclass
class PlotParameters:
    """ for storing the plot parameters """
    draw_points: bool = False
    draw_ellipsoids: bool = False
    draw_simplex_tree: bool = False
    r: float = 1
    n_bars: dict = field(default_factory=lambda: {0:5, 1:10, 2:10}) # number of bars in the reduced barcode
    x_axis_start: float = -0.05
    x_axis_end: float = 10

    def to_dict(self):
        return asdict(self)

    # @classmethod
    # def from_dict(cls, data: dict):
    #     return cls(**data)

    @classmethod
    def from_dict(cls, data: dict):
        # Use .get() for optional fields and provide defaults
        return cls(
            draw_points=data.get("draw_points", False),
            draw_ellipsoids=data.get("draw_ellipsoids", False),
            draw_simplex_tree=data.get("draw_simplex_tree", False),
            r=data.get("r", 1),
            n_bars=data.get("n_bars", {0: 5, 1: 10, 2: 10}),
            x_axis_start=data.get("x_axis_start", -0.05),
            x_axis_end=data.get("x_axis_end", 10)
        )





class Experiment:
    def __init__(self, dataset: Dataset, parameters: Parameters):
        self.dataset = dataset
        self.parameters = parameters
        self.results: Results = Results()
        self.plot_parameters: PlotParameters = PlotParameters()
        # self.turkevs_parameters: Optional[TurkevsParameters] = None

    def run(self):
        from ellipsoids.topological_computations import calculate_rips
        from ellipsoids.topological_computations import calculate_ellipsoids

        match self.parameters.complex_type:
            case ComplexType.BALL:
                self.results = calculate_rips(self.dataset, self.parameters)
            case ComplexType.ELLIPSOID:
                if not isinstance(self.parameters, EllipsoidParameters):
                    raise TypeError(f"Expected EllipsoidResults in '{self.run.__name__}', got {type(self.parameters).__name__}")
                self.results = calculate_ellipsoids(self.dataset, self.parameters)

    def print_execution_time(self):
        if self.results is None:
            raise RuntimeError(f"Experiment with parameters {self.parameters} has not been run yet, no execution time to print.")
        print(f"Execution time of {self.parameters.complex_type}-{self.parameters.complex_subtype} is {self.results.execution_time}")

    def _generate_filename(self, add_timestamp=True):
       filename = f"{self.dataset.data_type}-{self.dataset.n_points()}_{self.parameters.complex_type}-{self.parameters.complex_subtype}"
       if add_timestamp:
           timestamp = datetime.now().strftime("%Y%m%d_%H%M%S%f")
           filename = f"{filename}__{timestamp}"
       return filename

    def to_dict(self):
        experiment_data = {
            'dataset': self.dataset.to_dict(),
            'parameters': self.parameters.to_dict(),
            'results': self.results.to_dict(),
        }
        return experiment_data

    def _generate_filepath(self, folder, filename):
        from ellipsoids.data_handling import ensure_folder_exists
        ensure_folder_exists(folder)
        if filename is None:
            filename = f"{self._generate_filename()}.json"
        return os.path.join(folder, filename)


    def save_to_json(self, folder="data", filename=None):
        from ellipsoids.data_handling import CustomEncoder

        if not self.results:
            raise RuntimeError("Experiment has not been run yet, no results to save.")

        if filename == None:
            filename = self._generate_filename()

        experiment_dict = self.to_dict()
        json_string = json.dumps(experiment_dict, cls=CustomEncoder, indent=4)

        filepath = self._generate_filepath(folder, filename)
        if not filepath.endswith('.json'):
            filepath = f"{filepath}.json"

        with open(filepath, 'w') as outfile:
            outfile.write(json_string)
        print(f"Experiment data saved to {filepath}")

    def get_results(self):
        return self.results

    @classmethod
    def read_from_json(cls, filename: str):
        from ellipsoids.data_handling import read_variables
        json_dict = read_variables(filename)
        return cls.from_dict(json_dict)

    @classmethod
    def from_dict(cls, data: dict) -> "Experiment":
        dataset = Dataset.from_dict(data["dataset"])
        parameters = Parameters.from_dict(data["parameters"])
        results = Results.from_dict(data["results"])
        # plot_parameters = PlotParameters.from_dict(data["plot_parameters"])

        experiment = cls(dataset, parameters)
        experiment.results = results
        return experiment
