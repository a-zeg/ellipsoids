from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Optional, Any
from datetime import datetime
import os
import json
import hashlib

import numpy as np
import gudhi as gd



def spherisize(axes_lengths: np.ndarray, s:float=0):
    '''
    Returns axes lengths of a "spherisized" ellipsoid.
    For s=0, the function will output axes_lengths and for s=1 all elements of axes_lengths
    will be equal to the longest one.

    For example, for axes_ratios=np.array([3,2,1]), we get:
    - s=0: np.array([3,2,1])
    - s=0.5: np.array([3,2.5,2])
    - s=1: np.array([3,3,3])
    '''

    major_axis = np.max(axes_lengths)
    spherisized_axes_lengths = np.zeros(np.size(axes_lengths))
    for idx, axis in enumerate(axes_lengths):
        spherisized_axes_lengths[idx] = s*major_axis + (1-s)*axis
    return spherisized_axes_lengths



def scale_to_01(x: float, min: float = 0, max: float = 1):
    '''
    Scales the input x in the interval [min, max] to the interval [0,1]
    If x is bigger than max, returns 1.
    If x is smaller than min, returns 0.
    '''
    if x > max:
        return 1
    elif x < min:
        return 0
    else:
        return (x-min)/(max-min)



def spherisize_axes(axes_lengths: np.ndarray, r: float, r_spherisize: float = 10):
    '''
    At r=0 should have axes_ratio from the start.
    At r=r_spherisize, should get balls.
    '''

    return spherisize(axes_lengths, scale_to_01(r,max=r_spherisize))

def restore_numpy_array(field_data: Any) -> Optional[np.ndarray]:
    """
    Restores a NumPy array from a list (if necessary) or returns the original data.

    Args:
        field_data: The field data which might be a list or a NumPy array.

    Returns:
        np.ndarray: The restored NumPy array.
    """
    if field_data is None:
        return None
    if isinstance(field_data, list):
        return np.array(field_data)
    return field_data


class Ellipsoid:
    def __init__(self, center: np.ndarray, axes: np.ndarray, axes_lengths: np.ndarray):
        self.center = center
        self.axes = axes
        self.axes_lengths = axes_lengths
        self._cached_sigma = None

    def __eq__(self, other):
        if isinstance(other, Ellipsoid):
            return (np.array_equal(self.center, other.center) and
                    np.array_equal(self.axes, other.axes) and
                    np.array_equal(self.axes_lengths, other.axes_lengths))
        return False

    def get_sigma(self, r: float, r_spherisize: float):
        if r_spherisize == np.inf:
            if self._cached_sigma is None:
                self._cached_sigma = self.axes.T @ np.diag(self.axes_lengths ** 2) @ self.axes
            return self._cached_sigma
        else:
            lengths = spherisize_axes(self.axes_lengths, r, r_spherisize)
            return self.axes.T @ np.diag(lengths**2) @ self.axes

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



class TurkevsTransformation(Enum):
    STANDARD = ("std", "standard")
    TRANSLATION = ("trns", "translation")
    ROTATION = ("rot", "rotation")
    STRETCH = ("stretch", "stretch")
    SHEAR = ("shear", "shear")
    GAUSSIAN = ("gauss", "gaussian")
    OUTLIERS = ("out", "outliers")

    def __init__(self, shortname, fullname):
        self.shortname = shortname
        self.fullname = fullname

    @classmethod
    def from_shortname(cls, shortname: str):
        return next(t for t in cls if t.shortname == shortname)

    @classmethod
    def shortnames(cls):
        return [t.shortname for t in cls]



@dataclass
class TurkevsDatasetInfo:
    dataset_id: Optional[str]
    seed: int
    point_cloud_index: int
    transformation: TurkevsTransformation
    label: str                             # label for classification
    n_point_clouds: int

    def to_dict(self) -> dict:
        return {
            "dataset_id": self.dataset_id,
            "seed": self.seed,
            "point_cloud_index": self.point_cloud_index,
            "transformation": self.transformation.shortname,
            "label": self.label,
            "n_point_clouds": self.n_point_clouds
        }

    @classmethod
    def from_dict(cls, data: dict) -> "TurkevsDatasetInfo":
        return cls (
            dataset_id = data["dataset_id"],
            seed = data["seed"],
            point_cloud_index = data["point_cloud_index"],
            transformation = TurkevsTransformation.from_shortname(data["transformation"]),
            label = data["label"],
            n_point_clouds = data["n_point_clouds"]
        )

    def to_str(self) -> str:
        return f"turkevs_id={self.dataset_id}_point_cloud_index={self.point_cloud_index}_trnsf={self.transformation.shortname}"



@dataclass
class Dataset:
    points: np.ndarray
    data_type: str
    additional_info: Optional[Any] = None

    @property
    def n_points(self) -> int:
        return len(self.points)

    @property
    def ambient_dim(self) -> int:
        return len(self.points[0])

    def to_dict(self) -> dict:
        base = {
            "data_type": self.data_type,
            "points": self.points.tolist(),
        }
        if self.additional_info:
            base["additional_info"] = self.additional_info.to_dict()
        return base

    @classmethod
    def from_dict(cls, data: dict) -> "Dataset":
        points = restore_numpy_array(data["points"])
        data_type = data["data_type"]
        additional_info = None
        if "additional_info" in data and data_type=="turkevs":
            additional_info = TurkevsDatasetInfo.from_dict(data["additional_info"])
        return cls(points=points, data_type=data_type, additional_info=additional_info)



@dataclass
class DatasetSummary:
    data_type: str
    n_points: int
    additional_info: Optional[Any]

    def to_dict(self):
        return {
            "data_type": self.data_type,
            "n_points": self.n_points,
            "additional_info": self.additional_info.to_dict() if self.additional_info else None
        }

    @classmethod
    def from_dict(cls, data):
        data_type = data["data_type"]
        n_points = data.get("n_points", None)
        additional_info = TurkevsDatasetInfo.from_dict(data["additional_info"]) \
            if data.get("additional_info") and data_type=="turkevs" else None
        return cls(data_type=data_type, n_points = n_points, additional_info=additional_info)



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

    @staticmethod
    def from_any_dict(data: dict) -> "Parameters":
        complex_type = data.get("complex_type")
        if complex_type == ComplexType.ELLIPSOID.value:
            return EllipsoidParameters.from_dict(data)
        return Parameters.from_dict(data)

    def __eq__(self, other):
        if not isinstance(other, Parameters):
            return NotImplemented
        return self.to_dict() == other.to_dict()

    def __hash__(self):
        return hash(frozenset(self.to_dict().items()))



@dataclass
class EllipsoidParameters(Parameters):
    # WARNING: the first ComplexType in the next line is not just a type hint.
    # Without it, the default value won't be set correctly.
    complex_type: ComplexType = ComplexType.ELLIPSOID
    nbhd_size: int = 3
    axes_ratios: Optional[np.ndarray] = None
    r_spherisize: float = np.inf
    save_ellipsoid_list: bool = False
    use_cache: bool = True

    @property
    def use_pca_axes(self) -> bool:
        return self.axes_ratios is None

    def axes_ratios_to_str(self):
        return "pca" if self.use_pca_axes else str(self.axes_ratios)

    def to_filename(self) -> str:
        axes = self.axes_ratios_to_str()
        return f"nbhd_size={self.nbhd_size}_\
                 axes_ratios={axes}_\
                 r_spherisize={self.r_spherisize}_\
                 expansion_dim={self.expansion_dim}_\
                 {self.complex_type}"

    def to_dict(self):
        data = super().to_dict()
        data.update({
            "nbhd_size": self.nbhd_size,
            "axes_ratios": tuple(self.axes_ratios) if self.axes_ratios is not None else None,
            "r_spherisize": self.r_spherisize,
            "save_ellipsoid_list": self.save_ellipsoid_list,
            "use_cache": self.use_cache
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
            save_ellipsoid_list=data["save_ellipsoid_list"],
            use_cache=data.get("use_cache", None),
        )

    def __eq__(self, other):
        if not isinstance(other, Parameters):
            return NotImplemented
        return self.to_dict() == other.to_dict()

    def __hash__(self):
        return hash(frozenset(self.to_dict().items()))



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
    draw_points: bool = False
    draw_ellipsoids: bool = False
    draw_simplex_tree: bool = False
    filtration: float = 1
    n_bars: dict = field(default_factory=lambda: {0:5, 1:10, 2:10}) # number of bars in the reduced barcode
    x_axis_start: float = -0.05
    x_axis_end: float = 10

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict):
        return cls(
            draw_points=data.get("draw_points", False),
            draw_ellipsoids=data.get("draw_ellipsoids", False),
            draw_simplex_tree=data.get("draw_simplex_tree", False),
            filtration=data.get("r", 1),
            n_bars=data.get("n_bars", {0: 5, 1: 10, 2: 10}),
            x_axis_start=data.get("x_axis_start", -0.05),
            x_axis_end=data.get("x_axis_end", 10)
        )



@dataclass
class ExperimentSummary:
    dataset_summary: DatasetSummary
    parameters: Parameters
    barcode: list[tuple]
    execution_time: float

    def to_dict(self):
        return {
            "dataset_summary": self.dataset_summary.to_dict(),
            "parameters": self.parameters.to_dict(),
            "barcode": self.barcode,
            "execution_time": self.execution_time,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ExperimentSummary":
        return cls(
            dataset_summary=DatasetSummary.from_dict(data["dataset_summary"]),
            parameters=Parameters.from_any_dict(data["parameters"]),
            barcode=data["barcode"],
            execution_time=data["execution_time"],
        )

    def save_to_json(self, filename: str):
        from data_handling import save_to_json
        save_to_json(data=self.to_dict(), filename=filename, add_timestamp=False)

    def generate_filename(self):
        dataset_info = self.dataset_summary.data_type
        additional_info = self.dataset_summary.additional_info
        if isinstance(additional_info, TurkevsDatasetInfo):
            dataset_info = additional_info.to_str()
        filename = f"{dataset_info}_{self.parameters.complex_type}-{self.parameters.complex_subtype}"
        if isinstance(self.parameters, EllipsoidParameters):
            filename += f"_axes_ratios={self.parameters.axes_ratios_to_str()}"\
                f"_nbhd_size={self.parameters.nbhd_size}"\
                f"_r_s={self.parameters.r_spherisize}"
        return filename



class Experiment:
    def __init__(self, dataset: Dataset, parameters: Parameters):
        self.dataset = dataset
        self.parameters = parameters
        self.results: Results = Results()
        self.plot_parameters: PlotParameters = PlotParameters()

    def run(self):
        from ellipsoids.topological_computations import calculate_BALL_Results
        from ellipsoids.topological_computations import calculate_ELLIPSOID_Results

        match self.parameters.complex_type:
            case ComplexType.BALL:
                self.results = calculate_BALL_Results(self.dataset, self.parameters)
            case ComplexType.ELLIPSOID:
                if not isinstance(self.parameters, EllipsoidParameters):
                    raise TypeError(f"Expected EllipsoidParameters in '{self.run.__name__}', got {type(self.parameters).__name__}")
                self.results = calculate_ELLIPSOID_Results(self.dataset, self.parameters)

    def print_execution_time(self):
        if self.results is None:
            raise RuntimeError(f"Experiment with parameters {self.parameters} has not been run yet, no execution time to print.")
        print(f"Execution time of {self.parameters.complex_type}-{self.parameters.complex_subtype} is {self.results.execution_time}")

    def _generate_filename(self, add_timestamp=True):
       filename = f"{self.dataset.data_type}-{self.dataset.n_points}_{self.parameters.complex_type}-{self.parameters.complex_subtype}"
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

    def _generate_filepath(self, filename: Optional[str], folder="data", add_timestamp=True):
        if filename is None:
            filename = f"{self._generate_filename(add_timestamp=add_timestamp)}.json"
        return os.path.join(folder, filename)

    def save_to_json(self, folder="data", filename=None):
        from ellipsoids.data_handling import save_to_json
        if not self.results:
            raise RuntimeError("Experiment has not been run yet, no results to save.")
        filepath = self._generate_filepath(filename, folder)
        save_to_json(data=self.to_dict(), filename=filepath, add_timestamp=True)

    def serialize_summary(self):
        from ellipsoids.data_handling import CustomEncoder
        import json
        dataset_summary = DatasetSummary(
                data_type = self.dataset.data_type,
                n_points = self.dataset.n_points,
                additional_info = self.dataset.additional_info
            )
        experiment_summary = ExperimentSummary(
            dataset_summary = dataset_summary,
            parameters = self.parameters,
            barcode = self.results.barcode,
            execution_time = self.results.execution_time
        )
        return experiment_summary.to_dict()
        # data = experiment_summary.to_dict()
        # return json.dumps(data, cls=CustomEncoder, indent=2)

    def save_summary(self, folder):
        from ellipsoids.data_handling import save_to_json
        dataset_summary = DatasetSummary(
                data_type = self.dataset.data_type,
                n_points = self.dataset.n_points,
                additional_info = self.dataset.additional_info
            )
        experiment_summary = ExperimentSummary(
            dataset_summary = dataset_summary,
            parameters = self.parameters,
            barcode = self.results.barcode,
            execution_time = self.results.execution_time
            )
        filename = experiment_summary.generate_filename()
        filepath = os.path.join(folder,filename)
        save_to_json(data=experiment_summary.to_dict(), filename=filepath, add_timestamp=False)

    def get_results(self):
        return self.results

    @classmethod
    def read_from_json(cls, filename: str):
        from ellipsoids.data_handling import read_from_json
        json_dict = read_from_json(filename)
        return cls.from_dict(json_dict)

    @classmethod
    def from_dict(cls, data: dict) -> "Experiment":
        dataset = Dataset.from_dict(data["dataset"])
        parameters = Parameters.from_any_dict(data["parameters"])
        results = Results.from_dict(data["results"])

        experiment = cls(dataset, parameters)
        experiment.results = results
        experiment.plot_parameters = PlotParameters()
        return experiment



class ConversionType(Enum):
    RADIUS_TO_RIPS = "radius_to_RIPS_filtration"
    RIPS_TO_RADIUS = "RIPS_filtration_to_radius"
    RADIUS_TO_ALPHA = "radius_to_ALPHA_filtration"
    ALPHA_TO_RADIUS = "ALPHA_filtration_to_radius"



def convert(value: float, conversion_type: ConversionType):
    conversions = {
        ConversionType.RADIUS_TO_RIPS: lambda r: 2*r,
        ConversionType.RIPS_TO_RADIUS: lambda f: f/2,
        ConversionType.RADIUS_TO_ALPHA: lambda r: r**2,
        ConversionType.ALPHA_TO_RADIUS: lambda f: np.sqrt(f)
    }
    convert_func = conversions.get(conversion_type)
    if convert_func:
        return convert_func(value)
    else:
        raise ValueError(f"Invalid conversion type: {conversion_type}")









