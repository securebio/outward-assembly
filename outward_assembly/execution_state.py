from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any, Dict

import yaml

DEFAULT_COMPUTE_TIME_LIMIT = "5 hours"
DEFAULT_MAX_OUTER_ITERATIONS = 20
DEFAULT_DATASET_PRIORITY = 1
DEFAULT_READ_SUBSET_K = 27
DEFAULT_USE_BATCH = False


@dataclass
class ExecutionState:
    """Represents the current state of outer iterations of outward assembly, containing parameters for outward assembly that are fixed and configurable, as well as limits on the number of outer iterations that may occur based on the user defined parameters.

    Attributes:
        input_seed_path (str): Path to the input seed file.
        input_dataset_list (str): List of input datasets.
        adapter_path (str): Path to the adapter.
        output_path (str): Path for the output file.
        use_batch (bool): Flag indicating whether to use batch processing.
        work_dir (str): Working directory for the execution.
        read_subset_k (int): Selection of k for read filtering
        dataset_priority (int): Dataset priority for the current iteration.
        max_compute_time (timedelta): Maximum allowed compute time.
        max_outer_iterations (int): Maximum number of times Outward Assembly will run.
        automate (bool): Flag indicating whether to automate the assembly process.
        strategy (str): Strategy for the assembly process.
        start_time (datetime): Timestamp when the first Outward Assembly run started.
        current_outer_iterations (int): Counter for the current number of iterations that Outward Assembly has been repeated.
    """

    # General config
    input_seed_path: str
    input_dataset_list: str
    adapter_path: str
    output_dir: str
    output_filename: str
    use_batch: bool
    work_dir: str

    # Current state
    read_subset_k: int
    dataset_priority: int

    # Limits
    max_compute_time: timedelta
    max_outer_iterations: int

    # Automation related variables; these are only kept for the purpose of writing the config file
    automate: bool
    strategy: str

    # Tracking variables
    start_time: datetime
    current_outer_iterations: int = 0

    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> "ExecutionState":
        # Parse limits from config
        decision = config.get("decision", {})
        limits = decision.get("limits", {})
        assembly = config.get("assembly", {})

        # Convert time string to timedelta
        time_str = limits.get("compute_time", DEFAULT_COMPUTE_TIME_LIMIT)
        hours = float(time_str.split()[0])

        # Validate config file
        if decision.get("automate", True) and decision.get("strategy", ""):
            raise ValueError("Strategy must be specified if automate is True")

        return cls(
            input_seed_path=assembly.get("input_seed_path"),
            input_dataset_list=assembly.get("input_dataset_list"),
            adapter_path=assembly.get("adapter_path", None),
            output_dir=assembly.get("out_dir"),
            output_filename=assembly.get("output_filename"),
            use_batch=assembly.get("use_batch", DEFAULT_USE_BATCH),
            work_dir=assembly.get("work_dir"),
            read_subset_k=assembly.get("read_subset_k", DEFAULT_READ_SUBSET_K),
            dataset_priority=assembly.get("dataset_priority", DEFAULT_DATASET_PRIORITY),
            max_compute_time=timedelta(hours=hours),
            max_outer_iterations=limits.get("iterations", DEFAULT_MAX_OUTER_ITERATIONS),
            start_time=datetime.now(),
            automate=decision.get("automate", False),
            strategy=decision.get("strategy", None),
        )

    def get_adjustable_parameters(self) -> Dict[str, Any]:
        """Convert variables that the user may want to change between outward assembly runs into a dictionary."""
        return {
            "read_subset_k": self.read_subset_k,
            "dataset_priority": self.dataset_priority,
        }

    def update_adjustable_parameters(self, updated_state: Dict[str, Any]) -> None:
        """Update the execution state based on the actions that were sucessfully triggered due to the user's decision rules being met."""
        self.read_subset_k = updated_state.get("read_subset_k", self.read_subset_k)
        self.dataset_priority = updated_state.get("dataset_priority", self.dataset_priority)
        self.current_outer_iterations += 1

    def check_limits(self) -> bool:
        """Check if any limits have been exceeded."""
        if self.current_outer_iterations >= self.max_outer_iterations:
            return True

        elapsed_time = datetime.now() - self.start_time
        if elapsed_time > self.max_compute_time:
            return True

        return False

    def to_yaml_config(self) -> Dict[str, Any]:
        """Convert the execution state to a YAML-compatible dictionary format matching config.yaml structure."""
        return {
            "assembly": {
                "input_seed_path": self.input_seed_path,
                "input_dataset_list": self.input_dataset_list,
                "adapter_path": self.adapter_path,
                "work_dir": self.work_dir,
                "out_dir": self.output_dir,
                "output_filename": self.output_filename,
                "read_subset_k": self.read_subset_k,
                "use_batch": self.use_batch,
                "dataset_priority": self.dataset_priority,
            },
            "decision": {
                "automate": self.automate,
                "strategy": self.strategy,
                "limits": {
                    "compute_time": f"{self.max_compute_time.total_seconds() / 3600} hours",
                    "iterations": self.max_outer_iterations,
                },
            },
        }

    def write_config(self, yaml_path: str) -> None:
        """Write the current execution state to a YAML configuration file.

        Args:
            yaml_path: Path where the YAML file should be written
        """
        config = self.to_yaml_config()
        with open(yaml_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
