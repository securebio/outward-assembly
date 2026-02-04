# Experimental Features

This document describes experimental features of outward assembly. These features are:
- Not fully tested in production use
- Subject to change without notice
- May have incomplete documentation

For standard outward assembly usage, see [usage.md](usage.md).

---

## Automated Mode

*Automated* mode is designed for users who need to run outward assembly iteratively. This could include users who are unsure which samples contain reads that will successfully assemble with their seed sequence, and therefore want to search through their data iteratively.

With *automated* mode, you can strategically begin with a subset of your data and progressively expand your search if initial assembly results are insufficient. The system intelligently adjusts parameters between iterations based on your defined strategy, optimizing both computational resources and discovery potential.

*Warning: automated mode does not support multiline sequences in FASTQ files.*

### Usage

The primary entrypoint to *automated* outward assembly is `automate_assembly.py`. This script takes in a yaml file as input, which specifies the parameters used for running the pipeline.

*Disclaimer: right now, the `automate_assembly.py` script does not implement all the parameters of `outward_assembly`.*

Generally, using this script will look like the following:
1. Prepare a list of datasets along with priorities for each of them.
2. (Optional) Add your automation strategy to `outward_assembly/strategy.py`
   * The user may also decide to use this script without the automation turned on
3. Write your configuration in a yaml file.
4. Run `automate_assembly.py` and pass in your yaml file.

#### 1. Prepare a list of datasets along with priorities for each of them.
The YAML configuration file requires a list of S3 paths to datasets in SIZ format, each with an assigned priority.

This prioritization system allows you to begin with a smaller dataset subset and progressively include more data if initial assembly results are insufficient. When automation is enabled, the pipeline can automatically advance to datasets with the next priority level if the current assembly results don't meet the criteria defined in your strategy. Importantly, the pipeline will only look at the datasets within the current priority level, so if you want data in earlier priority levels to be considered, make sure to include them in the current priority level (in practice, this means that you will have the same data in multiple priority levels).

The input file should follow this CSV format:

```csv
s3_path,priority # Header
s3://random-test-data/here-is-some-data.fastq.zst,1 # Example row
```
#### 2. (Optional) Add your automation strategy to `outward_assembly/strategy.py`
Determining when to initiate another round of outward assembly can be challenging. While we plan to introduce a default automation configuration in the future, none currently exists. Users can create their own configurations by writing Python functions.

Users may utilize any of the variables outlined in `outward_assembly/strategy.py` to define conditions and actions for adjusting parameters between iterations. Once a strategy is created, users can reference the function name in their configuration file, allowing the pipeline to automatically import it.

An example strategy, `example_strategy`, is provided in `outward_assembly/strategy.py`.

*We are working to expand the number of conditions and actions available for automating outward assembly. In the future, we aim to include a default strategy for users who prefer not to create their own.*
#### 3. Write your configuration in a yaml file.
The YAML file is the primary configuration file for the outward assembly pipeline. It includes parameters for the assembly process, the automation strategy, and any compute restrictions. Below, all supported parameters are listed along with their defaults and whether they are optional.

```yaml
assembly:
  input_seed_path: <Path to seed sequence>
  input_dataset_list: <Path to dataset list from step 1>
  dataset_priority: <Dataset priority to start with> (default is 1)
  adapter_path: <Path to adapter sequence> (optional; default is None)
  work_dir: <Path to work directory>
  out_dir: <Path to output directory>
  output_filename: <Filename of output contigs>
  read_subset_k: <Kmer size for BBDuk> (default is 27)
  use_batch: <TRUE to use batch, False to use local> (default is FALSE)

decision: (optional)
  automate: <TRUE to use automation, FALSE to run outward assembly once> (default is FALSE)
  strategy: <Name of strategy to use> (required if automate is TRUE)
  limits: (optional)
    compute_time: <Max compute time in hours> (default is 5 hours)
    iterations: <Max iterations> (default is 20 iterations)
```

*Note: The `decision` parameter is optional. If not specified, the pipeline will run outward assembly once.*

#### 4. Run `automate_assembly.py` and pass in your yaml file.

The automation script can be run by executing the following command:

```bash
uv run automate_assembly.py --input_config <PATH TO YAML FILE>
```

### Working Directory Note

When running in *automated* mode, a new working directory is created each time the pipeline repeats. See [usage.md](usage.md#working-directory-structure) for the standard working directory structure.
