# Experimental Features

This document describes experimental features of outward assembly. These features are:
- Not fully tested in production use
- Subject to change without notice
- May have incomplete documentation

For standard outward assembly usage, see [usage.md](usage.md).

---

## Pileup Visualization

Visualizes how paired-end reads align to an assembled contig, highlighting a seed sequence of interest.

![Example pileup visualization](pileup_visualization/tiny_pileup.png)

<details>
<summary>More examples</summary>

**40 reads:**
![Medium pileup example](pileup_visualization/medium_pileup.png)

**150 reads:**
![Large pileup example](pileup_visualization/big_pileup.png)

</details>

### Requirements

```bash
uv sync --extra pileup  # Installs Pillow and pysam for pileup visualization
```

### Inputs

| Argument | Description |
|----------|-------------|
| `--fwd-fastq` | Forward reads (R1) FASTQ file |
| `--rev-fastq` | Reverse reads (R2) FASTQ file |
| `--contigs-fasta` | Assembled contigs FASTA file |
| `--contig-name` | Name of the contig to visualize |
| `--seed-sequence` | DNA sequence to highlight (e.g., the assembly seed) |
| `--output` | Output PNG path |

### Usage

```bash
python -m outward_assembly.pileup \
    --fwd-fastq reads_R1.fastq \
    --rev-fastq reads_R2.fastq \
    --contigs-fasta contigs.fasta \
    --contig-name "contig_1" \
    --seed-sequence "ATCGATCGATCG" \
    --output output.png
```

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--scale-x` | 2 | Horizontal pixel scale factor |
| `--scale-y` | 2 | Vertical pixel scale factor |
| `--padding` | 10 | Padding around image elements (pixels) |
| `--contig-height` | 10 | Height of contig ribbon (pixels) |
| `--threads` | 4 | Threads for bowtie2 alignment |

### Output

A PNG image showing:
- **Top ribbon**: The contig sequence (blue = non-seed, orange = seed region)
- **Stacked rows**: Each row is one read pair, sorted by alignment position
- **Left/right sidebars**: Magenta markers indicating insertions

#### Color Legend

| Color | Meaning |
|-------|---------|
| Blue | Matching bases (read pair contains entire seed) |
| Muted blue | Matching bases (read pair does not contain entire seed) |
| Orange | Matching bases within the seed region |
| Red | Mismatched bases |
| Black | Deletions (bases in contig but not in read) |
| Cyan | Soft-clipped bases |
| Pink line | Insertion marker (inline) |
| Magenta | Insertion count (edge sidebars) |
| Light grey | Unsequenced insert (gap between mates) |

#### Insertion and Deletion Semantics

The **contig** is the reference sequence. Each column in the visualization corresponds to one contig position.

- **Deletion**: Bases present in the contig but absent from the read. These are shown as black pixels inline, since they occupy contig positions.

- **Insertion**: Bases present in the read but absent from the contig. Since there's no contig column for these bases, insertions can't be shown inline at full resolution.

**How insertions are displayed:**

| Location | Appearance | Meaning |
|----------|------------|---------|
| **Inline** | Thin pink vertical line at right edge of base | Marks the position *after which* an insertion occurred. The base color is preserved; the pink line indicates "an insertion starts here." The size of the insertion is not shown. |
| **Edge sidebars** | Solid magenta pixels | The left/right margins show the *total count* of inserted bases for the upstream/downstream mate. |

### Python API

```python
from outward_assembly.pileup import create_pileup_visualization

image = create_pileup_visualization(
    fwd_fastq="reads_R1.fastq",
    rev_fastq="reads_R2.fastq",
    contigs_fasta="contigs.fasta",
    contig_name="contig_1",
    seed_sequence="ATCGATCGATCG",
    output_path="output.png",  # Optional: saves to file
    scale=(2, 2),
    padding=10,
    contig_height=10,
    threads=4
)
# Returns: numpy array of shape (height, width, 4) with RGBA values
```


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
  use_batch: <true to use batch, false to use local> (default is false)

decision: (optional)
  automate: <true to use automation, false to run outward assembly once> (default is false)
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

When running in *automated* mode, a new working directory is created each time the pipeline repeats.
