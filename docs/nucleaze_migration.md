# BBDuk to Nucleaze Migration

This document describes the migration from BBDuk (Java-based, part of BBTools/BBMap suite) to Nucleaze (Rust-based drop-in replacement).

## Motivation

- **Performance**: Nucleaze eliminates JVM spin-up time, which is significant when processing many small files
- **Dependencies**: Removes Java dependency, replacing it with a single static binary
- **Memory**: More predictable memory usage without JVM overhead

## Installation

### Prerequisites

Nucleaze requires Rust to build from source. Install Rust via rustup:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

**Note**: Nucleaze requires Rust 1.88+ due to the `sysinfo` dependency requiring this version.

### Building Nucleaze

```bash
git clone https://github.com/jackdougle/nucleaze.git
cd nucleaze
cargo build --release
cargo install --path .
```

This installs the `nucleaze` binary to `~/.cargo/bin/`. Ensure this is in your PATH.

### Docker Container

A Dockerfile is provided at `nextflow/docker/Dockerfile.nucleaze`:

```bash
cd nextflow/docker
docker build -f Dockerfile.nucleaze -t harmonb/outward-assembly-nucleaze:v0.1 .
```

The Dockerfile uses a multi-stage build:
1. **Builder stage**: Uses `rust:latest` to compile nucleaze from source
2. **Runtime stage**: Uses `ubuntu:22.04` with only `zstd` installed

## OS-Specific Issues

### macOS

- **Memory limits**: The `--maxmem` flag does not work on macOS because it uses `setrlimit()` which is not supported. The flag was removed from all nucleaze invocations.
- **PATH inheritance**: Subprocesses spawned via `subprocess.run()` with `shell=True` may not inherit the full PATH. The code uses `shutil.which("nucleaze")` to find the binary and passes the full path.
- **Megahit segfault on Apple Silicon**: Megahit crashes with SIGSEGV (exit code -11) on ARM64 Macs. **Workaround**: Run the pipeline inside a Docker Linux container instead of natively on macOS.
- **Subprocess xargs issues**: Some subprocess calls with long argument lists fail on macOS. These were replaced with `ThreadPoolExecutor` from `concurrent.futures`.
- **Biopython not in oa-tools**: You may need to install biopython in the oa-tools environment: `mamba install biopython -n oa-tools`

### Linux

- **Stdin piping unreliable**: When running nucleaze with `--in stdin.fq` or `--in stdin` in a piped command (e.g., `cat file.fq | nucleaze --in stdin.fq ...`), nucleaze fails with "Error processing read sequences: No such file or directory (os error 2)". This occurs in Docker containers and possibly native Linux. **Workaround**: Write input to a temp file first, then pass the file path to nucleaze. See `_subset_split_files_local()` in `pipeline_steps.py`.

### Nucleaze Bug Report: Stdin Handling

**Issue**: Nucleaze v1.4.0 fails to read from stdin when using `--in stdin.fq` or `--in stdin` in piped commands.

**Environment**:
- Docker container (linux/amd64 on macOS Apple Silicon via Rosetta)
- Ubuntu 22.04
- Nucleaze 1.4.0

**To reproduce**:
```bash
# This fails with "Error processing read sequences: No such file or directory (os error 2)"
cat reads.fq | nucleaze --in stdin.fq --outm out_1.fq --outm2 out_2.fq \
  --outu /dev/null --outu2 /dev/null --ref ref.fa --k 21 \
  --canonical --minhits 1 --interinput --threads 2

# Using /dev/stdin also fails with "reads file is empty"
cat reads.fq | nucleaze --in /dev/stdin --outm out_1.fq --outm2 out_2.fq \
  --outu /dev/null --outu2 /dev/null --ref ref.fa --k 21 \
  --canonical --minhits 1 --interinput --threads 2

# This works (file input instead of stdin)
nucleaze --in reads.fq --outm out_1.fq --outm2 out_2.fq \
  --outu /dev/null --outu2 /dev/null --ref ref.fa --k 21 \
  --canonical --minhits 1 --interinput --threads 2
```

**Suspected cause**: Threading issue where stdin is consumed or closed before worker threads can read it, or the stdin detection logic fails in certain environments.

**Current workaround**: Write piped input to a temp file, run nucleaze on the file, then delete the temp file.

## Parameter Mapping

| BBDuk Parameter | Nucleaze Parameter | Notes |
|-----------------|-------------------|-------|
| `in=` | `--in` | Input file |
| `out=` | `--out` | Output matching reads |
| `outm=` | `--outm` | Output matching reads (same as --out) |
| `outu=` | `--outu` | Output non-matching reads |
| `ref=` | `--ref` | Reference k-mer file |
| `k=` | `--k` | K-mer size |
| `mm=f` (exact matching) | `--canonical` | Use canonical k-mers (forward + reverse complement) |
| `hdist=0` | (default) | Hamming distance 0 is default |
| `rcomp=t` | `--canonical` | Include reverse complement |
| `int=t` | `--interinput` | Interleaved input |
| `threads=` | `--threads` | Number of threads |
| `ordered=t` | `--order` | Preserve read order |
| `-Xmx` | `--maxmem` | Memory limit (removed due to macOS issues) |

## Code Changes Summary

### Files Modified

1. **`outward_assembly/pipeline_steps.py`**
   - 3 BBDuk usages migrated to Nucleaze
   - Added `shutil.which()` for PATH-safe binary lookup
   - Replaced `xargs` with `ThreadPoolExecutor` for parallel execution
   - Added edge case handling for short adapter sequences (< k)

2. **`outward_assembly/kmer_freq_filter.py`**
   - 1 BBDuk usage migrated to Nucleaze
   - Added `shutil.which()` for binary lookup
   - Replaced `xargs` with `ThreadPoolExecutor`

3. **`nextflow/main.nf`**
   - Renamed `BBDUK` process to `NUCLEAZE`
   - Updated container label

4. **`nextflow/config/containers.config`**
   - Updated container label and image reference

5. **`oa_tools_env.yml`**
   - Removed `bbmap>=39,<40` dependency
   - Added comment about nucleaze installation via cargo

6. **`nextflow/docker/Dockerfile.nucleaze`** (new file)
   - Multi-stage Docker build for nucleaze container

## Edge Cases and Fixes

### Short Adapter Sequences

When adapter sequences are shorter than the k-mer size (k), nucleaze cannot produce any k-mers from the reference. The code now checks for this and skips the filtering step:

```python
# Check if adapters file has sequences long enough to produce k-mers
has_valid_adapter_kmers = False
for record in SeqIO.parse(adapters_path, "fasta"):
    if len(record.seq) >= k:
        has_valid_adapter_kmers = True
        break

if not has_valid_adapter_kmers:
    # No adapter k-mers possible, just copy input to output
    shutil.copy(candidate_kmers_path, final_query_kmers_path)
    ...
```

### Command Line Length Limits

When processing many files, `xargs` can hit command line length limits. The code now uses Python's `ThreadPoolExecutor` for parallel execution:

```python
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_cmd(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True)

with ThreadPoolExecutor(max_workers=num_parallel) as executor:
    futures = [executor.submit(run_cmd, cmd) for cmd in cmds]
    for future in as_completed(futures):
        future.result()  # Raises exception if command failed
```

## Testing

Run the test suite to verify the migration:

```bash
cd /path/to/outward-assembly
pytest
```

Expected results: 70 tests pass. 5 tests may fail due to missing external tools (nextflow, fastp) that are not part of the nucleaze migration.

## Running the Pipeline Locally (Without Nextflow)

The local profile runs everything on your machine and does not require Nextflow or AWS Batch.

### Setup

1. **Install Python dependencies** (using uv):
   ```bash
   uv sync --extra dev
   ```

2. **Create the tools environment** (using mamba/conda):
   ```bash
   mamba env create -n oa-tools -f oa_tools_env.yml --channel-priority flexible
   ```

3. **Install Nucleaze** (requires Rust 1.88+):
   ```bash
   git clone https://github.com/jackdougle/nucleaze.git
   cd nucleaze
   cargo build --release
   cargo install --path .
   ```

4. **Activate tools environment before running**:
   ```bash
   mamba activate oa-tools
   ```

### Running Outward Assembly

The primary entrypoint is the `outward_assembly` function in `outward_assembly/pipeline.py`. Example usage:

```python
from outward_assembly.pipeline import outward_assembly
from outward_assembly.io_helpers import s3_files_with_prefix

# Get list of S3 paths to your read files (SIZ format)
prefixes = [
    "outward-assembly-test-data/siz/simulated-abcbd-reads_1",
    "outward-assembly-test-data/siz/simulated-abcbd-reads_2",
]
paths = [p for prefix in prefixes for p in s3_files_with_prefix("nao-testing", prefix)]

# Run outward assembly locally
outward_assembly(
    seed_path="path/to/seed.fasta",
    input_s3_paths=paths,
    out_path="output/contigs.fasta",
    work_dir_parent="workdir",
    use_batch=False,  # Local profile - no Nextflow
)
```

### Running with Automation (YAML Config)

For iterative assembly with automation:

1. Create a dataset list CSV:
   ```csv
   s3_path,priority
   s3://bucket/data_div0001.fastq.zst,1
   s3://bucket/data_div0002.fastq.zst,1
   ```

2. Create a YAML config file:
   ```yaml
   assembly:
     input_seed_path: /path/to/seed.fasta
     input_dataset_list: /path/to/datasets.csv
     dataset_priority: 1
     work_dir: /path/to/workdir
     out_dir: /path/to/output
     output_filename: contigs.fasta
     read_subset_k: 27
     use_batch: false  # Local profile

   decision:
     automate: false  # Set to true for iterative automation
   ```

3. Run:
   ```bash
   mamba activate oa-tools
   uv run python automate_assembly.py --input_config config.yaml
   ```

### Running Tests

```bash
mamba activate oa-tools
uv run pytest
```

Expected results with nucleaze migration: 70 tests pass, 5 fail (due to missing external tools: nextflow, fastp).

## Performance Benchmark Results

### K-mer Filter Microbenchmark

Testing k-mer filtering only (941 read pairs):

| Tool      | Time (s) | Output Reads | Output Hash                      |
|-----------|----------|--------------|----------------------------------|
| Nucleaze  | 0.009    | 915          | c4c04e3d253c8b0cdc6ece1f10d003f4 |
| BBDuk     | 0.589    | 915          | c4c04e3d253c8b0cdc6ece1f10d003f4 |

**Result**: Nucleaze is ~65x faster for isolated k-mer filtering.

### Full Pipeline Benchmark (Docker)

Testing full outward assembly pipeline with megahit (1882 reads):

| Scenario   | BBDuk (main) | Nucleaze (migration) | Output Hash                      |
|------------|-------------|---------------------|----------------------------------|
| basic      | **24.16s**  | 27.68s              | 1bf62f2ba03350f5eeab17ebb745f3e6 |
| multi_seed | **24.31s**  | 26.94s              | 2e5382aab22bb60019e557b5b59d59d7 |

**Key findings**:
1. **Output is identical** - Both tools produce exactly the same results
2. **BBDuk slightly faster in this test** due to the temp file workaround for nucleaze's stdin bug adding ~3s I/O overhead
3. **Assembly time dominates** - megahit takes most of the time, making k-mer filtering differences less significant
4. **Expected with larger data**: Nucleaze should still win with many parallel invocations where JVM startup overhead compounds

Run benchmarks yourself:
```bash
# Build Docker image
docker build --platform linux/amd64 -f benchmarks/Dockerfile.benchmark -t oa-benchmark .

# Run full pipeline benchmark
docker run --platform linux/amd64 -it --rm \
  -v $(pwd):/workspace \
  -v ~/.aws:/root/.aws:ro \
  oa-benchmark \
  bash -c "cd /workspace && uv sync --extra dev && uv run python benchmarks/run_benchmark.py"
```

## Additional Setup Notes

### AWS CLI Required

Even when not using batch/Nextflow mode, AWS CLI is required to download test data from S3:
```bash
# Install AWS CLI
brew install awscli  # macOS
# or pip install awscli

# Configure credentials
aws configure
```

### Conda Shell Initialization

If `conda` and `mamba` commands don't work in your shell, add this to your `~/.zshrc` or `~/.bashrc`:
```bash
eval "$(/Users/lee/miniforge3/bin/conda shell.zsh hook)"
```

Or run this before using conda:
```bash
source /Users/lee/miniforge3/etc/profile.d/conda.sh
```

## TODO

- [ ] Publish Docker container to registry
- [ ] Pin nucleaze to a specific version/tag once available
- [ ] Update `docs/installation.md` to reference Nucleaze instead of BBMap/BBDuk
- [ ] Add quickstart script example that imports and uses the pipeline
- [ ] Document Docker workflow for macOS users (to work around Megahit crash)
