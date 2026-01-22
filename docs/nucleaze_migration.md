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

## Known Issues and Workarounds

### macOS megahit multi-threading bug

MEGAHIT v1.2.9 crashes with a segmentation fault (exit code -11) on macOS when using `--num-cpu-threads > 1`. This affects both Intel and Apple Silicon Macs.

**Upstream issue**: https://github.com/voutcn/megahit/issues/385
**Potential fix**: https://github.com/voutcn/megahit/pull/387

**Workaround**: The pipeline automatically forces `--num-cpu-threads 1` on macOS. See `_assemble_contigs()` in `pipeline_steps.py`. This makes assembly slower but functional.

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

**TODO**: Migrate some of the above into the installation docs

### Running Tests

```bash
mamba activate oa-tools
uv run pytest
```

Expected results with nucleaze migration: 70 tests pass, 5 fail (due to missing external tools: nextflow, fastp).

## Performance Benchmark Results

### K-mer Filter Microbenchmark

K-mer filtering is the step where Nucleaze **replaces BBDuk**. It finds reads sharing k-mers with the reference sequence. All downstream steps (assembly with megahit, contig subsetting) use unchanged code. Therefore, identical k-mer filter output guarantees identical final pipeline output.

Testing k-mer filtering (941 read pairs, macOS):

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
2. **Assembly time dominates** - megahit takes most of the time, making k-mer filtering differences less significant
3. **Expected with larger data**: Nucleaze should win with many parallel invocations where JVM startup overhead compounds

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