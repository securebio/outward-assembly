# Outward assembly
Welcome to the outward assembly repository! This repo provides Python tooling to assemble the genomic context around provided _seed sequences_ from a collection of read pairs. The seed sequences could be flagged chimeric junctions, kmers identified by reference-free growth detection, etc.

When the collection of reads is very large -- e.g. billions of reads from one delivery, or tens of billions of reads across deliveries -- doing a full metagenomic joint assembly is slow and computationally expensive. Instead of jointly assembling all reads, outward assembly attempts to iteratively grow contigs outward from the provided seeds. The basic algorithm is simple:
1. Start with the seeds as your initial contig.
2. Iteratively:
    * Find all read pairs that share a kmer with your contigs.
    * Assemble these read pairs.
    * Filter contigs to those that contain a seed.
3. Continue until either maximum iterations are reached or the assembly algorithm does not make progress from one iteration to another ("convergence").

Although the basic algorithm is simple, in practice, getting good assembly results is a bit more complex (see [algorithm details](docs/algorithm_details.md) for how we handle some of these complexities).

## Quick Start

### Installation

On an EC2 instance running Linux:

1. Install dependencies: `uv sync --extra dev`
2. Create tools environment: `mamba env create -n oa-tools -f oa_tools_env.yml --channel-priority flexible`
3. Activate tools environment: `mamba activate oa-tools`

See [installation docs](docs/installation.md) for more details.

### Running outward assembly

Outward assembly's main entry point is the Python function `outward_assembly`, and a command line interface is not provided.
Typically you'll run outward assembly by executing a short Python script ([example](docs/example-assembly-dir/script.py)) which configures and calls OA.
For example, to invoke this example script: `uv run docs/example-assembly-dir/script.py`.

The [usage docs](docs/usage.md) contain more details and guidance.

## Documentation

- [Installation](docs/installation.md)
- [Usage](docs/usage.md)
- [Algorithm details](docs/algorithm_details.md)  
- [Changelog](CHANGELOG.md)
- [Experimental features](docs/experimental.md)
