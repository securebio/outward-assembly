# Usage

This document covers how to run outward assembly.

Table of Contents:
1. [Overview](#overview)
2. [Running Outward Assembly](#running-outward-assembly)
3. [Working directory structure](#working-directory-structure)
4. [Tips](#tips)
5. [(Optional) AWS Batch for read search](#optional-aws-batch-for-read-search)

> **Note:** For experimental features including automated/iterative mode, see [experimental features](experimental.md).

## Overview

The primary entrypoint to outward assembly is the Python function `outward_assembly` in `outward_assembly/pipeline.py`; a command line interface does not (yet) exist. See the docstring of `outward_assembly` for a detailed description of keyword parameters. The required parameters are:
* The path to a seed sequence (fasta) to assemble outward from;
* A list of s3 paths of reads to assemble. Reads must be in [SIZ format](./algorithm_details.md#input-data);
* Path for output contigs.

By default, the entire pipeline runs on your local machine. This works well for most use cases; if you need more resources, try scaling up your machine (e.g. pick an EC2 instance with more cores and memory). For very large datasets, there's an option to [run the read search step on AWS Batch](#optional-aws-batch-for-read-search), trading startup time for potentially very wide parallelization.

Outward assembly generates lots of intermediate results in its working directory; these are especially useful if you want to monitor the progress of an assembly in progress, retrieve results from intermediate iterations, or debug a crash. See the `outward_assembly` arguments `work_dir_parent` and `cleanup`.

## Running Outward Assembly

Running outward assembly requires calling the `outward_assembly` function with the appropriate parameters. The only setup that the user needs to do is to generate a list of s3 paths to their reads.

We've provided the function `s3_files_with_prefix` in `io_helpers.py` to assist you in generating a list of s3 paths. E.g. to get paths to split files from three demux sets:

```python
prefixes = [
	"outward-assembly-test-data/siz/simulated-abcbd-reads_1",
	"outward-assembly-test-data/siz/simulated-abcbd-reads_2",
]
paths = [p for prefix in prefixes for p in s3_files_with_prefix("nao-testing", prefix)]
```

where `outward-assembly-test-data/siz/simulated-abcbd-reads_1` would correspond to all reads with that specific prefix for the first demux set, for example:

```
s3://nao-testing/outward-assembly-test-data/siz/simulated-abcbd-reads_1_div0001.fastq.zst
s3://nao-testing/outward-assembly-test-data/siz/simulated-abcbd-reads_1_div0002.fastq.zst
s3://nao-testing/outward-assembly-test-data/siz/simulated-abcbd-reads_1_div0003.fastq.zst
```

## Working directory structure
```
.
├── current_contigs.fasta # working set of contigs. Used for read filtering when adapter trimming is disabled.
├── input_s3_paths.txt # list of input paths, copied for debugging
├── log.txt # collects some command output; not well structured
├── megahit_out_iter<i>-<j> # iteration i, subiteration j
│   ├── chose_this_subiter # empty file created if this subiter's contigs were chosen
│   ├── contigs_filtered.fasta # final.contigs.fa filtered via overlap graph logic
│   ├── final.contigs.fa # megahit output of this subiter's assembly
│   └── # other megahit outputs
├── original_seed.fasta # 
├── query_kmers.fasta # used for filtering all input reads, only appears if adapter trimming is enabled
├── reads # reads_* from each iter copied here for debugging
│   └── iter_<i> 
├── reads_1.fastq # reads used this iteration, will be copied to reads/
├── reads_2.fastq
├── reads_ff_1.fastq # _ff reads only appear if frequency filtering is enabled
├── reads_ff_2.fastq
├── reads_untrimmed_1.fastq # _untrimmed reads only appear if adapter trimming is eabled
└── reads_untrimmed_2.fastq
└── config.yaml # the configuration file used to run the pipeline
```
Note that kmer counting logic occurs in a separate `kmers` directory which is created by `_high_frequency_kmers`.

## Tips
### Choosing a good seed
The seed serves two purposes in outward assembly:
1. It's the initial contig, i.e. we use kmers from the seed to filter reads in the first iteration.
2. It's used to filter contigs at the end of each iteration.

At present, there's no error tolerance  in the latter filtering step: contigs must contain your seed exactly (up to reverse complementing). Therefore your seed really must be error-free. If your seed has an error (relative to the likely genome that generated the reads you're interested in), then the following sad sequence occurs:
1. Outward assemble uses the seed to find read pairs containing seed kmers.
2. These read pairs are assembled in the first iteration.
3. No contig output in the first iteration exactly contains the seed.
4. Therefore, the algorithm did not progress in the first iteration and terminates early.

Ideally, your seed is the minimal sequence such that you're interested in a contig if and only if the contig contains the seed. In practice, 25-50bp seems to work well. If you're trying to accelerate the outward assembly pipeline by providing a longer sequence with more kmers, use the `warm_start_path` argument rather than elongating your seed.

## (Optional) AWS Batch for read search

For very large datasets (potentially tens to hundreds of billions of reads), you can run the read search step on AWS Batch instead of locally. This trades startup time for potentially very wide parallelization—useful when local execution would be impractically slow.

This requires additional setup—see [installation docs](installation.md#optional-aws-batch-for-read-search).

To enable, call `outward_assembly` with `use_batch=True`, passing values for `batch_workdir`, `batch_queue`, and `tower_token`.
