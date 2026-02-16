# Usage

This document covers how to run outward assembly.

Table of Contents:
1. [Overview](#overview)
2. [Configuring read input](#configuring-read-input)
3. [Working directory structure](#working-directory-structure)
4. [After outward assembly runs](#after-outward-assembly-runs)
5. [Tips](#tips)
6. [(Optional) AWS Batch for read search](#optional-aws-batch-for-read-search)

> **Note:** For experimental features including pileup visualizations and automated/iterative mode, see [experimental features](experimental.md).

## Overview

The primary entrypoint to outward assembly is the Python function `outward_assembly` in `outward_assembly/pipeline.py`; a command line interface does not (yet) exist. See the docstring of `outward_assembly` for a detailed description of keyword parameters. The required parameters are:
* The path to a seed sequence (fasta) to assemble outward from;
* A list of S3 paths of reads to assemble. Reads must be in [SIZ format](./algorithm_details.md#input-data);
* Path for output contigs.

Typically, you'll format your inputs (seed sequences, warm start sequences, adapter sequences, etc.) as fasta files and stick them into a directory along with a short Python script that calls `outward_assembly`. The docs contain a thoroughly-commented runnable [example](example-assembly-dir/script.py) of this workflow.

By default, the entire pipeline runs on your local machine. This works well for most use cases; if you need more resources, try scaling up your machine (e.g. pick an EC2 instance with more cores and memory). For very large datasets, there's an option to [run the read search step on AWS Batch](#optional-aws-batch-for-read-search), trading startup time for potentially very wide parallelization.

Outward assembly generates lots of intermediate results in its working directory; these are especially useful if you want to monitor the progress of an assembly in progress, retrieve results from intermediate iterations, or debug a crash. See the `outward_assembly` arguments `work_dir_parent` and `cleanup`.

## Configuring read input

Outward assembly requires streaming split, interleaved, zstd-compressed reads from S3; reads in other formats or locations are not (yet) supported. A script that calls `outward_assembly` typically begins by specifying the S3 paths of the read files to assemble; see `s3_files_with_prefix` in [io_helpers.py](../outward_assembly/io_helpers.py). Typical usage is to define a few prefixes (relative to bucket root) and collect all paths matching these prefixes:
```python
prefixes = [
	"deliveryA/siz/",
	"deliveryB/siz/"
]
paths = [p for prefix in prefixes for p in s3_files_with_prefix("my-bucket", prefix)]
```

## After outward assembly runs

You can see final output in the final output file you specified when calling `outward_assembly`. However, it's often useful to look around in the working directory to see how the algorithm progressed. A few things to check:
* How many iterations and subiterations were there?
* Which subiteration was chosen in each iteration? Look for empty `chose_this_subiter` sentinel files, which will appear in one sub-iteration per iteration.
	* It's common for the first few iterations to be able to use subiteration `-1`, but later iterations to need subiter `-2` or `-3`; this can correlate with a decline in contig quality. E.g. if you ran 10 iterations total, the first 6 of which used subiter `i-1`, compare the iteration 6 and iteration 10 (final) outputs.
* How many reads were selected in each iteration? Did any iterations pull in lots of high frequency kmer reads? Checking file sizes in the `reads/` dir is helpful.

See [Examining a working directory](docs/examining-work-dir.md) for a walkthrough of how to read the working directory of a completed outward assembly run.

## Tips

### Choosing a good seed

The seed serves two purposes in outward assembly:
1. It's the initial contig, i.e. we use kmers from the seed to filter reads in the first iteration.
2. It's used to filter contigs at the end of each iteration.

There's no error tolerance in the latter filtering step: contigs must contain your seed exactly (up to reverse complementing). Therefore your seed really must be error-free. If your seed has an error (relative to the likely genome that generated the reads you're interested in), then the following sad sequence occurs:
1. Outward assemble uses the seed to find read pairs containing seed kmers.
2. These read pairs are assembled in the first iteration.
3. No contig output in the first iteration exactly contains the seed.
4. Therefore, the algorithm did not progress in the first iteration and terminates early.

Your seed should be the minimal sequence such that you're interested in a contig if and only if the contig contains the seed. In practice, 24-50bp seems to work well.

### Multiple seeds

Outward assembly supports running with multiple seed sequences: just stick multiple sequences in your seed sequence fasta. The algorithm and logic are essentially unchanged from the one-seed case, except that at the end of each iteration, we retain any contig which fully contains _any_ of the seed sequences.

You might wish to run with multiple seeds if:
* You have many seed sequences that might come from the same genome, and you hope outward assembly will join them into a coherent contig.
* You're unsure which seed is correct (e.g. due to reads disagreeing about specific bases), so you assemble all plausible seeds.
* You've identified many seed seeds that should be assembled using the same underlying read data, and you want to avoid the computational cost of many independent outward assembly runs.
	* Note that running outward assembly with multiple seeds is not identical to running multiple single-seed outward assemblies. For example, including seed A might cause the read search step to pull in a specific read that affects the assembly of contigs containing seed B.

### Warm start

It's common to have sequences you think will be part of your final contigs -- perhaps up to some small errors -- but are not themselves the seed sequence. For example, your seed sequence might be a certain kmer of interest observed in several reads; the reads are likely to appear in the final contig, but aren't themselves the essential seed.

In these cases, you can use a _warm start_: point `warm_start_path` at these sequence(s) in fasta format, and all the kmers contained in these sequences will be used in the first iteration read search. This is different from extending the seed, because these seed sequences aren't used for contig filtering, so it's okay for them to contain errors. For example, suppose you have multiple reads that attest to a chimeric junction, but the reads disagree on bases near read ends or away from the jointly-attested junction. Then pick ~28bp around the junction as your seed, and include all the relevant reads in your warm start fasta.

There are no specific requirements for the warm start sequences: they can contain the seed, be adjacent to the seed, etc.; they're just a basket of kmers used to help find more reads in the first iteration read search, and they're only used in the first iteration.

### Choosing a good `read_subset_k`

The `read_subset_k` argument to `outward_assembly` is a key parameter governing the pipeline's behavior. A read in the corpus of input reads is "pulled in" for assembly if and only if it has an exact kmer match to the current contigs (or warm start, in iteration 1). Thus `read_subset_k` governs a sensitivity-specificity tradeoff:
* small `k` makes outward assembly more tolerant of sequencing errors and helps detect reads which overlap a contig by only a few bases.
* large `k` makes outward assembly less likely to pull in reads that aren't actually related to your seed sequence or unknown target genome.

In practice, start with a `read_subset_k` around 26. If you're failing to find good contigs, try a smaller `k`; if you're assembling lots of contigs that have nothing to do with your seed, try a larger `k`.

Make sure you either have a reliable warm start sequence or that `read_subset_k` is no larger than your seed. (If `read_subset_k > len(seed)`, there are no k-mers in your seed to search for!)

### What hardware to run on

Almost always, outward assembly runing time is dominated by the parallel BBDuk read searches: looking through a large haystack of reads to find the few needle-reads that contain kmers from our contigs. To make this search as fast as possible, you really want to be running on an EC2 instance in the same region as your data-containing S3 buckets. (This will also minimize data movement costs, since data S3 -> EC2 within region is free.)

The read search is compute and network bottlenecked, so consider compute-optimized instances like the c7a/c8a families. Outward assembly will run one BBDuk search process per 4 vCPU cores. *Very* roughly, streaming, decompressing, and searching a 1 million read pair SIZ chunk takes 4 cores about 15 seconds, so with a `n` core machine you can search `n`-million read pairs per minute.

## (Optional) AWS Batch for read search

For very large datasets (potentially tens to hundreds of billions of reads), you can run the read search step on AWS Batch instead of locally. This trades task startup time for potentially very wide parallelization -- useful when local execution would be impractically slow. In general, it's more efficient to run outward assembly on a larger machine than to use Batch, but Batch can be helpful for running thousands of read searches in parallel or using cheap spot instances.

This requires additional setup—see [installation docs](installation.md#optional-aws-batch-for-read-search).

To enable, call `outward_assembly` with `use_batch=True`, passing values for `batch_workdir`, `batch_queue`, and `tower_token`.
