# Algorithmic Details

While the basic logic of outward assembly -- filter reads, assemble, repeat -- is simple, we've introduced a few complications to avoid potential pitfalls and improve performance.

## Input data
The vast majority of outward assembly computational time is spent scanning reads to find those few reads with kmers matching the seed or contigs. We have a lot of reads! It's best to store these reads in S3 and stream them, rather than accessing them on local disk (i.e. an EBS volume on EC2):
* Even very expensive EBS volumes have limited throughput, whereas EC2 instances have comparatively high network bandwidth. Thus it's faster to stream large amounts of data from S3 than from a mounted EBS volume.
    * As of this writing in December 2024, EBS gp3 have default throughput of 125MB/s with up to 1GB/s at nontrivial cost. c7a EC2 instances start at ~1.5GB/s network throughput and large instances (which we use for big assembly problems) have around 6GB/s.
* We may wish to assemble from terabytes of reads at a time; an EBS volume that could hold all these would be expensive. Few instance types have large physically attached SSDs.
* Long term read storage is in S3, so streaming from S3 reduces data copying/movement.

Outward assembly only accepts input reads that are (a) stored in S3 (b) in SIZ format: **s**plit, **i**nterleaved, **z**std-compressed:
* For better parallelism, we split input reads into chunks of 1 million read pairs each.
* Read pairs are interleaved, rather than being stored in separate forward (`_1`) and reverse (`_2`) files. This is handy for streaming decompression and piping to the k-mer filtering tool.
* [zstd](https://github.com/facebook/zstd) is a lossless compression algorithm, similar in both concept and usage to gzip. Zstd has better compression ratio and decompression speed that gzip/pigz and is supported by Meta. Since we can stream decompress for outward assembly, we're not limited to the compression formats our tools support natively and can choose a more modern and efficient format.

## Overlapping contigs
Consider a genome with sections `ABCBD` with the seed sequence contained in the repeated `B` region. With sufficiently plentiful and long reads, an assembler might correctly recover the full genome, repeats and all. However, in practice -- and especially in outward assembly -- we'll run into ambiguities: after the contig has grown from the seed to the entire `B` region, it will be ambiguous if `B` should be preceeded by `A` or `C` and succeeded by `C` or `D`.

Typically assemblers (including megahit) deal with uncertainty -- branches in the assembly graph -- by forming new contigs for each branch [that can't be pruned]. In this case, as the assembler reaches assembly graph branch points at the end of the `B` region, it'll form new contigs that partially overlap `B`, like so:

![ABCBD contigs](../readme_images/abcbd_contigs.png)

Grey represents the unknown source genome, blue the seed sequence, and green the contigs. At this point, we've found the largest seed-containing contig we can, and the simple algorithm described up top terminates. But we shouldn't terminate -- the contigs overlapping `B` might provide important context about our seed sequence! The insight here is that a contig is important if it contains the seed _or_ if it _overlaps_ a contig that contains the seed.

So during the step "filter contigs to those that contain the seed," the full outward assembly algorithm:
* Creates an "overlap graph" where vertices are contigs and two vertices are connected by an edge if their respective contigs overlap.
* Partitions this overlap graph into connected components.
* Retains all contigs in a connected component with a contig that contains the seed sequence.

Ideally, this modified filtering step retains all contigs that could be stitched together with a seed-containing contig while discarding irrelevant unrelated contigs.

## Assembly strictness and coverage thresholds
Sequencing errors are common, and assemblers like megahit use coverage depth to decide when to apply error reduction procedures like tip pruning and bubble popping. Unfortunately, we're often assembling a sequence with very low coverage; maybe we've seen the seed just a handful of times -- if we've seen the seed a bunch and have deep coverage, we should have been trying to assemble one delivery ago.

Given the possibility of low coverage, we can't always set stringent coverage thresholds in assembly; that might prevent us from building longer contigs with the few (possibly correct) reads we have. We'd generally like to be as strict as possible about coverage while still managing to extend the contigs from last iteration.

Following this intuition, each iteration of outward assembly runs several assemblies, called subiters. These subiters are ordered such that `i-1` is the strictest assembly in iteration `i`, `i-2` the next strictest, etc. Each subiter produces a set of output contigs which we filter using the above graph procedure. The final iteration contigs come from the strictest subiter whose filtered contigs still progressed from last iteration; if there is no such subiter, we terminate early.

## High frequency kmers and chimeras
Technical chimeras are common sources of sequencing errors. When a technical chimera between `<thing you care about>` and `<something else>` arises, the `<something else>` is likely to be something very common like ribosomal RNA. This can lead to the following problem (using rRNA as an example):
* Iteration `i`: pull in a read covering a technical chimera where one half of the technical chimera is related to the seed and one half is rRNA.
* Iteration `i`: build a contig that contains the seed and some rRNA kmers.
* Iteration `i+1`: pull in all reads which contain kmers from iteration `i` contigs, so end up pulling in 5+% of all reads because rRNA is super common.
* Iteration `i+1`: with hundreds of millions of reads to assemble, assembly takes a prohibitively long time.

Outward assembly contains a simple escape hatch to prevent this situation: if any iteration pulls in more than `excess_read_thresh` reads, the algorithm exits early. But it's better if we never reach this situation in the first place! The key insight here is that frequenly occurring (like rRNA kmers) are highly unlikely to be part of the rare genome that contains the seed.

Outward assembly acccepts an optional path to a fasta file of high frequency kmers (or high frequency sequences of any length, though typically this file is generated from kmer counts). If provided, these high frequency kmers are used to generate a "frequency filtered" (`_ff`) set of reads that don't contain a high frequency kmer. These reads are used for one assembly each iteration. The assembly over frequency filtered reads is the first (most preferred) subiteration, i.e. outward assembly only uses non-frequency-filtered reads when assembly frequency-filtered reads fail to extend the previous iteration's contigs.

### Tradeoffs between when to filter
It's also possible to frequency filter reads up front and simply run the outward assembly algorithm over the subset of reads that don't contain high frequency kmers. A utility method `frequency_filter_reads` in `kmer_freq_filter.py` is provided. (Though as of this writing, not well tested or much used.)

Advantages:
* Most reads contain very common kmers, so even with a strict definition of "high frequency kmer," frequency filtering reads will often reduce the total data size at least by a factor of 2, meaning that the core iterations of outward assembly run much faster.
* No danger of mixing a high frequency kmer into your contigs since outward assembly never sees any high frequency kmers.

Disadvantages:
* Time/compute cost of up front filtering.
* If the true genome near the seed does contain a high frequency kmer, you can't recover it.
