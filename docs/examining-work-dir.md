# Examining an outward assembly working directory

Each invocation of outward assembly creates a working directory which records the state and progress of the algorithm at each iteration. Here's the structure at a high level:

```
.
├── current_contigs.fasta # working set of contigs. Used for read filtering when adapter trimming is disabled.
├── input_s3_paths.txt # list of input paths, copied for debugging
├── log.txt # collects some command output; not well structured
├── megahit_out_iter<i>-<j> # iteration i, subiteration j
│   ├── chose_this_subiter # empty file created if this subiter's contigs were chosen
│   ├── contigs_filtered.fasta # final.contigs.fa filtered to contigs relating to our seeds
│   ├── final.contigs.fa # megahit output of this subiter's assembly
│   └── # other megahit outputs
├── original_seed.fasta
├── query_kmers.fasta # used for filtering all input reads, only appears if adapter trimming is enabled
├── reads # reads_* from each iter copied here for debugging
│   └── iter_<i> 
├── reads_1.fastq # reads used this iteration, will be copied to reads/
├── reads_2.fastq
├── reads_ff_1.fastq # _ff reads only appear if frequency filtering is enabled
├── reads_ff_2.fastq
├── reads_untrimmed_1.fastq # _untrimmed reads only appear if adapter trimming is enabled
└── reads_untrimmed_2.fastq
└── config.yaml # the configuration file used to run the pipeline
```
(Note that kmer counting logic occurs in a separate `kmers` directory which is created by `high_freq_kmers_split_files`.)

The remainder of this document outlines how to use the working directory to answer specific questions about the algorithm run, using the [example script](example-assembly-dir/script.py)'s output as an example.

## Why did outward assembly stop, and which subiterations were chosen?

First, check how many iterations ran.
```
$ ls | grep "megahit_out"
megahit_out_iter1-1
megahit_out_iter1-2
megahit_out_iter1-3
megahit_out_iter2-1
megahit_out_iter2-2
megahit_out_iter2-3
megahit_out_iter3-1
megahit_out_iter3-2
megahit_out_iter3-3
```
Here we see outward assembly made it through 3 iterations (with subiterations `-1` through `-3` each). The example script has `max_iters=3`, so we didn't terminate early due to convergence or crash.

Then check to see if any subiteration in iteration 3 has a `chose_this_subiter` flag file. If yes, that means the contigs were extended in the last iteration, and we stopped because of `max_iters`. If no, that means the contigs weren't extended in the last iteration, and we stopped because of convergence.

```
$ ls megahit_out_iter3-*
megahit_out_iter3-1:
checkpoints.txt  contigs_filtered.fasta  done  final.contigs.fa  intermediate_contigs  log  options.json

megahit_out_iter3-2:
checkpoints.txt  contigs_filtered.fasta  done  final.contigs.fa  intermediate_contigs  log  options.json

megahit_out_iter3-3:
checkpoints.txt  contigs_filtered.fasta  done  final.contigs.fa  intermediate_contigs  log  options.json
Here we see that no subiteration was chosen in iteration 3, so we didn't extend contigs after iteration 2. The algorithm would have exited here even if we didn't have `max_iters=3`.
```

## How many reads did we find?

The reads found, filtered, and used for assembly in each iteration are stored in `reads/`. You can get a quick sense of how the algorithm progressed -- which iterations had lots of reads -- just by checking the size of these directories:

```
$ du -hd 1 reads/
312K	reads/iter_1
816K	reads/iter_2
1.2M	reads/iter_3
2.3M	reads/
```
Here we see linear-ish growth: 312K to 816K to 1.2M over the three iterations. On the other hand, if you saw a specific iteration where the size jumped dramatically, you should worry that iteration's read search found a bunch of (unrelated?) reads.

It's also helpful to check a specific iteration's reads to see exact read and base counts. For example, using `seqkit` (not an outward assembly dependency), we could check the read counts for iteration 2:
```
$ seqkit stats *
file                     format  type  num_seqs  sum_len  min_len  avg_len  max_len
reads_1.fastq            FASTQ   DNA        346   51,814      130    149.8      150
reads_2.fastq            FASTQ   DNA        346   51,848      140    149.8      150
reads_ff_1.fastq         FASTQ   DNA        346   51,814      130    149.8      150
reads_ff_2.fastq         FASTQ   DNA        346   51,848      140    149.8      150
reads_untrimmed_1.fastq  FASTQ   DNA        346   51,900      150      150      150
reads_untrimmed_2.fastq  FASTQ   DNA        346   51,900      150      150      150
```
Here `reads_` and `reads_ff` (`ff` for "[high] frequency filtered") have the same number of bases, so we see that frequency filtering didn't remove any reads this iteration. But both are a little shorter than `reads_untrimmed`, so we did strip a few bases with adapter trimming.

## What contigs did MEGAHIT build?

There are several reasons outward assembly fails to extend the previous iteration's contigs:
* MEGAHIT built more/longer contigs than in the previous iteration, but they no longer contain the seed sequence(s);
* MEGAHIT built the same contigs as it did last iteration;
* MEGAHIT built shorter contigs this iteration than last.

It's often helpful to examine the `megahit_out_iter<i>_<j>/final.contigs.fa` files to see MEGAHIT output per sub-iteration. In the case of the docs example script that didn't extend the contigs in the third iteration, we can first confirm that the filtered contigs didn't elongate from iteration 2 to iteration 3:
```
$ cat megahit_out_iter2-1/contigs_filtered.fasta
>k141_0 flag=1 multi=90.0000 len=1040
GGGGGGGGTTGGTAATGAAAAAAGATCTGTCGACGGAACACGAATTTGTTGAGGCCCGGG
ACTGGGAGTGTCGCATCGGATGGTGATCCATAAGTATAAGCTGCCGACGATGTTAGGGTC
GTGCGCACGTGATAGCTGTAGCTAAGCAAGACATGGACGCGCCCCATGCCAGCAGTTTTA
CTCCCACGGCGACGTCTGGACTACGCCTGTCTCTTTAGGAGATGCCCGATGCGACCTGCT
CCAACGACTCCCTCTGACGACTCCGTATTGAGAGGAAGGAGTCACCAATACATCTTAATA
GACAGCGCACTTATTAGTGCCCCGTCAGTGTAATCGACACTCAGCTTAGAGTTGTGTTAC
TAGTCTGATTCCATTAGGTAGCGTCGCATGTTTTGGTGCACCCAGAATTAATTGCCCGTC
TCGGTATGCAAAAATCTGTTCATCGAGCTGATCGCACTCTGTCGTTGCCAACTCAAATGC
TAGCTTAGTACCTATCCTCTTCGACGCATTGTTACGGGCCTAAAAGGTTCCAATCGAAAC
AGTCGATCGGATAGGGTACTGCGACATAAACAAAGCTCCCGTTGTCGACTTTGAGCACCG
GTAGTGAGTAATCCTTAGATCAAGCCTGGACCGGATGCATATCCTTGTGATAGACACTCA
GGGACAAGGTATCCCGACTTTATGCGCCTCGTAGGTGGAAAGTGTGCGCCCCTTGGCGTA
GCTTCTCCTACCACAAAACCGCAGCACCACAACACGCCTTCACCGGAATCCGGCTATGTC
AGTCGGAGGTCCAAGAGGCCTCTATATTAGCTGGGGGCTCGACTAGACGGCAAACCAGTT
ATTTTGCACATTAACAACTGTCGAGTACCCAGCTCCTTCACTGTTCTCTAACGGATCAGG
AATATAATAAGATCTACTGCGCCCCAGCCGCAGAACGATAAGGGCTGCATTTGCCTACCT
CTACAGCCATGTCGTACGGGAATCCAGCCTGCGGGGGGGGCCCCCCCCGCAGGCTGGATT
CCCGTACGACATGGCTGTAG
$ cat megahit_out_iter3-1/contigs_filtered.fasta
>k141_3 flag=0 multi=108.0000 len=1000
GGGGGGGGTTGGTAATGAAAAAAGATCTGTCGACGGAACACGAATTTGTTGAGGCCCGGG
ACTGGGAGTGTCGCATCGGATGGTGATCCATAAGTATAAGCTGCCGACGATGTTAGGGTC
GTGCGCACGTGATAGCTGTAGCTAAGCAAGACATGGACGCGCCCCATGCCAGCAGTTTTA
CTCCCACGGCGACGTCTGGACTACGCCTGTCTCTTTAGGAGATGCCCGATGCGACCTGCT
CCAACGACTCCCTCTGACGACTCCGTATTGAGAGGAAGGAGTCACCAATACATCTTAATA
GACAGCGCACTTATTAGTGCCCCGTCAGTGTAATCGACACTCAGCTTAGAGTTGTGTTAC
TAGTCTGATTCCATTAGGTAGCGTCGCATGTTTTGGTGCACCCAGAATTAATTGCCCGTC
TCGGTATGCAAAAATCTGTTCATCGAGCTGATCGCACTCTGTCGTTGCCAACTCAAATGC
TAGCTTAGTACCTATCCTCTTCGACGCATTGTTACGGGCCTAAAAGGTTCCAATCGAAAC
AGTCGATCGGATAGGGTACTGCGACATAAACAAAGCTCCCGTTGTCGACTTTGAGCACCG
GTAGTGAGTAATCCTTAGATCAAGCCTGGACCGGATGCATATCCTTGTGATAGACACTCA
GGGACAAGGTATCCCGACTTTATGCGCCTCGTAGGTGGAAAGTGTGCGCCCCTTGGCGTA
GCTTCTCCTACCACAAAACCGCAGCACCACAACACGCCTTCACCGGAATCCGGCTATGTC
AGTCGGAGGTCCAAGAGGCCTCTATATTAGCTGGGGGCTCGACTAGACGGCAAACCAGTT
ATTTTGCACATTAACAACTGTCGAGTACCCAGCTCCTTCACTGTTCTCTAACGGATCAGG
AATATAATAAGATCTACTGCGCCCCAGCCGCAGAACGATAAGGGCTGCATTTGCCTACCT
CTACAGCCATGTCGTACGGGAATCCAGCCTGCGGGGGGGG
```
Notice how the contig shrinks from 1040nt to 1000nt; the terminal `CCCCCCCCGCAGGCTGGATTCCCGTACGACATGGCTGTAG` is lost. Okay, why did this happen? Let's look at MEGAHIT's output in iteration 3-1, before any contig filtering:
```
$ cat megahit_out_iter3-1/final.contigs.fa
>k141_4 flag=0 multi=16.6643 len=421
GCCCACAAGCGTTGTTACGCATGATCTCGCTACATCCGATATTCGCGGCTACTGTGCTAAGCTGGAAATCCCCACGGCTTGGTCAAACGCCAGATGTGCCGCATGGTGCGCAGTGCGTATGCGTACGTCTGGGGCATGAGCACTAAGTCGCGAATACTCACGGAAGGGTTCCGGTGCACTGATTTAGGGGCAGCCAATGGGGCCGAGCTGACTATCGACAATGATGGGCGGAGGATCGCTTACGTCTATGCCAGCCTCGGGCTGTGCGAAGTGGGGGGGGCCCCCCCCGCAGGCTGGATTCCCGTACGACATGGCTGTAGAGGTAGGCAAATGCAGCCCTTATCGTTCTGCGGCTGGGGCGCAGTAGATCTTATTATATTCCTGATCCGTTAGAGAACAGTGAAGGAGCTGGGTACTCGAC
>k141_0 flag=0 multi=16.6148 len=411
GTCGAGTACCCAGCTCCTTCACTGTTCTCTAACGGATCAGGAATATAATAAGATCTACTGCGCCCCAGCCGCAGAACGATAAGGGCTGCATTTGCCTACCTCTACAGCCATGTCGTACGGGAATCCAGCCTGCGGGGGGGGTTTTTTTTGTTTGAGTTAGCCTGCGCTATCTACTGTTGTAGGCTAACTTTAGCTGTCACCAAACGCGAATACCCAACCTACGGCCTTAAATCTTATTTATTAGGAAAAGTTTTAATTTTGGTGAATCAACCCCATTAAAACGTACAGCCACTTGACCAACAACCACTCACCCTGTTGCTACGAGAAACCGCTCATGCTCACCGTTGAGTATCCGGGCACAGGTCTACTGCGAATTATTAGCAAGACTGTTCTGGTGCAAAGTTCGGGAAA
>k141_2 flag=0 multi=16.6148 len=411
CTACAGCTATCACGTGCGCACGACCCTAACATCGTCGGCAGCTTATACTTATGGATCACCATCCGATGCGACACTCCCAGTCCCGGGCCTCAACAAATTCGTGTTCCGTCGACAGATCTTTTTTCATTACCAACCCCCCCCTTTTTTTTGTTGCTATATCACTTAAGGTAGGCTATATATGTACATAAGCTGAGTCACTATAATCATCCAGATTTAACAACCCGGACCAAGGTCCCATTGGACGAAGACAGTAATAATACCGCTATTCTCAGTGATCTCAAAAGATGCTCTCAGGGTATGTAGTAGCTCCGCGCTTTATTCGGGCGGACATAATCGACAACCAGCTTTGTGGATTTGTGTATATCCCAAAGTATATTTCGTCACGCATTTATTGATTCTACGTAAAGACCT
>k141_1 flag=0 multi=15.6963 len=411
CTACAGCTATCACGTGCGCACGACCCTAACATCGTCGGCAGCTTATACTTATGGATCACCATCCGATGCGACACTCCCAGTCCCGGGCCTCAACAAATTCGTGTTCCGTCGACAGATCTTTTTTCATTACCAACCCCCCCCGGGGGGGGTAGCCGGAAAAGTCCATACCGCTATATGAGGCATATCAGCGCATTCCTGCGTTGAGGTACTAGACCGTGGCATCTGATAGTAGGGAGGTTTTGTCTTCACCTTCAAACGACTGCTCATGATATTAGGTTAATAGGTTTGTGGACGATATATCTGTTTCGTGCCTGCACGGTAATGGAGACCGTCTCCGCCACTTGCGCAAGCTTGTGTGCCGACACGAAGGCAGCTGGTAAAAGTGTGCCGCTCCTGCCCATTTGCCCCATC
>k141_3 flag=0 multi=108.0000 len=1000
CCCCCCCCGCAGGCTGGATTCCCGTACGACATGGCTGTAGAGGTAGGCAAATGCAGCCCTTATCGTTCTGCGGCTGGGGCGCAGTAGATCTTATTATATTCCTGATCCGTTAGAGAACAGTGAAGGAGCTGGGTACTCGACAGTTGTTAATGTGCAAAATAACTGGTTTGCCGTCTAGTCGAGCCCCCAGCTAATATAGAGGCCTCTTGGACCTCCGACTGACATAGCCGGATTCCGGTGAAGGCGTGTTGTGGTGCTGCGGTTTTGTGGTAGGAGAAGCTACGCCAAGGGGCGCACACTTTCCACCTACGAGGCGCATAAAGTCGGGATACCTTGTCCCTGAGTGTCTATCACAAGGATATGCATCCGGTCCAGGCTTGATCTAAGGATTACTCACTACCGGTGCTCAAAGTCGACAACGGGAGCTTTGTTTATGTCGCAGTACCCTATCCGATCGACTGTTTCGATTGGAACCTTTTAGGCCCGTAACAATGCGTCGAAGAGGATAGGTACTAAGCTAGCATTTGAGTTGGCAACGACAGAGTGCGATCAGCTCGATGAACAGATTTTTGCATACCGAGACGGGCAATTAATTCTGGGTGCACCAAAACATGCGACGCTACCTAATGGAATCAGACTAGTAACACAACTCTAAGCTGAGTGTCGATTACACTGACGGGGCACTAATAAGTGCGCTGTCTATTAAGATGTATTGGTGACTCCTTCCTCTCAATACGGAGTCGTCAGAGGGAGTCGTTGGAGCAGGTCGCATCGGGCATCTCCTAAAGAGACAGGCGTAGTCCAGACGTCGCCGTGGGAGTAAAACTGCTGGCATGGGGCGCGTCCATGTCTTGCTTAGCTACAGCTATCACGTGCGCACGACCCTAACATCGTCGGCAGCTTATACTTATGGATCACCATCCGATGCGACACTCCCAGTCCCGGGCCTCAACAAATTCGTGTTCCGTCGACAGATCTTTTTTCATTACCAACCCCCCCC
```
Notice how that the `CCCCCCCCGCAGGCTGGATTCCCGTACGACATGGCTGTAG` "missing" from our iteration 3 contig is found in MEGAHIT output contig `k141_4`! With the extra reads found in iteration 3, MEGAHIT's assembly took a different path through the `CCCCCCCCGCAGGCTGGATTCCCGTACGACATGGCTGTAG` assembly graph node than we had in iteration 2.
