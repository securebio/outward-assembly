# Nucleaze vs BBDuk Benchmark Results

Benchmarks run on macOS ARM64 (M1 Pro, 16GB RAM)
Both tools configured identically: k=23, canonical k-mers, exact matching, ordered output, 4 threads.

## K-mer Filter Scaling

Input: simulated paired-end reads (150bp) duplicated at various multipliers.
Reference: 10 k-mers from `tests/data/pipeline-end-to-end/fake_hfk.fasta`.

| Read pairs | Nucleaze (s) | BBDuk (s) | Speedup | Nucleaze peak mem | BBDuk peak mem | Mem ratio |
|-----------:|-------------:|----------:|--------:|------------------:|---------------:|----------:|
| 941 | 0.005 | 0.151 | 28x | 4 MB | 79 MB | 19x |
| 9,410 | 0.011 | 0.198 | 17x | 12 MB | 90 MB | 8x |
| 47,050 | 0.083 | 0.245 | 3x | 38 MB | 132 MB | 3x |
| 94,100 | 0.081 | 0.346 | 4x | 69 MB | 200 MB | 3x |
| 470,500 | 0.563 | 0.665 | 1.2x | 293 MB | 558 MB | 2x |
| 941,000 | 1.480 | 1.150 | ~1x | 544 MB | 993 MB | 2x |

Output hashes match between BBDuk and Nucleaze at all sizes tested.

## Full Pipeline

End-to-end pipeline run with 941 read pairs, 2 iterations, k-mer frequency filtering enabled.

| Branch | Scenario | Time (s) | Output hash |
|--------|----------|----------|-------------|
| main (BBDuk) | basic | 7.67 | 1bf62f2ba03350f5eeab17ebb745f3e6 |
| main (BBDuk) | multi_seed | 7.43 | 9f42727c290c073c476ccf289a63371b |
| nucleaze-migration | basic | 7.00 | 1bf62f2ba03350f5eeab17ebb745f3e6 |
| nucleaze-migration | multi_seed | 6.93 | 9f42727c290c073c476ccf289a63371b |

Output hashes match between branches
(the test dataset is small, so megahit assembly takes up most of the runtime)

## Reproducing

```bash
source /Users/lee/miniforge3/etc/profile.d/conda.sh && conda activate oa-tools
export PATH="$HOME/.cargo/bin:$PATH"

# K-mer filter scaling
uv run python benchmarks/kmer_filter_scaling.py

# Full pipeline
uv run python benchmarks/run_benchmark.py
```
