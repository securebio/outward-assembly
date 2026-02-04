##
# Example script demonstrating the use of outward assembly.
# This script is callable as-is (though make sure your AWS credentials are configured).
##

import outward_assembly as oa_module
from outward_assembly.io_helpers import s3_files_with_prefix
from outward_assembly.pipeline import outward_assembly
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

##
# Set up paths
##

script_dir = Path(__file__).parent
oa_root = Path(oa_module.__file__).parent.parent

# Every invocation of `outward_assembly` will create a randomly-named work directory.
# The parent directory where these per-invocation workdirs will be created must exist
# before calling outward_assembly.
work_dir_parent = script_dir / "work"
work_dir_parent.mkdir(exist_ok=True)

##
# Collect read files we're going to search through
##

prefixes = [
    "outward-assembly-test-data/pipeline-end-to-end/simulated",
]
paths = [
    p for prefix in prefixes for p in s3_files_with_prefix("nao-testing", prefix)
]
print(f"Got {len(paths)} paths")

##
# Compute high frequency kmers
##

# In real usage, your SIZ read files are probably chunks of larger data. E.g. you might
# have 1000 SIZ chunks which jointly are the output of one lane on one flow cell.
# Assuming reads are well mixed across SIZ chunks, you can safely compute high frequency
# kmers in just a few of these chunks.

high_freq_path = work_dir_parent / "kmers/high_freq_kmers.fasta"
if not Path(high_freq_path).is_file():
    from outward_assembly.kmer_freq_filter import high_freq_kmers_split_files
    from outward_assembly.io_helpers import process_s3_paths

    # high_freq_kmers_split_files defaults to k=31; make sure the kmer size used here
    # is at least as large as the read_subset_k used in outward_assembly
    high_freq_path = high_freq_kmers_split_files(
        # Adjust the step size in the below range to control the fraction of read files we compute
        # high frequency kmers in. In this example script we want to hit every file, but in
        # production work with thousands of SIZ chunks, you might do range(0, len(paths), 500)
        process_s3_paths([paths[i] for i in range(0, len(paths), 1)]),
        work_dir_parent,
        min_kmer_freq=200,
        num_parallel=2,  # or more on a big machine with plenty of EBS space, roughly cores/4
    )
    print("Got high frequency kmers")

##
# Actually invoke OA!
##

outward_assembly(
    s3_paths=paths,
    seed_path=script_dir / "seed_seq.fasta", # can contain one or multiple seeds
    output_path=script_dir / "final_contigs.fasta", # will error if this already exists
    high_freq_kmers_path=high_freq_path, # optional
    adapters_path=oa_root / "default_adapters.fasta", # or bring your own adapters
    warm_start_path=script_dir / "warm_start.fasta", # optional
    work_dir_parent=work_dir_parent,
    read_subset_k=27,
    max_iters=3,
    excess_read_thresh=1_000_000,
)
