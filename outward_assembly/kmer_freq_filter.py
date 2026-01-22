import logging
import math
import os
import shutil
import subprocess
import warnings
from multiprocessing import cpu_count
from pathlib import Path
from typing import Optional

from .io_helpers import PathLike, S3Files, s3_stream_cmd

# kmc and kmc_tools are expected to be available on PATH via the oa-tools conda environment


def _make_kmer_count_commands(
    s3_records: S3Files,
    kmers_dir: str | Path,
    *,
    k: int,
    min_kmer_freq: int,
    memory_GB: int = 5,
    threads: int = 4,
) -> list[str]:
    """Helper function to construct shell commands for counting kmers in each input file.

    Args:
        s3_records: Files to count kmers in
        kmers_dir: Directory for kmer counting output
        k: Kmer size
        min_kmer_freq: Minimum frequency threshold
        memory_GB: Memory limit for KMC in GB
        threads: Number of threads per KMC process

    Returns:
        List of shell commands as strings
    """
    # Find next power of 2 above 2 * min_kmer_freq using logarithms
    max_count = 2 ** math.ceil(math.log2(2 * min_kmer_freq))

    kmers_dir = Path(kmers_dir)
    commands = []

    for rec in s3_records:
        # Define temporary and output paths
        tmp_dir = kmers_dir / f"tmp_{rec.filename}"
        tmp_fastq = tmp_dir / "reads.fastq"
        kmc_tmp_dir = tmp_dir / "kmc_tmp"
        kmc_out_prefix = tmp_dir / "kmers"
        high_freq_out = kmers_dir / f"hfq_{rec.filename}.txt"

        # Build command sequence
        cmd = (
            f"mkdir -p {tmp_dir} {kmc_tmp_dir} && "
            f"{s3_stream_cmd(rec.s3_path)} | "
            f"zstd -d - -o {tmp_fastq} && "
            f"kmc -k{k} "  # k-mer length
            f"-cs{max_count} "  # max count before counter saturates
            f"-ci{min_kmer_freq} "  # min count to store k-mer
            f"-t{threads} "  # number of threads
            f"-m{memory_GB} "  # max memory usage in GB
            f"{tmp_fastq} {kmc_out_prefix} {kmc_tmp_dir} && "
            f"kmc_tools transform {kmc_out_prefix} dump {high_freq_out} && "
            f"rm -rf {tmp_dir}"
        )
        commands.append(cmd)

    return commands


def _high_freq_kmers_split_files(
    s3_records: S3Files,
    workdir: PathLike,
    *,
    num_parallel: Optional[int] = None,
    min_kmer_freq: int = 2000,
    k: int = 31,
    allow_tmp_workdir: bool = True,
) -> Path:
    """Get list of kmers that appear at least min_kmer_freq times in any single input file.

    Uses KMC3 for kmer counting and operates in parallel while managing disk space usage.
    All work is done in workdir/kmers/. Returns the path to a fasta of high frequency kmers.

    Args:
        s3_records: Files to count kmers in
        workdir: Working directory for temporary files
        num_parallel: Number of parallel processes (defaults to max(1, cpu_count()//4))
        min_kmer_freq: Minimum frequency threshold for high-frequency kmers
        k: Kmer size
        allow_tmp_workdir: allow a /tmp backed working directory (with warning)

    Returns:
        Path to fasta file containing high-frequency kmers

    Raises:
        ValueError: If workdir starts with /tmp or doesn't exist
    """
    workdir = Path(workdir)

    # Input validation
    if not workdir.exists():
        raise ValueError(f"Working directory {workdir} does not exist")
    if str(workdir).startswith("/tmp"):
        if not allow_tmp_workdir:
            raise ValueError(
                f"Working directory {workdir} starts with /tmp which may be memory-backed"
            )
        else:
            warnings.warn(
                "workdir is in /tmp which is typically memory-backed. "
                "Large read/kmer files here could consume system memory and potentially cause crashes. "
                "If this is not intentional, consider using a disk-backed location instead."
            )

    # Set default num_parallel if not provided
    if num_parallel is None:
        num_parallel = max(1, cpu_count() // 4)

    # Create kmers subdir - error if it exists as it may contain stale data
    kmers_dir = workdir / "kmers"
    kmers_dir.mkdir()

    # Create and execute kmer counting commands
    cmds = _make_kmer_count_commands(
        s3_records, kmers_dir, k=k, min_kmer_freq=min_kmer_freq
    )

    cmd_file = kmers_dir / "kmer_count_commands.txt"
    with open(cmd_file, "w") as f:
        f.write("\n".join(cmds))

    # Run commands in parallel with xargs
    logging.debug(f"Running KMC commands in parallel with {num_parallel} processes")
    subprocess.run(
        f"cat {cmd_file} | xargs -P {num_parallel} -I CMD bash -c 'CMD'",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Merge results
    result_path = kmers_dir / "high_freq_kmers.fasta"
    high_freq_files = [
        f for f in os.listdir(kmers_dir) if f.startswith("hfq_") and f.endswith(".txt")
    ]

    if high_freq_files:
        # Process high frequency kmers
        kmer_counts = {}
        for hff in high_freq_files:
            with open(kmers_dir / hff) as f:
                for line in f:
                    # Each line is <kmer> <whitespace> <count>
                    # e.g. "ACGTACGT 2341"
                    kmer, count = line.strip().split()
                    # Store max count seen for this kmer across all input files
                    kmer_counts[kmer] = max(kmer_counts.get(kmer, 0), int(count))

        # Write merged results
        with open(result_path, "w") as f:
            for kmer, count in kmer_counts.items():
                f.write(f">kmer_max_{count}_obs\n{kmer}\n")

        # Clean up individual kmer files
        for f in high_freq_files:
            (kmers_dir / f).unlink()
    else:
        logging.warning("No high-frequency kmers found in any input file")
        result_path.touch()

    # Clean up command file
    cmd_file.unlink()

    return result_path


def frequency_filter_reads(
    s3_records: S3Files,
    workdir: PathLike,
    out_dir: PathLike,
    *,
    num_parallel: Optional[int] = None,
    min_kmer_freq: int = 2000,
    k: int = 31,
) -> None:
    """Filter out reads containing high-frequency kmers from split interleaved reads.

    Takes zstd-compressed fastq files named reads_il_divXXXX.fastq.zst from reads_dir,
    identifies kmers appearing >= min_kmer_freq times in any single file, and outputs
    filtered reads (excluding those with high-freq kmers) to out_dir with same naming scheme.

    Args:
        s3_records: Files to count kmers in
        workdir: Working directory for temporary files
        out_dir: Output directory (presently may not be S3)
        num_parallel: Number of parallel processes (defaults to max(1, cpu_count()//4))
            Note that disk space required is on the order of <input file size> * num_parallel
        min_kmer_freq: Minimum frequency threshold for high-frequency kmers
        k: Kmer size

    To-do:
        Allow outputting to S3
    """
    workdir = Path(workdir)
    out_dir = Path(out_dir)
    if out_dir.exists():
        logging.warning(
            f"Output directory {out_dir} already exists - skipping frequency filtering"
        )
        return
    out_dir.mkdir(parents=True, exist_ok=False)

    # Set default num_parallel if not provided
    if num_parallel is None:
        num_parallel = max(1, cpu_count() // 4)

    # Get high-frequency kmers
    high_freq_kmers_path = _high_freq_kmers_split_files(
        s3_records, workdir, num_parallel=num_parallel, min_kmer_freq=min_kmer_freq, k=k
    )

    out_paths = [out_dir / Path(rec.s3_path).parts[-1] for rec in s3_records]

    # Find nucleaze binary (needed for xargs/sh which may not inherit PATH)
    nucleaze_bin = shutil.which("nucleaze")
    if nucleaze_bin is None:
        raise RuntimeError("nucleaze not found in PATH")

    # Create Nucleaze+compression commands
    # --outu outputs reads that do NOT match (i.e., reads without high-freq kmers)
    cmds = [
        f"{s3_stream_cmd(rec.s3_path)} | zstdcat - | "
        f"{nucleaze_bin} --in - "
        f"--outu {out_dir / rec.filename}_tmp_out.fq "
        f"--outm /dev/null "
        f"--ref {high_freq_kmers_path} --k {k} "
        f"--canonical --minhits 1 --interinput "
        f"--threads 3 && "
        f"zstd -q -T3 < {out_dir / rec.filename}_tmp_out.fq > {p_out} && "
        f"rm {out_dir / rec.filename}_tmp_out.fq"
        for rec, p_out in zip(s3_records, out_paths)
    ]

    # Run commands in parallel using concurrent.futures (avoids xargs command length limits)
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def run_cmd(cmd: str) -> subprocess.CompletedProcess:
        return subprocess.run(cmd, shell=True, check=True, capture_output=True)

    with ThreadPoolExecutor(max_workers=num_parallel) as executor:
        futures = [executor.submit(run_cmd, cmd) for cmd in cmds]
        for future in as_completed(futures):
            future.result()  # Raises exception if command failed

    logging.debug("Nucleaze filtering complete")
