import json
import logging
import shutil
import tempfile
import time
from pathlib import Path
from typing import List, Optional, TypedDict

from Bio import SeqIO
from Bio.Seq import Seq

from .io_helpers import (
    PathLike,
    S3Files,
    _count_lines,
    process_s3_paths,
)
from .pipeline_steps import (
    CURRENT_CONTIGS,
    READS_1_FASTQ,
    _adapter_trim_iter_reads,
    _assemble_contigs,
    _choose_best_subiter,
    _copy_iteration_reads,
    _fasta_longest_total,
    _frequency_filter_reads,
    _prepare_query_seqs,
    _record_inputs,
    _subset_contigs,
    _subset_split_files,
)

logger = logging.getLogger(__name__)


class InnerIterationMetrics(TypedDict):
    """
    Store metrics for each iteration (referred to formally as an "inner iteration") of the assembly process.

    The read_pair_counts description will vary based on the user's input. By default, it is the number of read pairs that are extracted from BBDuk. In the situation that the user passed in an adapter path or high frequency kmers, this is the number of reads that were adapter trimmed or frequency filtered, respectively.
    """

    iteration: int
    read_pair_count: int
    contig_count: int
    longest_contig_length: int
    total_contig_length: int


class AssemblyMetrics(TypedDict):
    """
    Store statistics for the entire assembly process. The variables represent the final values after all iterations have been run. Inner iterations are stored in the inner_iterations variable, which contains a list of InnerIterationMetrics dictionaries for each iteration.
    """

    total_time: float
    inner_iterations: List[InnerIterationMetrics]
    final_contig_count: int
    final_longest_contig_length: int
    final_total_contig_length: int
    final_read_pair_count: int
    work_dir: str


def _compute_iteration_metrics(iter: int, workdir: Path) -> InnerIterationMetrics:
    """Collect metrics for the current iteration.

    Args:
        iter (int): The current iteration number.
        workdir (PathLike): The working directory containing the read and contig files.

    Returns:
        InnerIterationMetrics: A dictionary containing the metrics for the current iteration.
    """
    # Note that _adapter_trim_iter_reads puts trimmed reads in READS_{1/2}_FASTQ, so if the
    # user has adapter trimming enabled and we call this function after trimming adapters,
    # we'll be counting adapter-trimmed reads -- which is presumably what we want to count.
    read_pair_count = _count_lines(workdir / READS_1_FASTQ) // 4

    # Get contig stats
    contig_stats = _fasta_longest_total(workdir / CURRENT_CONTIGS)

    return {
        "iteration": iter,
        "read_pair_count": read_pair_count,
        "contig_count": sum(1 for _ in SeqIO.parse(workdir / CURRENT_CONTIGS, "fasta")),
        "longest_contig_length": contig_stats.longest,
        "total_contig_length": contig_stats.total,
    }


def _compute_assembly_metrics(
    workdir: PathLike,
    inner_iterations: List[InnerIterationMetrics],
    start_time: float,
) -> AssemblyMetrics:
    """Collect final assembly statistics after all iterations.

    Args:
        workdir (PathLike): The working directory containing the read and contig files.
        inner_iterations (List[InnerIterationMetrics]): List of metrics from each iteration.
        start_time (float): The start time of the assembly process.

    Returns:
        AssemblyMetrics: A dictionary containing the final statistics for the assembly process.
    """
    # Initialize all required fields, with defaults for the case of empty inner_iterations
    summary: AssemblyMetrics = {
        "total_time": time.time() - start_time,
        "inner_iterations": inner_iterations,
        "final_contig_count": 0,
        "final_longest_contig_length": 0,
        "final_total_contig_length": 0,
        "final_read_pair_count": 0,
        "work_dir": str(workdir),
    }

    # Update with actual values if we have iteration data
    if len(inner_iterations) > 0:
        for k in [
            "contig_count",
            "longest_contig_length",
            "total_contig_length",
            "read_pair_count",
        ]:
            summary["final_" + k] = inner_iterations[-1][k]

    return summary


def _outward_main_loop(
    s3_records: S3Files,
    workdir: PathLike,
    seed_seqs: List[Seq],
    overlap_contig_filtering: bool,
    read_subset_k: int,
    max_iters: int,
    high_freq_kmers_path: Optional[PathLike],
    freq_filter_k: int,
    adapters_path: Optional[PathLike],
    excess_read_thresh: int,
    warm_start_path: Optional[PathLike],
    use_batch: bool,
    batch_workdir: Optional[PathLike],
    batch_queue: Optional[str],
    tower_token: Optional[str] = None,
) -> AssemblyMetrics:
    """Main filter reads => assemble => filter contigs loop.

    Please see the docstring for outward_assembly for parameter descriptions.
    """
    workdir = Path(workdir)
    if not workdir.is_dir():
        raise ValueError(f"Working directory does not exist: {workdir}")
    if not (workdir / CURRENT_CONTIGS).is_file():
        raise ValueError(f"Current contigs file not found: {workdir / CURRENT_CONTIGS}")

    start_time = time.time()
    inner_iterations: List[InnerIterationMetrics] = []

    for iter in range(1, max_iters + 1):
        logger.debug(f"Iteration {iter}")

        # Step 1: subset reads to those containing kmers from current contigs

        kmer_ref_path = _prepare_query_seqs(
            workdir=workdir,
            warm_start_path=None if iter > 1 else warm_start_path,
            adapters_path=adapters_path,
            read_subset_k=read_subset_k,
        )

        _subset_split_files(
            s3_records,
            kmer_ref_path,
            read_subset_k=read_subset_k,
            workdir=workdir,
            use_batch=use_batch,
            batch_workdir=batch_workdir,
            batch_queue=batch_queue,
            tower_token=tower_token,
        )

        logger.debug("Subsetted reads")
        if (
            iter_read_count := _count_lines(workdir / READS_1_FASTQ) // 4
            > excess_read_thresh
        ):
            logger.debug(f"Too many reads ({iter_read_count}) at iter {iter}.")
            break

        # Step 2: prepare and assemble reads
        if adapters_path is not None:
            _adapter_trim_iter_reads(workdir, Path(adapters_path))
        if high_freq_kmers_path is not None:
            _frequency_filter_reads(workdir, Path(high_freq_kmers_path), freq_filter_k)
        _copy_iteration_reads(workdir, iter)
        _assemble_contigs(workdir, iter, high_freq_kmers_path is not None)
        logger.debug("Ran megahit")

        # Step 3: subset assembled contigs to those containing seed
        _subset_contigs(workdir, iter, seed_seqs, include_overlaps=overlap_contig_filtering)
        logger.debug("Found seed-containing contigs")

        # Step 4: check if we made progress this iteration and compute metrics.
        # We compute metrics after calling _choose_best_subiter, since
        # _choose_best_subiter updates the current contigs fasta, which is
        # read by _compute_iteration_metrics.
        progressed_this_iteration = _choose_best_subiter(workdir, iter)
        inner_iterations.append(_compute_iteration_metrics(iter, workdir))
        if progressed_this_iteration:
            if iter < max_iters:
                logger.debug(f"Progressed in iteration {iter}; continuing")
            else:
                logger.debug("Progressed but at max iterations")
        else:
            logger.debug("No contig improvement found; exiting")
            break

    # Collect final stats
    assembly_metrics = _compute_assembly_metrics(workdir, inner_iterations, start_time)

    # Save stats to JSON in work directory
    with open(workdir / "assembly_metrics.json", "w") as f:
        json.dump(assembly_metrics, f, indent=2)

    return assembly_metrics


def outward_assembly(
    s3_paths: List[str],
    seed_path: PathLike,
    output_path: PathLike,
    overlap_contig_filtering: bool = False,
    *,
    use_batch: bool = False,
    batch_workdir: Optional[PathLike] = None,
    batch_queue: Optional[str] = None,
    tower_token: Optional[str] = None,
    high_freq_kmers_path: Optional[PathLike] = None,
    freq_filter_k: int = 21,
    adapters_path: Optional[PathLike] = None,
    max_iters: int = 20,
    read_subset_k: int = 21,
    work_dir_parent: Optional[PathLike] = None,
    cleanup: bool = False,
    overwrite_output: bool = False,
    excess_read_thresh: int = 100_000,
    warm_start_path: Optional[PathLike] = None,
) -> AssemblyMetrics:
    """Assembly algorithm: iterative outward assembly from a seed.

    Uses BBDuk for read filtering and Megahit for assembly. Can optionally filter out
    reads containing high frequency kmers and/or adapturs during each iteration.

    Args:
        s3_paths: List of s3 paths to input read files
        seed_path: Fasta containing seed sequence
        output_path: Where to copy final contigs fasta
        overlap_contig_filtering: Whether to include contigs that do not contain but are connected to seed
        use_batch: Whether to use batch mode for outward assembly or local implementation
        batch_workdir: If using batch mode, you must provide an s3 path to a work directory for the batch job
        batch_queue: If using batch mode, you must provide a batch queue name
        tower_token: Seqera access token
        high_freq_kmers_path: Path to fasta file containing high frequency kmers to use
            for per-iteration frequency filtering, or None to disable per-iteration
            filtering. If you want to frequecy filter all your reads once up front, see
            the methods in kmer_freq_filter.py and then pass the filtered read files
            as S3 paths here (likely with high_freq_kmers_path=None).
        freq_filter_k: Kmer sized used for frequency filtering reads. Ignored if
            high_freq_kmers_path is None.
        adapters_path: Fasta file of adapters for per-iteration read trimming, or None
            to disable per-iteration trimming (e.g. if you started with cleaned reads)
        max_iters: Maximum assembly iterations before stopping
        read_subset_k: Kmer size for subsetting reads (default 21)
        work_dir_parent: Parent dir for temp files (default: home directory)
        cleanup: Whether to remove temporary files after completion (default False)
        overwrite_output: Whether to overwrite an existing output file (default False)
        excess_read_thresh: Exit early if we ever pull in this many reads in an
            iteration (default 100k)
        warm_start_path: Optional path to fasta file containing "warm start" sequences
            used to grab more reads in the first iteration of assembly. Prefer warm start
            over providing an over-long and potentially erroneous seed.
    Returns:
        AssemblyMetrics: A dictionary containing the summarized metrics for the assembly process.
    Raises:
        ValueError: If input directories/files don't exist or have unexpected structure
    """
    # Munge and validate inputs
    s3_records = process_s3_paths(s3_paths)
    seed_path = Path(seed_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not seed_path.is_file():
        raise ValueError(f"Seed sequence file not found: {seed_path}")
    if output_path.is_file() and not overwrite_output:
        raise ValueError(f"Output file already exists: {output_path}")
    if high_freq_kmers_path is not None:
        high_freq_kmers_path = Path(high_freq_kmers_path)
        if not high_freq_kmers_path.is_file():
            raise ValueError(f"High frequency kmers file not found: {high_freq_kmers_path}")

    # Set up temporary work directory
    if work_dir_parent is None:
        work_dir_parent = Path.home()
    else:
        work_dir_parent = Path(work_dir_parent)
        if not work_dir_parent.is_dir():
            raise ValueError(f"Parent directory does not exist: {work_dir_parent}")
    workdir = Path(tempfile.mkdtemp(dir=work_dir_parent))

    if use_batch and batch_workdir is None:
        raise ValueError("batch_workdir must be provided if use_batch is True")
    if use_batch and not str(batch_workdir).startswith("s3://"):
        raise ValueError("batch_workdir must be an s3 path")
    if use_batch and batch_queue is None:
        raise ValueError("batch_queue must be provided if use_batch is True")

    # Let's assemble
    try:
        logger.info(
            f"Outward assembly:\n  seed path: {seed_path}\n  workdir: {workdir}\n  output path: {output_path}"
        )
        _record_inputs(workdir, seed_path, s3_paths)

        # Copy seed sequences to start as current contigs
        current_contigs = workdir / CURRENT_CONTIGS
        shutil.copy2(seed_path, current_contigs)

        # Read all seed sequences from multi-fasta
        seed_seqs = [record.seq for record in SeqIO.parse(seed_path, "fasta")]

        if not seed_seqs:
            raise ValueError("No seed sequences found in seed fasta file")

        # Run the main assembly loop
        assembly_metrics = _outward_main_loop(
            s3_records=s3_records,
            workdir=workdir,
            seed_seqs=seed_seqs,
            overlap_contig_filtering=overlap_contig_filtering,
            read_subset_k=read_subset_k,
            max_iters=max_iters,
            high_freq_kmers_path=high_freq_kmers_path,
            freq_filter_k=freq_filter_k,
            adapters_path=adapters_path,
            excess_read_thresh=excess_read_thresh,
            warm_start_path=warm_start_path,
            use_batch=use_batch,
            batch_workdir=batch_workdir,
            batch_queue=batch_queue,
            tower_token=tower_token,
        )

        # Copy final contigs to output location

        shutil.copy2(current_contigs, output_path)

        return assembly_metrics

    finally:
        if cleanup:
            shutil.rmtree(workdir)
        else:
            logger.debug(f"Keeping working directory at {workdir}")
