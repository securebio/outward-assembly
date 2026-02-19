import os
import shutil
import subprocess
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional

import ahocorasick
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .basic_seq_operations import SeqOrientation
from .io_helpers import PathLike, S3Files, concat_and_tag_fastq, s3_stream_cmd
from .overlap_graph import get_overlapping_sequence_ids

# File names used across functions
CURRENT_CONTIGS = "current_contigs.fasta"
WARM_START = "seed_with_warm_start.fasta"
FINAL_CONTIGS = "contigs.fasta"
CONTIG_KMERS = "contig_kmers.fasta"
CANDIDATE_KMERS = "candidate_query_kmers.fasta"
QUERY_KMERS = "query_kmers.fasta"
READS_PREFIX = "reads_il"
READS_1_FASTQ = "reads_1.fastq"
READS_2_FASTQ = "reads_2.fastq"
READS_FILTERED_1_FASTQ = "reads_ff_1.fastq"
READS_FILTERED_2_FASTQ = "reads_ff_2.fastq"
READS_UNTRIMMED_1_FASTQ = "reads_untrimmed_1.fastq"
READS_UNTRIMMED_2_FASTQ = "reads_untrimmed_2.fastq"
MEGAHIT_OUT_PREFIX = "megahit_out_iter"
MEGAHIT_FINAL_CONTIGS = "final.contigs.fa"
MEGAHIT_FILTERED_CONTIGS = "contigs_filtered.fasta"
CHOSEN_SUBITER_FLAG = "chose_this_subiter"
LOG_FILE = "log.txt"

NF_PROFILE_ENV_VAR = "NEXTFLOW_PROFILE"


class FastaStats(NamedTuple):
    longest: int
    total: int


def _reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def _assemble_contigs(workdir: PathLike, iter: int, freq_filter: bool) -> None:
    """Assemble iteration contigs in working directory.

    Runs one or more subiteration assemblies with lower subiters corresponding to
    stricter or generally preferred assembly parameters. Generate subdirs
    workdir/megahit_out_iter<iter>-<subiter> for each assembly.

    Args:
        workdir: Working directory containing read files
        iter: Current iteration number
        freq_filter: Whether to run an assembly with frequency-filtered reads
    """
    workdir = Path(workdir)
    cmds = []
    subiter = 1

    if freq_filter:
        # First assembly using filtered reads
        # fmt: off
        cmd = [
            "megahit",
            "-1", str(workdir / READS_FILTERED_1_FASTQ),
            "-2", str(workdir / READS_FILTERED_2_FASTQ),
            "-o", str(workdir / f"{MEGAHIT_OUT_PREFIX}{iter}-{subiter}"),
            "--k-list", "15,21,29,39,59,79,99,119,141",
        ]
        # fmt: on
        cmds.append(cmd)
        subiter += 1

    # Standard assembly
    # fmt: off
    cmd = [
        "megahit",
        "-1", str(workdir / READS_1_FASTQ),
        "-2", str(workdir / READS_2_FASTQ),
        "-o", str(workdir / f"{MEGAHIT_OUT_PREFIX}{iter}-{subiter}"),
        "--k-list", "15,21,29,39,59,79,99,119,141",
    ]
    # fmt: on
    cmds.append(cmd)
    subiter += 1

    # More permissive assembly
    # fmt: off
    cmd = [
        "megahit",
        "-1", str(workdir / READS_1_FASTQ),
        "-2", str(workdir / READS_2_FASTQ),
        "-o", str(workdir / f"{MEGAHIT_OUT_PREFIX}{iter}-{subiter}"),
        "--min-count", "1",
        "--prune-level", "1",
        "--prune-depth", "1",
        "--min-contig-len", "100",
        "--k-list", "15,21,29,39,59,79,99,119,141",
    ]
    # fmt: on
    cmds.append(cmd)

    log_path = workdir / LOG_FILE
    for cmd in cmds:
        with open(log_path, "a") as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)


def _contig_ids_by_seed_ahocorasick(
    records: List[SeqRecord],
    seed_seqs: List[Seq],
) -> Dict[int, SeqOrientation]:
    """
    Find which contigs contain which seeds using Aho-Corasick algorithm.

    Args:
        records: List of SeqRecord objects (contigs)
        seed_seqs: List of seed sequences to search for

    Returns:
        Dict mapping contig index to orientation (FORWARD or REVERSE)

    Note:
        When a contig matches multiple seeds (or matches both forward and reverse
        complement), the orientation is determined by the leftmost match position
        in the contig.
    """

    if not seed_seqs or not records:
        return {}

    automaton = ahocorasick.Automaton()

    for seed in seed_seqs:
        seed_str = str(seed).upper()
        seed_rc = _reverse_complement(seed_str)

        automaton.add_word(seed_str, SeqOrientation.FORWARD)

        if seed_rc != seed_str:
            automaton.add_word(seed_rc, SeqOrientation.REVERSE)

    automaton.make_automaton()

    matches: Dict[int, SeqOrientation] = {}

    for contig_idx, record in enumerate(records):
        contig_seq = str(record.seq).upper()

        for _end_pos, orientation in automaton.iter(contig_seq):
            matches[contig_idx] = orientation
            break

    return matches


def _subset_contigs(
    workdir: PathLike,
    iter: int,
    seed_seqs: List[Seq],
    include_overlaps: bool = True,
    overlap_n0: int = 7,
    overlap_n1: int = 31,
) -> None:
    """Subset assembled contigs for the iteration to those containing seed sequences.

    For each workdir/megahit_out_iter<iter>-<subiter>/final.contigs.fa, creates subset
    workdir/megahit_out_iter<iter>-<subiter>/contigs_filtered.fasta containing just the
    contigs that contain the seed (or its reverse complement).

    Args:
        workdir: Working directory containing megahit output
        iter: Current iteration number
        seed_seqs: Seed sequences to search for
        include_overlaps: Whether to include contigs that do not contain a seed themselves,
            but are connected via overlapping sequences with a contig that does
        overlap_n0: Minimum overlap length for exact matches
        overlap_n1: Minimum overlap length when allowing 1 error
    """
    workdir = Path(workdir)
    subiter_dirs = [
        d
        for d in workdir.iterdir()
        if d.is_dir() and d.name.startswith(f"{MEGAHIT_OUT_PREFIX}{iter}-")
    ]
    # Find the assembly output of each megahit subiteration
    for subiter_dir in subiter_dirs:
        contigs_path = subiter_dir / MEGAHIT_FINAL_CONTIGS
        if not contigs_path.is_file():
            continue

        subset_path = subiter_dir / MEGAHIT_FILTERED_CONTIGS
        records: List[SeqRecord] = list(SeqIO.parse(contigs_path, "fasta"))

        filtered_records = []

        # Get the indices of all contigs that have seeds in them, along with their orientation
        # with respect to the seed.
        subsetted_ids_and_orientations: Dict[int, SeqOrientation] = (
            _contig_ids_by_seed_ahocorasick(records, seed_seqs)
        )

        if include_overlaps:
            seqs = [rec.seq for rec in records]
            subsetted_ids_and_orientations = get_overlapping_sequence_ids(
                seqs, subsetted_ids_and_orientations, overlap_n0, overlap_n1
            )

        for idx, orientation in subsetted_ids_and_orientations.items():
            record = records[idx]
            if record.seq is not None and orientation == SeqOrientation.REVERSE:
                # Contig is in the reverse direction with respect to the seed. We want to report
                # it as forward with respect to the seed, so take the reverse compliment
                record.seq = record.seq.reverse_complement()
            filtered_records.append(record)

        SeqIO.write(filtered_records, subset_path, "fasta")


def _choose_best_subiter(
    workdir: PathLike, iter: int, longest_thresh: int = 1, total_thresh: int = 30
) -> bool:
    """Choose a subiter's output to be current contigs for next iteration.

    A set of contigs is considered an improvement if either:
    - The longest contig grew by at least longest_thresh bases
    - The longest contig didn't get shorter and total contig length
      improved by at least total_thresh bases

    If an improvement is found, updates the current_contigs.fasta file with
    the best contigs and marks that subiter as chosen.

    Args:
        workdir: Working directory
        iter: Current iteration number
        longest_thresh: Minimum increase in longest contig length
        total_thresh: Minimum increase in total contig length

    Returns:
        Whether any subiter provided adequate improvement
    """
    workdir = Path(workdir)
    current_contigs = workdir / CURRENT_CONTIGS
    if not current_contigs.is_file():
        return False

    cur_stats = _fasta_longest_total(current_contigs)

    subiter_dirs = sorted(
        d
        for d in workdir.iterdir()
        if d.is_dir() and d.name.startswith(f"{MEGAHIT_OUT_PREFIX}{iter}-")
    )

    for subiter_dir in subiter_dirs:
        contigs_path = subiter_dir / MEGAHIT_FILTERED_CONTIGS
        if not contigs_path.is_file():
            continue

        subiter_stats = _fasta_longest_total(contigs_path)

        if subiter_stats.longest >= cur_stats.longest + longest_thresh or (
            subiter_stats.longest == cur_stats.longest
            and subiter_stats.total >= cur_stats.total + total_thresh
        ):
            # Found a good subiteration!
            # Copy rather than move to preserve original assembly output
            shutil.copy2(contigs_path, current_contigs)
            (subiter_dir / CHOSEN_SUBITER_FLAG).touch()
            return True

    return False


def _fasta_longest_total(path: PathLike) -> FastaStats:
    """Get longest and total length of sequences in a fasta file.

    Args:
        path: Path to fasta file

    Returns:
        Named tuple containing longest and total sequence lengths
    """
    path = Path(path)
    if not path.is_file():
        return FastaStats(0, 0)

    lengths = [len(rec.seq) for rec in SeqIO.parse(path, "fasta")]
    if lengths:
        return FastaStats(max(lengths), sum(lengths))
    return FastaStats(0, 0)


def _subset_split_files(
    s3_records: S3Files,
    ref_fasta_path: PathLike,
    read_subset_k: int,
    workdir: PathLike,
    use_batch: bool,
    *,
    batch_workdir: Optional[PathLike] = None,
    batch_queue: Optional[str] = None,
    tower_token: Optional[str] = None,
    num_parallel: Optional[int] = None,
    n_threads: int = 3,
    ordered: bool = True,
) -> None:
    """Determine whether to use batch or local Nucleaze to subset reads from split files.

    Args:
        s3_records: Processed S3 paths of read files
        ref_fasta_path: Path to reference fasta for filtering
        read_subset_k: Kmer size for matching
        workdir: Working directory for output
        use_batch: Whether to use batch mode
        batch_workdir: S3 path to batch work directory; Argument is ignored if use_batch is False
        batch_queue: Batch queue name; Argument is ignored if use_batch is False
        num_parallel: Number of parallel Nucleaze processes; Argument is ignored if use_batch is True
        n_threads: Threads per Nucleaze process; Argument is ignored if use_batch is True
        ordered: Force output order to match input read order; Argument is ignored if use_batch is True
    """
    if use_batch:
        _subset_split_files_batch(
            s3_records,
            ref_fasta_path,
            read_subset_k=read_subset_k,
            workdir=workdir,
            batch_workdir=batch_workdir,
            batch_queue=batch_queue,
            tower_token=tower_token,
        )
    else:
        _subset_split_files_local(
            s3_records,
            ref_fasta_path,
            read_subset_k=read_subset_k,
            workdir=workdir,
            num_parallel=num_parallel,
            n_threads=n_threads,
            ordered=ordered,
        )


def _subset_split_files_local(
    s3_records: S3Files,
    ref_fasta_path: PathLike,
    read_subset_k: int,
    workdir: PathLike,
    *,
    num_parallel: Optional[int] = None,
    n_threads: int = 3,
    ordered: bool = True,
) -> None:
    """Use parallel Nucleaze to subset reads from split files sharing kmers with ref. By
    default, the order of filtered reads will match the order of inputs, which makes
    the overall assembly algorithm closer to determinsitic.

    Args:
        s3_records: Processed S3 paths of read files
        ref_fasta_path: Path to reference fasta for filtering
        read_subset_k: Kmer size for matching
        workdir: Working directory for output
        num_parallel: Number of parallel Nucleaze processes
        n_threads: Threads per Nucleaze process
        ordered: Force output order to match input read order
    """
    ref_fasta_path = Path(ref_fasta_path)
    workdir = Path(workdir)
    if not ref_fasta_path.is_file():
        raise ValueError(f"Reference fasta does not exist: {ref_fasta_path}")
    if not workdir.is_dir():
        raise ValueError(f"Working directory does not exist: {workdir}")

    if num_parallel is None:
        num_parallel = max(1, cpu_count() // 4)

    nucleaze_bin = shutil.which("nucleaze")
    if nucleaze_bin is None:
        raise RuntimeError("nucleaze not found in PATH")

    # Build Nucleaze commands - equivalent to BBDuk with:
    #   rcomp=t -> --canonical (use canonical k-mer form)
    #   minkmerhits=1 -> --minhits 1 (default)
    #   mm=f -> (default, exact matching only)
    #   interleaved=t -> --interinput
    #   ordered=t -> --order
    # Note: nucleaze requires both --outm/--outm2 AND --outu/--outu2 when using paired input
    # We discard unmatched reads to /dev/null
    cmds = [
        f"{s3_stream_cmd(rec.s3_path)} | "
        f"zstdcat - | "
        f"{nucleaze_bin} --in - "
        f"--outm {workdir / rec.filename}_1.fastq --outm2 {workdir / rec.filename}_2.fastq "
        f"--outu /dev/null --outu2 /dev/null "
        f"--ref {ref_fasta_path} --k {read_subset_k} "
        f"--canonical --minhits 1 --interinput "
        f"{'--order ' if ordered else ''}"
        f"--threads {n_threads}"
        for rec in s3_records
    ]

    def run_cmd(cmd: str) -> subprocess.CompletedProcess:
        return subprocess.run(cmd, shell=True, check=True, capture_output=True)

    with ThreadPoolExecutor(max_workers=num_parallel) as executor:
        futures = [executor.submit(run_cmd, cmd) for cmd in cmds]
        for future in as_completed(futures):
            future.result()  # Raises exception if command failed

    # Concatenate per-split hits
    for read_num in (1, 2):
        output_path = workdir / f"reads_{read_num}.fastq"
        split_files = [workdir / f"{rec.filename}_{read_num}.fastq" for rec in s3_records]
        concat_and_tag_fastq(split_files, output_path)
        for split_file in split_files:
            split_file.unlink()


def _create_nextflow_config(
    workdir: Path,
    batch_workdir: PathLike,
    s3_records_file: Path,
    ref_fasta_path: Path,
    read_subset_k: int,
    batch_queue: str,
    tower_token: str,
) -> Path:
    """Create Nextflow configuration file for a specific run.

    Args:
        workdir: Working directory for this run. Overwrites <workdir>/run_config.nf
        batch_workdir: S3 path for Nextflow work directory
        s3_records_file: Path to file containing S3 paths
        ref_fasta_path: Path to reference fasta
        read_subset_k: K-mer size for filtering
        batch_queue: AWS Batch queue name
        tower_token: Seqera Tower token

    Returns:
        Path to created configuration file
    """
    nextflow_dir = Path(__file__).resolve().parent.parent / "nextflow"

    config_content = textwrap.dedent(
        f"""
        // Dynamic configuration for this run
        params {{
            base_dir = "{batch_workdir}"
            s3_files = "{s3_records_file}"
            ref_fasta_path = "{ref_fasta_path}"
            kmer = {read_subset_k}
        }}
        tower.accessToken = "{tower_token}"
        
        // Process configuration
        process.queue = '{batch_queue}'
        
        // Include all static configuration files from the OA repo
        includeConfig "{nextflow_dir}/static_configs.config"
        """
    ).strip()

    config_path = workdir / "run_config.nf"
    with open(config_path, "w") as f:
        f.write(config_content)

    return config_path


def _subset_split_files_batch(
    s3_records: S3Files,
    ref_fasta_path: PathLike,
    read_subset_k: int,
    workdir: PathLike,
    batch_workdir: PathLike,
    batch_queue: str,
    tower_token: str,
) -> None:
    """Use Nextflow (AWS Batch) to subset reads from split files sharing kmers with ref. By
    default, the order of filtered reads will match the order of inputs, which makes
    the overall assembly algorithm closer to deterministic.

    Args:
        s3_records: Processed S3 paths of read files
        ref_fasta_path: Path to reference fasta for filtering
        read_subset_k: Kmer size for matching
        workdir: Working directory for output
        batch_workdir: S3 path to batch work directory
        batch_queue: Batch queue name
        tower_token: Seqera Tower access token
    """
    ref_fasta_path = Path(ref_fasta_path)
    workdir = Path(workdir)
    if not ref_fasta_path.is_file():
        raise ValueError(f"Reference fasta does not exist: {ref_fasta_path}")
    if not workdir.is_dir():
        raise ValueError(f"Working directory does not exist: {workdir}")

    # Write s3_records to a file as filename tab s3_path
    s3_records_file = workdir / "s3_records.txt"
    with open(s3_records_file, "w") as f:
        for rec in s3_records:
            f.write(f"{rec.filename}\t{rec.s3_path}\n")

    # Create dynamic configuration
    dynamic_config = _create_nextflow_config(
        workdir=workdir,
        batch_workdir=batch_workdir,
        s3_records_file=s3_records_file,
        ref_fasta_path=ref_fasta_path,
        read_subset_k=read_subset_k,
        batch_queue=batch_queue,
        tower_token=tower_token,
    )

    # Get paths to Nextflow components
    nextflow_dir = Path(__file__).resolve().parent.parent / "nextflow"
    nextflow_main = nextflow_dir / "main.nf"

    # Run Nextflow with single dynamic config (which includes static configs)
    profile = os.environ.get(NF_PROFILE_ENV_VAR, "standard")
    nextflow_cmd = [
        "nextflow",
        "run",
        str(nextflow_main),
        "-c",
        str(dynamic_config),
        "-profile",
        profile,
    ]

    result = subprocess.run(nextflow_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Nextflow stdout: {result.stdout}")
        print(f"Nextflow stderr: {result.stderr}")
        raise ValueError(f"Nextflow run failed with exit code {result.returncode}")
    else:
        print("Nextflow run succeeded")

    # Copy reads from s3 back to workdir; catch errors
    forward_read = subprocess.run(
        [
            "aws",
            "s3",
            "cp",
            f"{batch_workdir}/output/results/reads_1.fastq",
            workdir / "reads_1.fastq",
        ]
    )
    reverse_read = subprocess.run(
        [
            "aws",
            "s3",
            "cp",
            f"{batch_workdir}/output/results/reads_2.fastq",
            workdir / "reads_2.fastq",
        ]
    )

    if forward_read.returncode != 0 or reverse_read.returncode != 0:
        raise ValueError("Failed to copy reads from s3")
    else:
        print("Reads copied successfully")


def _copy_iteration_reads(workdir: PathLike, iter_num: int) -> None:
    """Copy read files from an iteration into debug directory structure.

    Creates workdir/reads/iter_<N> and copies both standard and frequency-filtered
    read files if they exist.

    Args:
        workdir: Working directory containing read files
        iter_num: Iteration number for directory naming
    """
    workdir = Path(workdir)
    reads_dir = workdir / "reads"
    reads_dir.mkdir(exist_ok=True)
    iter_dir = reads_dir / f"iter_{iter_num}"
    iter_dir.mkdir()

    read_files = [
        READS_1_FASTQ,
        READS_2_FASTQ,
        READS_UNTRIMMED_1_FASTQ,
        READS_UNTRIMMED_2_FASTQ,
        READS_FILTERED_1_FASTQ,
        READS_FILTERED_2_FASTQ,
    ]
    for file in read_files:
        if (workdir / file).exists():
            shutil.copy2(workdir / file, iter_dir / file)


def _frequency_filter_reads(workdir: Path, high_freq_kmers_path: Path, k: int) -> None:
    """Create filtered reads for this iteration's assemblies"""
    nucleaze_bin = shutil.which("nucleaze")
    if nucleaze_bin is None:
        raise RuntimeError("nucleaze not found in PATH")

    with open(workdir / LOG_FILE, "a") as log:
        # Use Nucleaze to filter out reads containing high-frequency kmers
        # --outu outputs reads that do NOT match the reference (unmatched reads)
        # Note: nucleaze requires --outm/--outm2 when using --in2, so we discard matched reads to /dev/null
        subprocess.run(
            [
                nucleaze_bin,
                "--in",
                str(workdir / READS_1_FASTQ),
                "--in2",
                str(workdir / READS_2_FASTQ),
                "--outu",
                str(workdir / READS_FILTERED_1_FASTQ),
                "--outu2",
                str(workdir / READS_FILTERED_2_FASTQ),
                "--outm",
                "/dev/null",
                "--outm2",
                "/dev/null",
                "--ref",
                str(high_freq_kmers_path),
                "--k",
                str(k),
                "--canonical",
                "--minhits",
                "1",
                "--order",
                "--threads",
                "4",
            ],
            stdout=log,
            stderr=log,
            check=True,
        )


def _extract_fasta_kmers(fasta_path: Path, out_path: Path, k: int) -> None:
    """Extract all k-mers from sequences in a fasta file."""
    kmers = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        kmers.extend(seq[i : i + k] for i in range(len(seq) - k + 1))
    with open(out_path, "w") as f:
        for i, kmer in enumerate(kmers):
            f.write(f">kmer_{i + 1}\n{kmer}\n")


def _adapter_sanitize_kmer_queries(
    workdir: Path, contigs_path: Path, adapters_path: Path, k: int
) -> Path:
    """Remove adapter kmers from set of kmers we'll search for in all reads.
    Return the path to the sanitized query kmers."""
    # We don't want to remove entire contigs that contain an adapter kmer, so first we
    # generate all the kmers from the contigs, then filter that set of kmers to not
    # contain adapter kmers.
    candidate_kmers_path = workdir / CANDIDATE_KMERS
    final_query_kmers_path = workdir / QUERY_KMERS

    _extract_fasta_kmers(contigs_path, candidate_kmers_path, k)

    # Check if adapters file has sequences long enough to produce k-mers
    # If not, skip filtering and just use all candidate kmers
    has_valid_adapter_kmers = False
    for record in SeqIO.parse(adapters_path, "fasta"):
        if len(record.seq) >= k:
            has_valid_adapter_kmers = True
            break

    if not has_valid_adapter_kmers:
        # No adapter k-mers possible, just copy input to output
        shutil.copy(candidate_kmers_path, final_query_kmers_path)
        candidate_kmers_path.unlink()
        return final_query_kmers_path

    # Use Nucleaze to remove adapter kmers from candidate kmers
    # --outu outputs unmatched sequences (those without adapter kmers)
    nucleaze_bin = shutil.which("nucleaze")
    if nucleaze_bin is None:
        raise RuntimeError("nucleaze not found in PATH")

    cmd_parts = [
        nucleaze_bin,
        "--in",
        str(candidate_kmers_path),
        "--outu",
        str(final_query_kmers_path),  # we want the unmatched candidate kmers
        "--ref",
        str(adapters_path),
        "--k",
        str(k),
        "--canonical",
        "--minhits",
        "1",
        "--order",
        "--threads",
        "4",
    ]
    with open(workdir / LOG_FILE, "a") as log:
        subprocess.run(cmd_parts, shell=False, check=True, stdout=log, stderr=log)
    candidate_kmers_path.unlink()
    return final_query_kmers_path


def _prepare_query_seqs(
    workdir: Path,
    warm_start_path: Optional[PathLike],
    adapters_path: Optional[PathLike],
    read_subset_k: int,
) -> Path:
    """Prepare query sequences used for read filtering:
        * If provided, add warm start contigs to current contigs, otherwise just use the
            current contigs. Typically only done in first iteration when current contigs
            is just the seed.
        * If adapters are provided, sanitize query sequences so they do not contain
        adapter kmers.

    Return the path to the query sequences for read filtering (which may just be the
    unmodified current contigs)."""
    if warm_start_path is not None:
        # Combine seed (current_contigs) with warm start sequences
        warm_start_combined = workdir / WARM_START
        with open(workdir / CURRENT_CONTIGS, "r") as f:
            seed_content = f.read()
        with open(warm_start_path, "r") as f:
            warm_content = f.read()
        with open(warm_start_combined, "w") as outfile:
            outfile.write(seed_content.rstrip("\n"))
            outfile.write("\n")  # Exactly one newline between files
            outfile.write(warm_content)
        query_seqs_path = warm_start_combined
    else:
        query_seqs_path = workdir / CURRENT_CONTIGS
    if adapters_path is not None:
        kmer_ref_path = _adapter_sanitize_kmer_queries(
            workdir, query_seqs_path, Path(adapters_path), read_subset_k
        )
    else:
        kmer_ref_path = query_seqs_path
    return kmer_ref_path


def _adapter_trim_iter_reads(workdir: Path, adapters_path: Path):
    """Remove adapters from reads we'll assemble this iteration. At present this does
    adapter trimming as well as a gentler version of the cleaning typical of fastp in
    the MGS pipeline, ideally so outward assembly can safely run on raw data."""
    # subsetted reads are in the regular reads_{1/2}.fastq; rename so we can clean
    (workdir / READS_1_FASTQ).rename(workdir / READS_UNTRIMMED_1_FASTQ)
    (workdir / READS_2_FASTQ).rename(workdir / READS_UNTRIMMED_2_FASTQ)
    # fmt: off
    with open(workdir / LOG_FILE, "a") as log:
        subprocess.run(
            [
                "fastp",
                "--in1", str(workdir / READS_UNTRIMMED_1_FASTQ),
                "--in2", str(workdir / READS_UNTRIMMED_2_FASTQ),
                "--out1", str(workdir / READS_1_FASTQ),
                "--out2", str(workdir / READS_2_FASTQ),
                "--adapter_fasta", str(adapters_path),
                "--cut_front",
                "--cut_tail",
                "--trim_poly_x",
                "--cut_mean_quality", "10",
                "--average_qual", "10",
                "--qualified_quality_phred", "10",
                "--low_complexity_filter",
                "--dont_eval_duplication",
                "--thread", "4",
            ],
            stdout=log,
            stderr=log,
            check=True,
        )
    # fmt: on


def _record_inputs(workdir: Path, seed_path: Path, s3_paths: List[str]) -> None:
    """Record inputs in a working directory.
    To-do: save other arguments to outward_assembly?"""
    shutil.copy2(seed_path, workdir / "original_seed.fasta")
    with open(workdir / "input_s3_paths.txt", "w+") as f:
        f.write("\n".join(s3_paths))
