"""
Verify that kmer filtering with nucleaze produces the same results as the original bbduk implementation.
"""

import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest

from outward_assembly.io_helpers import _count_lines


def _hash_file(path: Path) -> str:
    """Get MD5 hash of file contents."""
    return hashlib.md5(path.read_bytes()).hexdigest()


def _deinterleave_fastq(interleaved: Path, out_1: Path, out_2: Path) -> None:
    """Deinterleave a paired-end FASTQ file into two separate files."""
    with open(interleaved) as f, open(out_1, "w") as f1, open(out_2, "w") as f2:
        while True:
            lines = [f.readline() for _ in range(4)]
            if not lines[0]:
                break
            f1.writelines(lines)
            lines = [f.readline() for _ in range(4)]
            f2.writelines(lines)


def _run_nucleaze(
    reads_1: Path,
    reads_2: Path,
    ref_path: Path,
    out_1: Path,
    out_2: Path,
    k: int = 23,
) -> None:
    """Run nucleaze k-mer filtering."""
    nucleaze_bin = shutil.which("nucleaze")
    if nucleaze_bin is None:
        pytest.skip("nucleaze not found in PATH")

    cmd = [
        nucleaze_bin,
        "--in", str(reads_1),
        "--in2", str(reads_2),
        "--outu", str(out_1),
        "--outu2", str(out_2),
        "--outm", "/dev/null",
        "--outm2", "/dev/null",
        "--ref", str(ref_path),
        "--k", str(k),
        "--canonical",
        "--minhits", "1",
        "--order",
        "--threads", "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"nucleaze failed: {result.stderr}")


def _run_bbduk(
    reads_1: Path,
    reads_2: Path,
    ref_path: Path,
    out_1: Path,
    out_2: Path,
    k: int = 23,
) -> None:
    """Run BBDuk k-mer filtering with equivalent parameters to nucleaze."""
    bbduk_bin = shutil.which("bbduk.sh")
    if bbduk_bin is None:
        pytest.skip("bbduk.sh not found in PATH")

    cmd = [
        bbduk_bin,
        f"in1={reads_1}",
        f"in2={reads_2}",
        f"outu1={out_1}",
        f"outu2={out_2}",
        f"ref={ref_path}",
        f"k={k}",
        "mm=f",
        "hdist=0",
        "rcomp=t",
        "ordered=t",
    ]
    # BBDuk writes all output (including normal logging) to stderr
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Look for actual error indicators in stderr
        stderr = result.stderr
        if "Exception" in stderr or "Error" in stderr:
            raise RuntimeError(f"bbduk failed (exit {result.returncode}): {stderr}")
        raise RuntimeError(f"bbduk failed with exit code {result.returncode}")


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.requires_tools
@pytest.mark.requires_aws
def test_nucleaze_bbduk_output_equivalence(temp_workdir):
    """Test that nucleaze produces identical output to BBDuk."""
    # Download test data
    s3_path = "s3://nao-testing/outward-assembly-test-data/pipeline-end-to-end/simulated-abcbd-reads.fastq.zst"
    zst_file = temp_workdir / "reads.fastq.zst"
    subprocess.run(
        ["aws", "s3", "cp", s3_path, str(zst_file)],
        check=True,
        capture_output=True,
    )

    # Decompress
    fastq_file = temp_workdir / "reads.fastq"
    subprocess.run(
        ["zstd", "-d", str(zst_file), "-o", str(fastq_file)],
        check=True,
        capture_output=True,
    )

    # Deinterleave into paired files
    reads_1 = temp_workdir / "reads_1.fastq"
    reads_2 = temp_workdir / "reads_2.fastq"
    _deinterleave_fastq(fastq_file, reads_1, reads_2)

    # Reference k-mers file
    ref_path = Path(__file__).parent.parent / "data/pipeline-end-to-end/fake_hfk.fasta"

    # nucleaze
    nucleaze_out_1 = temp_workdir / "nucleaze_out_1.fastq"
    nucleaze_out_2 = temp_workdir / "nucleaze_out_2.fastq"
    _run_nucleaze(reads_1, reads_2, ref_path, nucleaze_out_1, nucleaze_out_2)

    # BBDuk
    bbduk_out_1 = temp_workdir / "bbduk_out_1.fastq"
    bbduk_out_2 = temp_workdir / "bbduk_out_2.fastq"
    _run_bbduk(reads_1, reads_2, ref_path, bbduk_out_1, bbduk_out_2)

    # Verify read counts match
    nucleaze_reads = _count_lines(nucleaze_out_1) // 4
    bbduk_reads = _count_lines(bbduk_out_1) // 4
    assert nucleaze_reads == bbduk_reads, (
        f"Read count mismatch: nucleaze={nucleaze_reads}, bbduk={bbduk_reads}"
    )

    # Verify output hashes match
    assert _hash_file(nucleaze_out_1) == _hash_file(bbduk_out_1), "Output mismatch for read 1"
    assert _hash_file(nucleaze_out_2) == _hash_file(bbduk_out_2), "Output mismatch for read 2"
