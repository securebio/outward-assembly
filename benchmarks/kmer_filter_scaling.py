#!/usr/bin/env python3
"""Scaling benchmark comparing BBDuk vs Nucleaze k-mer filtering at different input sizes.

Duplicates the test reads to create larger inputs and measures how performance scales.

Usage:
    source /Users/lee/miniforge3/etc/profile.d/conda.sh && conda activate oa-tools
    export PATH="$HOME/.cargo/bin:$PATH"
    uv run python benchmarks/kmer_filter_scaling.py
"""

import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def count_fastq_reads(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path) as f:
        return sum(1 for _ in f) // 4


def download_and_prepare(tmpdir: Path) -> tuple[Path, Path]:
    """Download and deinterleave test data."""
    s3_path = "s3://nao-testing/outward-assembly-test-data/pipeline-end-to-end/simulated-abcbd-reads.fastq.zst"
    zst_file = tmpdir / "reads.fastq.zst"

    print("  Downloading test data...")
    subprocess.run(["aws", "s3", "cp", s3_path, str(zst_file)], check=True, capture_output=True)

    fastq_file = tmpdir / "reads.fastq"
    subprocess.run(["zstd", "-d", str(zst_file), "-o", str(fastq_file)], check=True, capture_output=True)

    reads_1 = tmpdir / "base_reads_1.fastq"
    reads_2 = tmpdir / "base_reads_2.fastq"

    with open(fastq_file) as f, open(reads_1, "w") as f1, open(reads_2, "w") as f2:
        while True:
            lines = [f.readline() for _ in range(4)]
            if not lines[0]:
                break
            f1.writelines(lines)
            lines = [f.readline() for _ in range(4)]
            f2.writelines(lines)

    return reads_1, reads_2


def duplicate_reads(base_1: Path, base_2: Path, multiplier: int, tmpdir: Path) -> tuple[Path, Path]:
    """Create larger input by duplicating reads with unique names."""
    out_1 = tmpdir / f"reads_{multiplier}x_1.fastq"
    out_2 = tmpdir / f"reads_{multiplier}x_2.fastq"

    base_lines_1 = base_1.read_text().splitlines()
    base_lines_2 = base_2.read_text().splitlines()

    with open(out_1, "w") as f1, open(out_2, "w") as f2:
        for rep in range(multiplier):
            for i in range(0, len(base_lines_1), 4):
                # Rename header to avoid duplicate read names
                header_1 = base_lines_1[i].split()[0] + f"_rep{rep}"
                f1.write(header_1 + "\n")
                f1.write(base_lines_1[i + 1] + "\n")
                f1.write("+\n")
                f1.write(base_lines_1[i + 3] + "\n")

                header_2 = base_lines_2[i].split()[0] + f"_rep{rep}"
                f2.write(header_2 + "\n")
                f2.write(base_lines_2[i + 1] + "\n")
                f2.write("+\n")
                f2.write(base_lines_2[i + 3] + "\n")

    return out_1, out_2


def get_peak_memory_mb(cmd: list[str]) -> tuple[subprocess.CompletedProcess, float]:
    """Run command via /usr/bin/time to capture peak RSS in MB."""
    # Use /usr/bin/time -l on macOS to get peak RSS
    time_cmd = ["/usr/bin/time", "-l"] + cmd
    result = subprocess.run(time_cmd, capture_output=True, text=True)
    peak_mb = 0.0
    for line in result.stderr.splitlines():
        if "maximum resident set size" in line.lower():
            # macOS: bytes, Linux: KB
            val = int(line.strip().split()[0])
            peak_mb = val / (1024 * 1024)  # bytes to MB on macOS
            break
    return result, peak_mb


def run_nucleaze(reads_1, reads_2, ref_path, out_1, out_2, k=23):
    nucleaze_bin = shutil.which("nucleaze")
    if not nucleaze_bin:
        return None
    cmd = [
        nucleaze_bin,
        "--in", str(reads_1), "--in2", str(reads_2),
        "--outu", str(out_1), "--outu2", str(out_2),
        "--outm", "/dev/null", "--outm2", "/dev/null",
        "--ref", str(ref_path), "--k", str(k),
        "--canonical", "--minhits", "1", "--order", "--threads", "4",
    ]
    start = time.perf_counter()
    result, peak_mb = get_peak_memory_mb(cmd)
    elapsed = time.perf_counter() - start
    return {"elapsed_seconds": elapsed, "success": result.returncode == 0, "peak_memory_mb": peak_mb}


def run_bbduk(reads_1, reads_2, ref_path, out_1, out_2, k=23):
    bbduk_bin = shutil.which("bbduk.sh")
    if not bbduk_bin:
        return None
    # BBDuk is a shell script that launches Java, so /usr/bin/time captures the shell.
    # Instead, measure via the java process directly by wrapping.
    cmd = [
        bbduk_bin,
        f"in={reads_1}", f"in2={reads_2}",
        f"out={out_1}", f"out2={out_2}",
        f"ref={ref_path}", f"k={k}",
        "mm=f", "hdist=0", "rcomp=t", "ordered=t",
        "-Xmx2g",
    ]
    start = time.perf_counter()
    result, peak_mb = get_peak_memory_mb(cmd)
    elapsed = time.perf_counter() - start
    return {"elapsed_seconds": elapsed, "success": result.returncode == 0, "peak_memory_mb": peak_mb}


def main():
    repo_root = Path(__file__).parent.parent
    ref_path = repo_root / "tests/data/pipeline-end-to-end/fake_hfk.fasta"
    multipliers = [1, 10, 50, 100, 500, 1000]

    print("=== K-mer Filter Scaling Benchmark ===")
    print(f"Multipliers: {multipliers}")
    print()

    results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        print("Preparing base test data...")
        base_1, base_2 = download_and_prepare(tmpdir)
        base_reads = count_fastq_reads(base_1)
        print(f"  Base: {base_reads} read pairs")
        print()

        for mult in multipliers:
            num_reads = base_reads * mult
            print(f"--- {mult}x = {num_reads:,} read pairs ---")

            if mult == 1:
                r1, r2 = base_1, base_2
            else:
                print(f"  Generating {mult}x input...")
                r1, r2 = duplicate_reads(base_1, base_2, mult, tmpdir)

            actual_reads = count_fastq_reads(r1)
            row = {"multiplier": mult, "input_reads": actual_reads}

            # Nucleaze
            nuc_out_1 = tmpdir / f"nuc_{mult}_1.fq"
            nuc_out_2 = tmpdir / f"nuc_{mult}_2.fq"
            nuc = run_nucleaze(r1, r2, ref_path, nuc_out_1, nuc_out_2)
            if nuc and nuc["success"]:
                row["nucleaze_seconds"] = nuc["elapsed_seconds"]
                row["nucleaze_peak_mb"] = nuc["peak_memory_mb"]
                row["nucleaze_output_reads"] = count_fastq_reads(nuc_out_1)
                print(f"  Nucleaze: {nuc['elapsed_seconds']:.3f}s, {nuc['peak_memory_mb']:.0f} MB peak")
            else:
                print(f"  Nucleaze: FAILED or not found")

            # Clean up nucleaze output to save disk
            for f in [nuc_out_1, nuc_out_2]:
                if f.exists():
                    f.unlink()

            # BBDuk
            bb_out_1 = tmpdir / f"bb_{mult}_1.fq"
            bb_out_2 = tmpdir / f"bb_{mult}_2.fq"
            bb = run_bbduk(r1, r2, ref_path, bb_out_1, bb_out_2)
            if bb and bb["success"]:
                row["bbduk_seconds"] = bb["elapsed_seconds"]
                row["bbduk_peak_mb"] = bb["peak_memory_mb"]
                row["bbduk_output_reads"] = count_fastq_reads(bb_out_1)
                print(f"  BBDuk:    {bb['elapsed_seconds']:.3f}s, {bb['peak_memory_mb']:.0f} MB peak")
            else:
                print(f"  BBDuk:    FAILED or not found")

            # Clean up bbduk output
            for f in [bb_out_1, bb_out_2]:
                if f.exists():
                    f.unlink()

            if "nucleaze_seconds" in row and "bbduk_seconds" in row:
                speedup = row["bbduk_seconds"] / row["nucleaze_seconds"]
                mem_ratio = row["bbduk_peak_mb"] / max(row["nucleaze_peak_mb"], 1)
                print(f"  Speedup:  {speedup:.1f}x faster, {mem_ratio:.0f}x less memory")

            # Clean up duplicated input
            if mult > 1:
                r1.unlink()
                r2.unlink()

            results.append(row)
            print()

    # Print summary table
    print("=== Summary ===")
    print(f"{'Reads':>10}  {'Nucleaze (s)':>12}  {'BBDuk (s)':>12}  {'Speedup':>8}  {'Nuc MB':>8}  {'BBDuk MB':>8}  {'Mem ratio':>10}")
    print("-" * 85)
    for row in results:
        reads = f"{row['input_reads']:,}"
        nuc_t = f"{row.get('nucleaze_seconds', 0):.3f}" if "nucleaze_seconds" in row else "N/A"
        bb_t = f"{row.get('bbduk_seconds', 0):.3f}" if "bbduk_seconds" in row else "N/A"
        nuc_m = f"{row.get('nucleaze_peak_mb', 0):.0f}" if "nucleaze_peak_mb" in row else "N/A"
        bb_m = f"{row.get('bbduk_peak_mb', 0):.0f}" if "bbduk_peak_mb" in row else "N/A"
        if "nucleaze_seconds" in row and "bbduk_seconds" in row:
            speedup = f"{row['bbduk_seconds'] / row['nucleaze_seconds']:.1f}x"
        else:
            speedup = "N/A"
        if "nucleaze_peak_mb" in row and "bbduk_peak_mb" in row and row["nucleaze_peak_mb"] > 0:
            mem_r = f"{row['bbduk_peak_mb'] / row['nucleaze_peak_mb']:.0f}x"
        else:
            mem_r = "N/A"
        print(f"{reads:>10}  {nuc_t:>12}  {bb_t:>12}  {speedup:>8}  {nuc_m:>8}  {bb_m:>8}  {mem_r:>10}")

    # Save JSON
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = repo_root / f"benchmarks/results/scaling_{timestamp}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
