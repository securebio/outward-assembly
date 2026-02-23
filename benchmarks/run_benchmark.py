#!/usr/bin/env python3
"""Benchmark script to compare pipeline performance between BBDuk and Nucleaze.

This script runs the outward assembly pipeline locally and measures:
- Execution time
- Output correctness (contig sequences)

Usage:
    # From repository root, with oa-tools environment activated:
    mamba activate oa-tools
    uv run python benchmarks/run_benchmark.py

Results are saved to benchmarks/results/<branch_name>_<timestamp>/
"""

import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from outward_assembly.io_helpers import process_s3_paths
from outward_assembly.pipeline import outward_assembly


def get_git_branch() -> str:
    """Get current git branch name."""
    result = subprocess.run(
        ["git", "branch", "--show-current"],
        capture_output=True,
        text=True,
        cwd=Path(__file__).parent.parent,
    )
    return result.stdout.strip() or "unknown"


def hash_file(path: Path) -> str:
    """Get MD5 hash of file contents."""
    if not path.exists():
        return "FILE_NOT_FOUND"
    return hashlib.md5(path.read_bytes()).hexdigest()


def run_pipeline_benchmark(
    output_dir: Path,
    s3_paths: list[str],
    seed_path: Path,
    high_freq_kmers_path: Path | None = None,
    adapters_path: Path | None = None,
    warm_start_path: Path | None = None,
    max_iters: int = 2,
) -> dict:
    """Run the pipeline and collect timing/output metrics."""
    work_dir = output_dir / "workdir"
    work_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / "final_contigs.fasta"

    kwargs = {
        "s3_paths": s3_paths,
        "seed_path": seed_path,
        "output_path": output_path,
        "max_iters": max_iters,
        "work_dir_parent": work_dir,
        "use_batch": False,  # Local profile
        "cleanup": False,  # Keep intermediate files for inspection
    }

    if high_freq_kmers_path:
        kwargs["high_freq_kmers_path"] = high_freq_kmers_path
        kwargs["freq_filter_k"] = 23

    if adapters_path:
        kwargs["adapters_path"] = adapters_path

    if warm_start_path:
        kwargs["warm_start_path"] = warm_start_path

    # Run pipeline with timing
    start_time = time.perf_counter()
    try:
        outward_assembly(**kwargs)
        success = True
        error = None
    except Exception as e:
        success = False
        error = str(e)
    end_time = time.perf_counter()

    elapsed = end_time - start_time

    # Collect results
    results = {
        "elapsed_seconds": elapsed,
        "success": success,
        "error": error,
        "output_exists": output_path.exists(),
        "output_hash": hash_file(output_path),
        "output_size_bytes": output_path.stat().st_size if output_path.exists() else 0,
    }

    # Read output contigs if they exist
    if output_path.exists():
        results["output_content"] = output_path.read_text()

    return results


def main():
    branch = get_git_branch()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Setup paths - use /tmp for shorter paths (avoids xargs length limits on main)
    repo_root = Path(__file__).parent.parent
    data_dir = repo_root / "tests/data/pipeline-end-to-end"
    results_dir = Path(f"/tmp/oa-bench/{branch}_{timestamp}")
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Outward Assembly Benchmark ===")
    print(f"Branch: {branch}")
    print(f"Results: {results_dir}")
    print()

    # Test data - using the pipeline end-to-end test data
    s3_paths = [
        "s3://nao-testing/outward-assembly-test-data/pipeline-end-to-end/simulated-abcbd-reads.fastq.zst"
    ]

    # Run benchmark scenarios
    # Note: adapters_path omitted to skip fastp dependency (not relevant to BBDuk vs Nucleaze comparison)
    scenarios = {
        "basic": {
            "seed_path": data_dir / "seed_seq.fasta",
            "high_freq_kmers_path": data_dir / "fake_hfk.fasta",
            "adapters_path": None,  # Skip adapter trimming (requires fastp)
            "warm_start_path": None,  # Skip warm start (requires adapters)
            "max_iters": 2,
        },
        "multi_seed": {
            "seed_path": data_dir / "multi_seed_seq.fasta",
            "high_freq_kmers_path": data_dir / "fake_hfk.fasta",
            "adapters_path": None,  # Skip adapter trimming
            "warm_start_path": None,  # Skip warm start
            "max_iters": 2,
        },
    }

    all_results = {
        "branch": branch,
        "timestamp": timestamp,
        "scenarios": {},
    }

    for scenario_name, scenario_kwargs in scenarios.items():
        print(f"Running scenario: {scenario_name}...")
        scenario_dir = results_dir / scenario_name
        scenario_dir.mkdir(parents=True, exist_ok=True)

        results = run_pipeline_benchmark(
            output_dir=scenario_dir,
            s3_paths=s3_paths,
            **scenario_kwargs,
        )

        all_results["scenarios"][scenario_name] = results

        status = "SUCCESS" if results["success"] else "FAILED"
        print(f"  {status} in {results['elapsed_seconds']:.2f}s")
        if results["error"]:
            print(f"  Error: {results['error']}")
        print()

    # Save results summary
    summary_path = results_dir / "summary.json"
    with open(summary_path, "w") as f:
        # Don't include full output content in summary
        summary = {
            "branch": all_results["branch"],
            "timestamp": all_results["timestamp"],
            "scenarios": {
                name: {k: v for k, v in data.items() if k != "output_content"}
                for name, data in all_results["scenarios"].items()
            },
        }
        json.dump(summary, f, indent=2)

    print(f"Results saved to: {summary_path}")

    # Print summary table
    print("\n=== Summary ===")
    print(f"{'Scenario':<15} {'Status':<10} {'Time (s)':<12} {'Output Hash':<34}")
    print("-" * 75)
    for name, data in all_results["scenarios"].items():
        status = "OK" if data["success"] else "FAIL"
        time_str = f"{data['elapsed_seconds']:.2f}"
        hash_str = data["output_hash"][:32] if data["output_hash"] else "N/A"
        print(f"{name:<15} {status:<10} {time_str:<12} {hash_str:<34}")


if __name__ == "__main__":
    main()
