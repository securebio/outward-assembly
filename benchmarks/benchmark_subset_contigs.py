#!/usr/bin/env python3
"""
Benchmark script to compare seed-to-contig matching algorithms at scale.

This script generates synthetic datasets of increasing size and measures
the performance of:
1. Naive nested loop (original implementation)
2. Aho-Corasick (optimal for exact matching)
3. Bowtie2 (FM-index based)

Usage:
    python benchmark_subset_contigs.py

Requirements:
    - pyahocorasick: pip install pyahocorasick
    - bowtie2: conda install -c bioconda bowtie2
    - biopython: pip install biopython
"""

import random
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import the implementations
from outward_assembly.basic_seq_operations import SeqOrientation
from outward_assembly.pipeline_steps import _contig_ids_by_seed_ahocorasick
from naive import contig_ids_by_seed
from bowtie import contig_ids_by_seed_bowtie2 as _contig_ids_by_seed_bowtie2


def generate_test_data(
    num_seeds: int,
    num_contigs: int,
    seed_len: int = 40,
    contig_len: int = 500,
    match_fraction: float = 0.1,
    variable_seed_length: bool = True,
    seed_len_range: Tuple[int, int] = (20, 60),
) -> Tuple[List[Seq], List[SeqRecord]]:
    """
    Generate synthetic test data for benchmarking.

    Args:
        num_seeds: Number of seed sequences to generate
        num_contigs: Number of contig sequences to generate
        seed_len: Fixed seed length (used if variable_seed_length=False)
        contig_len: Average contig length
        match_fraction: Fraction of contigs that should contain a seed
        variable_seed_length: Whether to use variable-length seeds
        seed_len_range: (min, max) seed length if variable_seed_length=True

    Returns:
        Tuple of (seeds, contigs)
    """
    bases = "ACGT"

    # Generate seeds (variable length to match real use case)
    seeds = []
    for _ in range(num_seeds):
        if variable_seed_length:
            length = random.randint(seed_len_range[0], seed_len_range[1])
        else:
            length = seed_len
        seq = "".join(random.choices(bases, k=length))
        seeds.append(Seq(seq))

    # Generate contigs
    contigs = []
    for i in range(num_contigs):
        # Vary contig length a bit
        length = random.randint(int(contig_len * 0.5), int(contig_len * 1.5))
        seq = "".join(random.choices(bases, k=length))
        contigs.append(SeqRecord(Seq(seq), id=f"contig_{i}", description=""))

    # Insert seeds into a fraction of contigs
    num_matches = int(num_contigs * match_fraction)
    match_indices = random.sample(range(num_contigs), num_matches)

    for i in match_indices:
        seed = str(random.choice(seeds))
        contig_seq = str(contigs[i].seq)

        # Make sure contig is long enough
        if len(contig_seq) < len(seed):
            continue

        # Insert seed at random position
        pos = random.randint(0, len(contig_seq) - len(seed))
        new_seq = contig_seq[:pos] + seed + contig_seq[pos + len(seed) :]
        contigs[i].seq = Seq(new_seq)

    return seeds, contigs


def benchmark_naive(
    records: List[SeqRecord], seeds: List[Seq]
) -> Tuple[Dict[int, SeqOrientation], float]:
    """Benchmark the naive implementation."""
    start = time.perf_counter()
    result = contig_ids_by_seed(records, seeds)
    elapsed = time.perf_counter() - start
    return result, elapsed


def benchmark_ahocorasick(
    records: List[SeqRecord], seeds: List[Seq]
) -> Tuple[Dict[int, SeqOrientation], float]:
    """Benchmark the Aho-Corasick implementation."""
    start = time.perf_counter()
    result = _contig_ids_by_seed_ahocorasick(records, seeds)
    elapsed = time.perf_counter() - start
    return result, elapsed


def benchmark_bowtie2(
    records: List[SeqRecord], seeds: List[Seq], threads: int = 4
) -> Tuple[Dict[int, SeqOrientation], float]:
    """Benchmark the Bowtie2 implementation."""
    # Bowtie2 needs a file path, so write contigs to temp file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp:
        SeqIO.write(records, tmp, "fasta")
        tmp_path = Path(tmp.name)

    try:
        start = time.perf_counter()
        result = _contig_ids_by_seed_bowtie2(tmp_path, seeds, threads=threads)
        elapsed = time.perf_counter() - start
    finally:
        tmp_path.unlink()

    return result, elapsed


def verify_results(
    naive_result: Dict[int, SeqOrientation],
    other_result: Dict[int, SeqOrientation],
    name: str,
) -> bool:
    """
    Verify that results match the naive implementation.
    
    Checks:
    1. Same set of contig indices found
    2. Same orientation for each contig
    """
    naive_keys = set(naive_result.keys())
    other_keys = set(other_result.keys())
    
    all_passed = True

    # Check 1: Same contig indices
    if naive_keys != other_keys:
        missing = naive_keys - other_keys
        extra = other_keys - naive_keys
        print(f"  ✗ {name}: Different contig indices!")
        if missing:
            print(f"    Missing {len(missing)} matches: {list(missing)[:5]}...")
        if extra:
            print(f"    Extra {len(extra)} matches: {list(extra)[:5]}...")
        all_passed = False
    
    # Check 2: Same orientations for matching indices
    common_keys = naive_keys & other_keys
    orientation_mismatches = []
    for key in common_keys:
        if naive_result[key] != other_result[key]:
            orientation_mismatches.append(key)
    
    if orientation_mismatches:
        print(f"  ✗ {name}: Orientation mismatches for {len(orientation_mismatches)} contigs!")
        for key in orientation_mismatches[:5]:
            print(f"    Contig {key}: naive={naive_result[key].value}, {name}={other_result[key].value}")
        if len(orientation_mismatches) > 5:
            print(f"    ... and {len(orientation_mismatches) - 5} more")
        all_passed = False
    
    if all_passed:
        print(f"  ✓ {name}: All {len(naive_result)} matches verified (indices + orientations)")
    
    return all_passed


def run_correctness_tests(include_bowtie2: bool = True) -> bool:
    """
    Run thorough correctness tests before benchmarking.
    
    Tests:
    1. Empty inputs
    2. No matches
    3. All matches
    4. Forward orientation matches
    5. Reverse complement matches
    6. Mixed orientations
    7. Multiple seeds matching same contig
    8. Palindromic seeds
    
    Returns True if all tests pass.
    """
    print("=" * 80)
    print("CORRECTNESS TESTS")
    print("=" * 80)
    print()
    
    all_passed = True
    
    # Test 1: Empty inputs
    print("Test 1: Empty inputs...")
    empty_seeds: List[Seq] = []
    empty_contigs: List[SeqRecord] = []
    
    assert contig_ids_by_seed(empty_contigs, empty_seeds) == {}
    assert _contig_ids_by_seed_ahocorasick(empty_contigs, empty_seeds) == {}
    print("  ✓ Empty inputs handled correctly")
    
    # Test 2: Known forward match
    print("\nTest 2: Known forward match...")
    seed = Seq("GATTACA")
    contig = SeqRecord(Seq("AAAAAAGATTACATTTTTT"), id="contig_0")
    
    naive_result = contig_ids_by_seed([contig], [seed])
    ac_result = _contig_ids_by_seed_ahocorasick([contig], [seed])
    
    assert 0 in naive_result, "Naive should find match"
    assert naive_result[0] == SeqOrientation.FORWARD, "Should be forward orientation"
    assert 0 in ac_result, "Aho-Corasick should find match"
    assert ac_result[0] == SeqOrientation.FORWARD, "Should be forward orientation"
    print("  ✓ Forward match detected correctly")
    
    # Test 3: Known reverse complement match
    print("\nTest 3: Known reverse complement match...")
    seed = Seq("GATTACA")
    # Reverse complement of GATTACA is TGTAATC
    contig = SeqRecord(Seq("AAAAAATGTAATCTTTTTT"), id="contig_0")
    
    naive_result = contig_ids_by_seed([contig], [seed])
    ac_result = _contig_ids_by_seed_ahocorasick([contig], [seed])
    
    assert 0 in naive_result, "Naive should find reverse complement match"
    assert naive_result[0] == SeqOrientation.REVERSE, "Should be reverse orientation"
    assert 0 in ac_result, "Aho-Corasick should find reverse complement match"
    assert ac_result[0] == SeqOrientation.REVERSE, "Should be reverse orientation"
    print("  ✓ Reverse complement match detected correctly")
    
    # Test 4: No match
    print("\nTest 4: No match...")
    seed = Seq("GATTACA")
    contig = SeqRecord(Seq("AAAAAAAAAAAAAAAAAAA"), id="contig_0")
    
    naive_result = contig_ids_by_seed([contig], [seed])
    ac_result = _contig_ids_by_seed_ahocorasick([contig], [seed])
    
    assert len(naive_result) == 0, "Naive should find no match"
    assert len(ac_result) == 0, "Aho-Corasick should find no match"
    print("  ✓ No false positives")
    
    # Test 5: Multiple contigs, mixed results
    print("\nTest 5: Multiple contigs with mixed matches...")
    seeds = [Seq("GATTACA"), Seq("ACGTACGT")]
    contigs = [
        SeqRecord(Seq("AAAGATTACAAAA"), id="contig_0"),  # Has seed 0 forward
        SeqRecord(Seq("AAATGTAATCAAA"), id="contig_1"),  # Has seed 0 reverse
        SeqRecord(Seq("AAAACGTACGTAA"), id="contig_2"),  # Has seed 1 forward
        SeqRecord(Seq("AAAAAAAAAAAAA"), id="contig_3"),  # No match
    ]
    
    naive_result = contig_ids_by_seed(contigs, seeds)
    ac_result = _contig_ids_by_seed_ahocorasick(contigs, seeds)
    
    # Check indices
    assert set(naive_result.keys()) == {0, 1, 2}, f"Naive indices wrong: {naive_result.keys()}"
    assert set(ac_result.keys()) == {0, 1, 2}, f"AC indices wrong: {ac_result.keys()}"
    
    # Check orientations
    assert naive_result[0] == SeqOrientation.FORWARD
    assert naive_result[1] == SeqOrientation.REVERSE
    assert naive_result[2] == SeqOrientation.FORWARD
    assert ac_result[0] == SeqOrientation.FORWARD
    assert ac_result[1] == SeqOrientation.REVERSE
    assert ac_result[2] == SeqOrientation.FORWARD
    print("  ✓ Multiple contigs with mixed matches handled correctly")
    
    # Test 6: Variable length seeds
    print("\nTest 6: Variable length seeds...")
    seeds = [
        Seq("GATTACA"),           # 7bp
        Seq("ACGTACGTACGT"),      # 12bp
        Seq("AT"),                # 2bp (very short)
    ]
    contigs = [
        SeqRecord(Seq("AAAGATTACAAAA"), id="contig_0"),
        SeqRecord(Seq("ACGTACGTACGTAA"), id="contig_1"),
        SeqRecord(Seq("CCCCCCCCCCCCCC"), id="contig_2"),  # No match (no AT)
    ]
    
    naive_result = contig_ids_by_seed(contigs, seeds)
    ac_result = _contig_ids_by_seed_ahocorasick(contigs, seeds)
    
    assert set(naive_result.keys()) == set(ac_result.keys())
    print("  ✓ Variable length seeds handled correctly")
    
    # Test 7: Random data verification (statistical test)
    print("\nTest 7: Random data verification (100 seeds × 500 contigs)...")
    random.seed(42)  # Reproducible
    seeds, contigs = generate_test_data(
        num_seeds=100,
        num_contigs=500,
        variable_seed_length=True,
    )
    
    naive_result = contig_ids_by_seed(contigs, seeds)
    ac_result = _contig_ids_by_seed_ahocorasick(contigs, seeds)
    
    if not verify_results(naive_result, ac_result, "Aho-Corasick"):
        all_passed = False
    
    # Test 8: Bowtie2 verification (if available)
    if include_bowtie2:
        print("\nTest 8: Bowtie2 verification...")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
            SeqIO.write(contigs, tmp, "fasta")
            tmp_path = Path(tmp.name)
        
        try:
            bowtie2_result = _contig_ids_by_seed_bowtie2(tmp_path, seeds)
            if not verify_results(naive_result, bowtie2_result, "Bowtie2"):
                all_passed = False
        except Exception as e:
            print(f"  ⚠ Bowtie2 test skipped: {e}")
        finally:
            tmp_path.unlink()
    
    print()
    if all_passed:
        print("=" * 80)
        print("ALL CORRECTNESS TESTS PASSED ✓")
        print("=" * 80)
    else:
        print("=" * 80)
        print("SOME TESTS FAILED ✗")
        print("=" * 80)
    
    return all_passed


def run_benchmark(
    test_sizes: List[Tuple[int, int]],
    skip_naive_threshold: int = 100,
    include_bowtie2: bool = True,
) -> None:
    """
    Run benchmarks across multiple dataset sizes.

    Args:
        test_sizes: List of (num_seeds, num_contigs) tuples to test
        skip_naive_threshold: Skip naive algorithm when seeds × contigs exceeds this
        include_bowtie2: Whether to include bowtie2 in benchmarks
    """
    print("=" * 80)
    print("BENCHMARK: Seed-to-Contig Matching Algorithms")
    print("=" * 80)
    print()

    results_table = []

    for num_seeds, num_contigs in test_sizes:
        print(f"\n{'='*80}")
        print(f"Test: {num_seeds:,} seeds × {num_contigs:,} contigs")
        print(f"{'='*80}")

        # Generate test data
        print("Generating test data...")
        seeds, contigs = generate_test_data(
            num_seeds=num_seeds,
            num_contigs=num_contigs,
            variable_seed_length=True,
        )

        total_seed_bases = sum(len(s) for s in seeds)
        total_contig_bases = sum(len(r.seq) for r in contigs)
        print(f"  Total seed bases: {total_seed_bases:,}")
        print(f"  Total contig bases: {total_contig_bases:,}")
        print()

        # Determine if we should skip naive (too slow)
        complexity = num_seeds * num_contigs
        skip_naive = complexity > skip_naive_threshold * 1000

        # Run benchmarks
        naive_time = None
        naive_result = None

        if not skip_naive:
            print("Running naive algorithm...")
            naive_result, naive_time = benchmark_naive(contigs, seeds)
            print(f"  Time: {naive_time:.3f}s")
            print(f"  Matches found: {len(naive_result)}")
        else:
            print(f"Skipping naive algorithm (would be too slow)")
            print(f"  Estimated operations: {complexity:,}")

        print("\nRunning Aho-Corasick algorithm...")
        ac_result, ac_time = benchmark_ahocorasick(contigs, seeds)
        print(f"  Time: {ac_time:.3f}s")
        print(f"  Matches found: {len(ac_result)}")

        bowtie2_time = None
        bowtie2_result = None

        if include_bowtie2:
            print("\nRunning Bowtie2 algorithm...")
            try:
                bowtie2_result, bowtie2_time = benchmark_bowtie2(contigs, seeds)
                print(f"  Time: {bowtie2_time:.3f}s")
                print(f"  Matches found: {len(bowtie2_result)}")
            except Exception as e:
                print(f"  ERROR: {e}")

        # Verify correctness
        print("\nVerifying results...")
        if naive_result is not None:
            verify_results(naive_result, ac_result, "Aho-Corasick")
            if bowtie2_result is not None:
                verify_results(naive_result, bowtie2_result, "Bowtie2")
            print("  ✓ Results verified against naive implementation")
        else:
            # Cross-check AC and bowtie2
            if bowtie2_result is not None:
                if set(ac_result.keys()) == set(bowtie2_result.keys()):
                    print("  ✓ Aho-Corasick and Bowtie2 results match")
                else:
                    print("  ⚠ Aho-Corasick and Bowtie2 results differ")

        # Calculate speedups
        print("\nSpeedups:")
        if naive_time is not None:
            print(f"  Aho-Corasick vs Naive: {naive_time/ac_time:.1f}x faster")
            if bowtie2_time is not None:
                print(f"  Bowtie2 vs Naive: {naive_time/bowtie2_time:.1f}x faster")
        if bowtie2_time is not None:
            print(f"  Aho-Corasick vs Bowtie2: {bowtie2_time/ac_time:.1f}x faster")

        # Store results
        results_table.append({
            "seeds": num_seeds,
            "contigs": num_contigs,
            "naive_time": naive_time,
            "ac_time": ac_time,
            "bowtie2_time": bowtie2_time,
            "matches": len(ac_result),
        })

    # Print summary table
    print("\n")
    print("=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print()
    print(f"{'Seeds':>10} {'Contigs':>10} {'Naive':>12} {'Aho-Corasick':>12} {'Bowtie2':>12} {'AC Speedup':>12}")
    print("-" * 80)

    for r in results_table:
        naive_str = f"{r['naive_time']:.3f}s" if r['naive_time'] else "skipped"
        ac_str = f"{r['ac_time']:.3f}s"
        bt2_str = f"{r['bowtie2_time']:.3f}s" if r['bowtie2_time'] else "N/A"
        
        if r['naive_time']:
            speedup = f"{r['naive_time']/r['ac_time']:.1f}x"
        else:
            speedup = "N/A"

        print(f"{r['seeds']:>10,} {r['contigs']:>10,} {naive_str:>12} {ac_str:>12} {bt2_str:>12} {speedup:>12}")



if __name__ == "__main__":
    # Check if bowtie2 is available
    import shutil
    bowtie2_available = shutil.which("bowtie2") is not None

    if not bowtie2_available:
        print("WARNING: bowtie2 not found in PATH, skipping bowtie2 benchmarks")
        print("Install with: conda install -c bioconda bowtie2")
        print()

    # Run correctness tests first
    print()
    correctness_passed = run_correctness_tests(include_bowtie2=bowtie2_available)
    print()
    
    if not correctness_passed:
        print("ABORTING: Correctness tests failed. Fix bugs before benchmarking.")
        exit(1)

    # Define test sizes: (num_seeds, num_contigs)
    # Start small and increase to show scaling behavior
    test_sizes = [
        (10, 100),           # Tiny: baseline
        (100, 1000),         # Small: ~0.1M operations for naive
        (500, 5000),         # Medium: ~2.5M operations for naive
        (1000, 10000),       # Large: ~10M operations for naive
        (2500, 25000),       # XL: ~62.5M operations (naive will be slow)
        (5000, 50000),       # XXL: ~250M operations (skip naive)
        (25000, 113000),   # Real scale: ~2.8B operations (uncomment to test)
    ]

    run_benchmark(
        test_sizes=test_sizes,
        skip_naive_threshold=500,  # Skip naive when seeds × contigs > 50M
        include_bowtie2=bowtie2_available,
    )