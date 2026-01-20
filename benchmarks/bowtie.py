"""
Optimized _subset_contigs using bowtie2 for seed-to-contig matching.

This replaces the O(seeds × contigs) nested loop with:
1. O(n) index building (bowtie2-build)
2. O(seeds × log(index_size)) alignment (bowtie2)
3. O(matches) SAM parsing

For 25k seeds × 113k contigs, this reduces runtime from ~2 hours to ~seconds.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from outward_assembly.basic_seq_operations import SeqOrientation


# Constants (matching the original code)
MEGAHIT_OUT_PREFIX = "megahit_out_iter"
MEGAHIT_FINAL_CONTIGS = "final.contigs.fa"
MEGAHIT_FILTERED_CONTIGS = "contigs_filtered.fasta"


def contig_ids_by_seed_bowtie2(
    contigs_path: Path,
    seed_seqs: List[Seq],
    threads: int = 4,
) -> Dict[int, SeqOrientation]:
    """
    Find which contigs contain which seeds using bowtie2.
    
    This is the core optimization: instead of O(seeds × contigs) substring
    searches in Python, we use bowtie2's FM-index for O(seeds × log(n)) lookups.
    
    Args:
        contigs_path: Path to the contigs FASTA file
        seed_seqs: List of seed sequences to search for
        threads: Number of threads for bowtie2
        
    Returns:
        Dict mapping contig index to orientation (FORWARD or REVERSE)
        
    How it works:
    -------------
    1. BUILD INDEX: bowtie2-build creates an FM-index of all contigs.
       - FM-index is a compressed suffix array that enables fast substring search
       - Building is O(n) where n = total contig sequence length
       - This is the "expensive upfront work" that pays off during search
       
    2. ALIGN SEEDS: bowtie2 aligns each seed against the index.
       - Each seed query is O(m) where m = seed length (independent of contig count!)
       - We use --end-to-end to require the full seed to match
       - We use --score-min C,0,0 to require perfect matches (no mismatches)
       - The -a flag reports ALL alignments (a seed might match multiple contigs)
       
    3. PARSE SAM: Extract contig indices and orientations from alignment output.
       - SAM flag 16 indicates reverse complement alignment
       - We build a name→index mapping to convert contig names back to indices
    """
    
    if not seed_seqs:
        return {}
    
    # Build contig name → index mapping
    # We need this because bowtie2 reports contig names, not indices
    records: List[SeqRecord] = list(SeqIO.parse(contigs_path, "fasta"))
    contig_name_to_idx: Dict[str, int] = {
        rec.id: idx for idx, rec in enumerate(records)
    }
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # =====================================================================
        # STEP 1: Build bowtie2 index
        # =====================================================================
        # This creates several files: contigs_index.1.bt2, .2.bt2, .3.bt2, .4.bt2,
        # .rev.1.bt2, .rev.2.bt2
        # 
        # These contain the FM-index (Burrows-Wheeler transform + auxiliary data)
        # which allows O(m) substring queries regardless of reference size.
        
        index_prefix = tmpdir / "contigs_index"
        
        build_cmd = [
            "bowtie2-build",
            "--threads", str(threads),
            "--quiet",              # suppress progress messages
            str(contigs_path),
            str(index_prefix)
        ]
        
        result = subprocess.run(build_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bowtie2-build failed: {result.stderr}")
        
        # =====================================================================
        # STEP 2: Write seeds to temporary FASTA
        # =====================================================================
        # bowtie2 needs input as a file, not a Python list
        
        seeds_path = tmpdir / "seeds.fasta"
        seed_records = [
            SeqRecord(Seq(str(seq)), id=f"seed_{i}", description="")
            for i, seq in enumerate(seed_seqs)
        ]
        SeqIO.write(seed_records, seeds_path, "fasta")
        
        # =====================================================================
        # STEP 3: Run bowtie2 alignment
        # =====================================================================
        # Key flags explained:
        #
        # -x INDEX_PREFIX    : Use the index we just built
        # -f                 : Input is FASTA (not FASTQ)
        # -a                 : Report ALL alignments, not just the best
        #                      (a seed might appear in multiple contigs)
        # --end-to-end       : Require entire seed to align (no soft clipping)
        # --score-min C,0,0  : Minimum score = 0, meaning perfect match required
        #                      (any mismatch gives negative score, failing this threshold)
        # --no-unal          : Don't output unaligned seeds (reduces output size)
        # --no-hd            : No SAM header (we don't need it)
        # --no-sq            : No @SQ lines (we don't need them)
        # -p THREADS         : Use multiple threads for alignment
        
        sam_path = tmpdir / "alignments.sam"
        
        align_cmd = [
            "bowtie2",
            "-x", str(index_prefix),
            "-f", str(seeds_path),
            "-a",                       # all alignments
            "--end-to-end",             # full seed must match
            "--score-min", "C,0,0",     # perfect matches only
            "--no-unal",                # skip unaligned
            "--no-hd",                  # no header
            "--no-sq",                  # no @SQ lines
            "-p", str(threads),
            "-S", str(sam_path)
        ]
        
        result = subprocess.run(align_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bowtie2 failed: {result.stderr}")
        
        # =====================================================================
        # STEP 4: Parse SAM output
        # =====================================================================
        # SAM format (tab-separated):
        # Col 1: Query name (seed_0, seed_1, ...)
        # Col 2: FLAG (bitwise flags, 16 = reverse complement)
        # Col 3: Reference name (contig name)
        # Col 4: Position (1-based)
        # Col 5: MAPQ (mapping quality)
        # Col 6+: CIGAR, mate info, etc.
        #
        # We only need columns 2 (FLAG) and 3 (contig name).
        # FLAG & 16 tells us if the seed aligned to the reverse complement.
        
        matches: Dict[int, SeqOrientation] = {}
        
        with open(sam_path) as f:
            for line in f:
                # Skip any header lines that might have slipped through
                if line.startswith('@'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                # seed_name = fields[0]  # e.g., "seed_0" - we don't need this
                flag = int(fields[1])
                contig_name = fields[2]
                
                # '*' means no alignment (shouldn't happen with --no-unal, but be safe)
                if contig_name == '*':
                    continue
                
                # Look up the contig index
                if contig_name not in contig_name_to_idx:
                    continue  # shouldn't happen, but be defensive
                    
                contig_idx = contig_name_to_idx[contig_name]
                
                # FLAG bit 16 (0x10) = SEQ is reverse complemented
                # If this bit is set, the seed matched the reverse complement of the contig,
                # meaning the contig is in REVERSE orientation relative to the seed.
                orientation = SeqOrientation.REVERSE if (flag & 16) else SeqOrientation.FORWARD
                
                # Note: if a contig matches multiple seeds, we might overwrite.
                # The original code had the same behavior - last match wins.
                # If you need to handle multiple seeds per contig differently,
                # you'd need to modify this logic.
                matches[contig_idx] = orientation
        
        return matches


def _subset_contigs_bowtie2(
    workdir: Path,
    iter: int,
    seed_seqs: List[Seq],
    include_overlaps: bool = True,
    overlap_n0: int = 7,
    overlap_n1: int = 31,
    threads: int = 4,
) -> None:
    """
    Subset assembled contigs to those containing seed sequences, using bowtie2.
    
    This is a drop-in replacement for _subset_contigs that uses bowtie2 for
    the seed-to-contig matching step, reducing runtime from O(seeds × contigs)
    to O(seeds × log(contigs)).
    
    Args:
        workdir: Working directory containing megahit output
        iter: Current iteration number
        seed_seqs: Seed sequences to search for
        include_overlaps: Whether to include contigs connected via overlaps
        overlap_n0: Minimum overlap length for exact matches
        overlap_n1: Minimum overlap length when allowing 1 error
        threads: Number of threads for bowtie2
        
    Output:
        Creates workdir/megahit_out_iter<iter>-<subiter>/contigs_filtered.fasta
        containing contigs that contain seeds (or overlap with seed-containing contigs),
        oriented forward with respect to the seed.
    """
    workdir = Path(workdir)
    
    # Find all subiter directories for this iteration
    # e.g., megahit_out_iter1-0, megahit_out_iter1-1, etc.
    subiter_dirs = [
        d
        for d in workdir.iterdir()
        if d.is_dir() and d.name.startswith(f"{MEGAHIT_OUT_PREFIX}{iter}-")
    ]
    
    for subiter_dir in subiter_dirs:
        contigs_path = subiter_dir / MEGAHIT_FINAL_CONTIGS
        if not contigs_path.is_file():
            continue
            
        subset_path = subiter_dir / MEGAHIT_FILTERED_CONTIGS
        
        # Load all contigs (we need them for overlap detection and output)
        records: List[SeqRecord] = list(SeqIO.parse(contigs_path, "fasta"))
        
        if not records:
            # No contigs to process
            SeqIO.write([], subset_path, "fasta")
            continue
        
        # =====================================================================
        # CORE OPTIMIZATION: Use bowtie2 instead of nested loops
        # =====================================================================
        subsetted_ids_and_orientations: Dict[int, SeqOrientation] = contig_ids_by_seed_bowtie2(
            contigs_path, seed_seqs, threads=threads
        )
        
        # =====================================================================
        # Overlap expansion (unchanged from original)
        # =====================================================================
        # This finds contigs that don't contain a seed themselves, but are
        # connected via overlapping sequences to a contig that does.
        # This is important for assembly because ambiguities create branches.
        if include_overlaps:
            # Import the original overlap function
            # (This could also be optimized, but it's not the main bottleneck)
            from outward_assembly.basic_seq_operations import get_overlapping_sequence_ids
            
            seqs = [rec.seq for rec in records]
            subsetted_ids_and_orientations = get_overlapping_sequence_ids(
                seqs, subsetted_ids_and_orientations, overlap_n0, overlap_n1
            )
        
        # =====================================================================
        # Orient and write filtered contigs
        # =====================================================================
        filtered_records = []
        for idx, orientation in subsetted_ids_and_orientations.items():
            record = records[idx]
            if record.seq is not None and orientation == SeqOrientation.REVERSE:
                # Contig is reverse relative to seed - flip it so output is
                # consistently forward with respect to the seed
                record.seq = record.seq.reverse_complement()
            filtered_records.append(record)
        
        SeqIO.write(filtered_records, subset_path, "fasta")


# =============================================================================
# Example usage and testing
# =============================================================================

if __name__ == "__main__":
    import time
    
    # Quick sanity check
    print("Testing bowtie2-based contig subsetting...")
    
    # Check that bowtie2 is available
    result = subprocess.run(["bowtie2", "--version"], capture_output=True, text=True)
    if result.returncode != 0:
        print("ERROR: bowtie2 not found. Install with: conda install bowtie2")
        exit(1)
    
    version_line = result.stdout.split('\n')[0]
    print(f"Found: {version_line}")
    
    # Create test data
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create fake contigs
        contigs_path = tmpdir / "contigs.fasta"
        test_contigs = [
            SeqRecord(Seq("AAAAAAAAAA" + "GATTACA" + "TTTTTTTTTT"), id="contig_0"),
            SeqRecord(Seq("CCCCCCCCCC" + "GATTACA" + "GGGGGGGGGG"), id="contig_1"),
            SeqRecord(Seq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), id="contig_2"),  # no seed
        ]
        SeqIO.write(test_contigs, contigs_path, "fasta")
        
        # Test seeds
        seeds = [Seq("GATTACA")]
        
        # Run the optimized function
        start = time.time()
        result = contig_ids_by_seed_bowtie2(contigs_path, seeds)
        elapsed = time.time() - start
        
        print(f"\nResults (took {elapsed:.3f}s):")
        for idx, orientation in result.items():
            print(f"  Contig {idx}: {orientation.value}")
        
        # Verify
        assert 0 in result, "Should find seed in contig_0"
        assert 1 in result, "Should find seed in contig_1"
        assert 2 not in result, "Should NOT find seed in contig_2"
        
        print("\n✓ All tests passed!")