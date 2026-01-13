#!/usr/bin/env python3
"""
Pileup Visualization Module

Visualizes how read pairs attest to assembled contigs

Usage:
    python pileup_visualization.py reads_R1.fastq reads_R2.fastq contigs.fasta contig_name seed_seq output.png
"""

import os
import subprocess
import tempfile
import numpy as np
import pysam
from dataclasses import dataclass
from collections import defaultdict
from typing import List, Tuple, Dict, Set, Optional
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class PileupReadPairRow:
    """
    Represents a single row in the pileup visualization (one read pair).

    Attributes:
        colors: (N, 4) RGBA array for the visual span of the read pair
        start_coord: Contig coordinate where this row begins (can be negative)
        ins_upstream: Number of insertions in the upstream mate
        ins_downstream: Number of insertions in the downstream mate
        contains_seed: Whether this read pair overlaps the seed region
    """
    colors: np.ndarray
    start_coord: int
    ins_upstream: int = 0
    ins_downstream: int = 0
    contains_seed: bool = False


# Default color palette
DEFAULT_COLORS = {
    "match": np.array([74, 144, 226, 255], dtype=np.uint8),       # #4A90E2 - Blue (seed-containing reads)
    "match_noseed": np.array([140, 170, 210, 255], dtype=np.uint8), # Muted blue (non-seed reads)
    "seed": np.array([255, 140, 0, 255], dtype=np.uint8),         # #FF8C00 - Orange
    "mismatch": np.array([80, 80, 80, 255], dtype=np.uint8),      # #505050 - Dark grey
    "deletion": np.array([0, 0, 0, 255], dtype=np.uint8),         # Black
    "soft_clip": np.array([100, 200, 220, 255], dtype=np.uint8),  # Cyan
    "insertion_marker": np.array([255, 0, 255, 255], dtype=np.uint8), # Magenta
    "insert_gap": np.array([224, 224, 224, 255], dtype=np.uint8), # #E0E0E0 - Light grey
    "bg": np.array([255, 255, 255, 255], dtype=np.uint8),         # #FFFFFF - White
}


# =============================================================================
# Alignment
# =============================================================================

def run_bowtie2(
    fwd_fastq: str,
    rev_fastq: str,
    contigs_fasta: str,
    output_dir: Optional[str] = None,
    threads: int = 4
) -> str:
    """
    Run bowtie2 alignment of paired reads against contigs.
    
    Args:
        fwd_fastq: Path to forward reads FASTQ
        rev_fastq: Path to reverse reads FASTQ
        contigs_fasta: Path to contigs FASTA
        output_dir: Directory for output files (default: same as contigs)
        threads: Number of threads for bowtie2
        
    Returns:
        Path to output BAM file (sorted and indexed)
    """
    if output_dir is None:
        output_dir = os.path.dirname(contigs_fasta) or "."
    
    base_name = Path(contigs_fasta).stem
    index_base = os.path.join(output_dir, f"{base_name}_bt2")
    sam_path = os.path.join(output_dir, f"{base_name}_aligned.sam")
    bam_path = os.path.join(output_dir, f"{base_name}_aligned.bam")
    sorted_bam_path = os.path.join(output_dir, f"{base_name}_aligned.sorted.bam")
    
    # Build index if needed
    if not os.path.exists(f"{index_base}.1.bt2"):
        print(f"Building bowtie2 index...")
        subprocess.run([
            "bowtie2-build", "-q", contigs_fasta, index_base
        ], check=True)
    
    # Run alignment if needed
    if not os.path.exists(sorted_bam_path):
        print(f"Running bowtie2 alignment...")
        subprocess.run([
            "bowtie2",
            "-x", index_base,
            "-1", fwd_fastq,
            "-2", rev_fastq,
            "--local",
            "--threads", str(threads),
            "-S", sam_path
        ], check=True)
        
        # Convert to BAM
        print("Converting to BAM...")
        subprocess.run([
            "samtools", "view", "-bS", "-o", bam_path, sam_path
        ], check=True)
        
        # Sort BAM
        print("Sorting BAM...")
        subprocess.run([
            "samtools", "sort", "-o", sorted_bam_path, bam_path
        ], check=True)
        
        # Index BAM
        print("Indexing BAM...")
        subprocess.run([
            "samtools", "index", sorted_bam_path
        ], check=True)
        
        # Clean up intermediates
        os.remove(sam_path)
        os.remove(bam_path)
    
    return sorted_bam_path


# =============================================================================
# Sequence Utilities
# =============================================================================

def load_contig(contigs_fasta: str, contig_name: str) -> str:
    """Load a specific contig sequence from a FASTA file."""
    with open(contigs_fasta) as f:
        for name, seq in SimpleFastaParser(f):
            # Handle cases where name might have description after space
            if name == contig_name or name.split()[0] == contig_name:
                return seq.upper()
    raise ValueError(f"Contig '{contig_name}' not found in {contigs_fasta}")


def find_seed_positions(contig_seq: str, seed_seq: str) -> Set[int]:
    """
    Find all positions in contig covered by the seed sequence.
    
    Searches for both the seed and its reverse complement.
    
    Args:
        contig_seq: The contig sequence
        seed_seq: The seed sequence to find
        
    Returns:
        Set of all contig positions (0-indexed) covered by seed occurrences
    """
    positions = set()
    contig_seq = contig_seq.upper()
    seed_seq = seed_seq.upper()
    seed_rc = reverse_complement(seed_seq)
    
    # Search for forward seed
    start = 0
    while True:
        idx = contig_seq.find(seed_seq, start)
        if idx == -1:
            break
        for pos in range(idx, idx + len(seed_seq)):
            positions.add(pos)
        start = idx + 1
    
    # Search for reverse complement (if different)
    if seed_rc != seed_seq:
        start = 0
        while True:
            idx = contig_seq.find(seed_rc, start)
            if idx == -1:
                break
            for pos in range(idx, idx + len(seed_rc)):
                positions.add(pos)
            start = idx + 1
    
    if not positions:
        print(f"Warning: Seed sequence not found in contig")
    
    return positions


# =============================================================================
# CIGAR Processing
# =============================================================================

# CIGAR operation codes (from SAM spec)
CIGAR_M = 0   # Match/mismatch (consumes query and reference)
CIGAR_I = 1   # Insertion (consumes query only)
CIGAR_D = 2   # Deletion (consumes reference only)
CIGAR_S = 4   # Soft clip (consumes query only)
CIGAR_H = 5   # Hard clip (consumes neither)
CIGAR_EQ = 7  # Sequence match (consumes query and reference)
CIGAR_X = 8   # Sequence mismatch (consumes query and reference)


def process_mate(
    read: pysam.AlignedSegment,
    contig_seq: str,
    seed_positions: Set[int],
    colors: Dict[str, np.ndarray]
) -> Tuple[List[np.ndarray], int, int, int, bool]:
    """
    Convert one mate's alignment to a list of pixel colors.

    Args:
        read: pysam AlignedSegment
        contig_seq: The contig sequence string
        seed_positions: Set of positions that are part of the seed
        colors: Color palette dict

    Returns:
        (color_list, start_coord, end_coord, insertion_count, hit_seed)
        start_coord/end_coord are in contig coordinates (can be negative)
        hit_seed is True if this mate overlaps any seed position
    """
    color_list = []
    insertions = 0
    hit_seed = False

    query_pos = 0  # Position in read sequence
    ref_pos = read.reference_start  # Position on contig
    start_coord = ref_pos  # Will be adjusted for leading soft clip

    if read.cigartuples is None:
        return [], ref_pos, ref_pos, 0, False

    for op_idx, (op, length) in enumerate(read.cigartuples):

        if op == CIGAR_M:  # Match/mismatch - need to check actual bases
            for _ in range(length):
                read_base = read.query_sequence[query_pos].upper()
                if 0 <= ref_pos < len(contig_seq):
                    ref_base = contig_seq[ref_pos].upper()
                else:
                    ref_base = 'N'

                if ref_pos in seed_positions:
                    color_list.append(colors["seed"])
                    hit_seed = True
                elif read_base == ref_base:
                    color_list.append(colors["match"])
                else:
                    color_list.append(colors["mismatch"])

                query_pos += 1
                ref_pos += 1

        elif op == CIGAR_EQ:  # Sequence match (guaranteed match)
            for _ in range(length):
                if ref_pos in seed_positions:
                    color_list.append(colors["seed"])
                    hit_seed = True
                else:
                    color_list.append(colors["match"])
                query_pos += 1
                ref_pos += 1

        elif op == CIGAR_X:  # Sequence mismatch (guaranteed mismatch)
            color_list.extend([colors["mismatch"]] * length)
            query_pos += length
            ref_pos += length

        elif op == CIGAR_I:  # Insertion (bases in read, not in reference)
            insertions += length
            query_pos += length
            # Mark the insertion location by coloring the previous aligned base
            if color_list:
                color_list[-1] = colors["insertion_marker"]

        elif op == CIGAR_D:  # Deletion (bases in reference, not in read)
            color_list.extend([colors["deletion"]] * length)
            ref_pos += length

        elif op == CIGAR_S:  # Soft clip
            if op_idx == 0:
                # Leading soft clip: bases are BEFORE reference_start
                start_coord = ref_pos - length
                color_list = [colors["soft_clip"]] * length + color_list
            else:
                # Trailing soft clip: bases are AFTER aligned region
                color_list.extend([colors["soft_clip"]] * length)
            query_pos += length

        elif op == CIGAR_H:  # Hard clip
            pass  # Bases not in query_sequence

    end_coord = start_coord + len(color_list)
    return color_list, start_coord, end_coord, insertions, hit_seed


def process_read_pair(
    upstream: pysam.AlignedSegment,
    downstream: pysam.AlignedSegment,
    contig_seq: str,
    seed_positions: Set[int],
    colors: Dict[str, np.ndarray]
) -> Optional[PileupReadPairRow]:
    """
    Combine two mates into a single PileupReadPairRow.

    Args:
        upstream: The mate with the lower reference_start
        downstream: The mate with the higher reference_start
        contig_seq: The contig sequence
        seed_positions: Set of seed positions
        colors: Color palette

    Returns:
        PileupReadPairRow or None if processing fails
    """
    up_colors, up_start, up_end, up_ins, up_hit_seed = process_mate(
        upstream, contig_seq, seed_positions, colors
    )
    down_colors, down_start, down_end, down_ins, down_hit_seed = process_mate(
        downstream, contig_seq, seed_positions, colors
    )

    contains_seed = up_hit_seed or down_hit_seed

    if not up_colors and not down_colors:
        return None

    # Helper to swap match -> match_noseed for non-seed-containing reads
    def adjust_colors_for_noseed(color_list):
        if contains_seed:
            return color_list
        # Replace match color with match_noseed
        match_color = colors["match"]
        noseed_color = colors["match_noseed"]
        return [noseed_color if np.array_equal(c, match_color) else c for c in color_list]

    up_colors = adjust_colors_for_noseed(up_colors)
    down_colors = adjust_colors_for_noseed(down_colors)

    # Handle single-mate case
    if not up_colors:
        return PileupReadPairRow(
            colors=np.array(down_colors, dtype=np.uint8),
            start_coord=down_start,
            ins_upstream=down_ins,
            ins_downstream=0,
            contains_seed=contains_seed
        )
    if not down_colors:
        return PileupReadPairRow(
            colors=np.array(up_colors, dtype=np.uint8),
            start_coord=up_start,
            ins_upstream=up_ins,
            ins_downstream=0,
            contains_seed=contains_seed
        )

    # Combine mates with gap or overlap handling
    gap = down_start - up_end

    if gap > 0:
        # Gap between mates - fill with insert_gap color
        gap_colors = [colors["insert_gap"]] * gap
        combined = up_colors + gap_colors + down_colors
    elif gap < 0:
        # Mates overlap - skip overlapping portion of downstream
        overlap = -gap
        if overlap < len(down_colors):
            combined = up_colors + down_colors[overlap:]
        else:
            # Downstream entirely within upstream
            combined = up_colors
    else:
        # Exactly adjacent
        combined = up_colors + down_colors

    return PileupReadPairRow(
        colors=np.array(combined, dtype=np.uint8),
        start_coord=up_start,
        ins_upstream=up_ins,
        ins_downstream=down_ins,
        contains_seed=contains_seed
    )


# =============================================================================
# BAM Parsing
# =============================================================================

def get_read_pair_rows(
    bam_path: str,
    contig_name: str,
    contig_seq: str,
    seed_positions: Set[int],
    colors: Dict[str, np.ndarray]
) -> List[PileupReadPairRow]:
    """
    Parse a BAM file to generate PileupReadPairRow objects for a contig.
    
    Args:
        bam_path: Path to sorted, indexed BAM file
        contig_name: Name of contig to visualize
        contig_seq: The contig sequence
        seed_positions: Set of seed positions
        colors: Color palette
        
    Returns:
        List of PileupReadPairRow objects, sorted by start position
    """
    samfile = pysam.AlignmentFile(bam_path, "rb")
    read_pairs = defaultdict(list)
    
    # Collect reads mapping to our contig
    try:
        reads = samfile.fetch(contig_name)
    except ValueError:
        # Contig not in BAM - try without index
        print(f"Warning: Could not fetch by region, scanning entire BAM")
        reads = (r for r in samfile if r.reference_name == contig_name)
    
    for read in reads:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        read_pairs[read.query_name].append(read)
    
    samfile.close()
    
    rows = []
    for qname, mates in read_pairs.items():
        if len(mates) == 2:
            # Sort by reference_start to determine upstream/downstream
            mates.sort(key=lambda x: x.reference_start)
            row = process_read_pair(mates[0], mates[1], contig_seq, seed_positions, colors)
        elif len(mates) == 1:
            # Single mate (other mate unmapped or mapped elsewhere)
            mate = mates[0]
            color_list, start, end, ins, hit_seed = process_mate(
                mate, contig_seq, seed_positions, colors
            )
            if color_list:
                row = PileupReadPairRow(
                    colors=np.array(color_list, dtype=np.uint8),
                    start_coord=start,
                    ins_upstream=ins,
                    ins_downstream=0,
                    contains_seed=hit_seed,
                )
            else:
                row = None
        else:
            continue
        
        if row is not None:
            rows.append(row)
    
    rows.sort(key=lambda r: r.start_coord)
    
    return rows


# =============================================================================
# Rendering
# =============================================================================

class PileupRenderer:
    """Renders PileupReadPairRow objects into a visualization."""
    
    def __init__(self, colors: Optional[Dict[str, np.ndarray]] = None):
        self.colors = DEFAULT_COLORS.copy()
        if colors:
            self.colors.update(colors)
    
    def render(
        self,
        contig_len: int,
        seed_positions: Set[int],
        rows: List[PileupReadPairRow],
        scale: Tuple[int, int] = (2, 2),
        padding: int = 10,
        contig_height: int = 1
    ) -> np.ndarray:
        """
        Render the pileup visualization.

        Args:
            contig_len: Length of the contig
            seed_positions: Set of seed positions (for coloring contig ribbon)
            rows: List of PileupReadPairRow objects
            scale: (x_scale, y_scale) for pixel scaling
            padding: Padding in pixels
            contig_height: Height of contig ribbon in pixels (default 3 for thin ribbon)

        Returns:
            RGBA numpy array of the rendered image
        """
        if not rows:
            print("Warning: No rows to render")
            return np.full((100, 100, 4), self.colors["bg"], dtype=np.uint8)

        # Sort rows by start_coord for classic pileup diagonal pattern
        rows = sorted(rows, key=lambda r: r.start_coord)

        # Calculate bounds
        min_x = min(0, min(r.start_coord for r in rows))
        max_x = max(contig_len, max(r.start_coord + len(r.colors) for r in rows))
        genome_width = max_x - min_x

        max_ins_up = max(r.ins_upstream for r in rows)
        max_ins_down = max(r.ins_downstream for r in rows)

        # Calculate canvas dimensions - just 1 pixel per insertion (magenta)
        ins_up_width = max_ins_up
        ins_down_width = max_ins_down
        pad_left = padding if max_ins_up > 0 else 0
        pad_right = padding if max_ins_down > 0 else 0

        total_width = ins_up_width + pad_left + genome_width + pad_right + ins_down_width
        total_height = padding + contig_height + padding + len(rows) + padding

        # Create canvas
        canvas = np.full((total_height, total_width, 4), self.colors["bg"], dtype=np.uint8)

        # Genome section offset
        genome_x_offset = ins_up_width + pad_left

        # Draw contig ribbon
        contig_y_start = padding
        contig_y_end = padding + contig_height
        contig_x_in_genome = 0 - min_x  # Where contig position 0 is in genome section

        # Draw contig ribbon - first fill with match color, then overlay seed positions
        contig_x_start = genome_x_offset + contig_x_in_genome
        canvas[contig_y_start:contig_y_end, contig_x_start:contig_x_start + contig_len] = self.colors["match"]
        for x in seed_positions:
            canvas_x = contig_x_start + x
            if 0 <= canvas_x < total_width:
                canvas[contig_y_start:contig_y_end, canvas_x] = self.colors["seed"]

        # Draw read rows
        reads_y_start = padding + contig_height + padding

        for i, row in enumerate(rows):
            y = reads_y_start + i

            # Draw main alignment using numpy slicing
            x_start = genome_x_offset + (row.start_coord - min_x)
            x_end = x_start + len(row.colors)
            # Clip to canvas bounds
            src_start = max(0, -x_start)
            src_end = len(row.colors) - max(0, x_end - total_width)
            dst_start = max(0, x_start)
            dst_end = min(total_width, x_end)
            if dst_start < dst_end:
                canvas[y, dst_start:dst_end] = row.colors[src_start:src_end]

            # Draw upstream insertions (magenta pixels)
            if row.ins_upstream > 0:
                ins_count = min(row.ins_upstream, ins_up_width)
                canvas[y, 0:ins_count] = self.colors["insertion_marker"]

            # Draw downstream insertions (magenta pixels)
            if row.ins_downstream > 0:
                ins_count = min(row.ins_downstream, ins_down_width)
                down_x_start = total_width - ins_down_width
                canvas[y, down_x_start:down_x_start + ins_count] = self.colors["insertion_marker"]

        # Scale image
        scaled = np.repeat(canvas, scale[0], axis=1)
        scaled = np.repeat(scaled, scale[1], axis=0)

        return scaled
    
    def add_legend(self, image: np.ndarray) -> np.ndarray:
        """
        Add a color legend below the visualization.

        Args:
            image: The rendered pileup image (numpy array)

        Returns:
            New numpy array with legend appended at the bottom
        """
        from PIL import Image as PILImage

        img_height, img_width = image.shape[:2]

        legend_items = [
            ("Seed read", self.colors["match"]),
            ("Non-seed read", self.colors["match_noseed"]),
            ("Seed bases", self.colors["seed"]),
            ("Mismatch", self.colors["mismatch"]),
            ("Deletion", self.colors["deletion"]),
            ("Soft Clip", self.colors["soft_clip"]),
            ("Insertion", self.colors["insertion_marker"]),
            ("Unsequenced", self.colors["insert_gap"]),
        ]

        # Layout parameters
        box_size = 12
        item_spacing = 100
        row_height = 20
        margin = 15

        # Calculate how many items per row
        items_per_row = max(1, (img_width - 2 * margin) // item_spacing)
        num_rows = (len(legend_items) + items_per_row - 1) // items_per_row
        legend_height = num_rows * row_height + 2 * margin

        # Create new canvas
        new_height = img_height + legend_height
        canvas = np.full((new_height, img_width, 4), self.colors["bg"], dtype=np.uint8)

        # Copy original image
        canvas[:img_height, :, :] = image

        # Draw legend using PIL for text
        pil_image = PILImage.fromarray(canvas)
        draw = ImageDraw.Draw(pil_image)

        for idx, (label, color) in enumerate(legend_items):
            row = idx // items_per_row
            col = idx % items_per_row

            x = margin + col * item_spacing
            y = img_height + margin + row * row_height

            # Draw color box with border
            draw.rectangle(
                [x, y, x + box_size, y + box_size],
                fill=tuple(color),
                outline=(150, 150, 150, 255)
            )
            # Draw label
            draw.text((x + box_size + 5, y - 1), label, fill=(60, 60, 60, 255))

        return np.array(pil_image)

    def save(self, image: np.ndarray, path: str):
        """Save numpy array image to file."""
        from PIL import Image as PILImage
        PILImage.fromarray(image).save(path)
        print(f"Saved: {path}")


# =============================================================================
# Main Entry Point
# =============================================================================

def create_pileup_visualization(
    fwd_fastq: str,
    rev_fastq: str,
    contigs_fasta: str,
    contig_name: str,
    seed_sequence: str,
    output_path: Optional[str] = None,
    colormap: Optional[Dict] = None,
    scale: Tuple[int, int] = (2, 2),
    padding: int = 10,
    contig_height: int = 10,
    threads: int = 4,
    work_dir: Optional[str] = None
) -> np.ndarray:
    """
    Create a pileup visualization showing read pair alignments to a contig.
    
    Args:
        fwd_fastq: Path to forward reads FASTQ file
        rev_fastq: Path to reverse reads FASTQ file
        contigs_fasta: Path to assembled contigs FASTA file
        contig_name: Name of the specific contig to visualize
        seed_sequence: Seed sequence to highlight in visualization
        output_path: Optional path to save PNG image
        colormap: Optional color overrides (dict of name -> RGBA array)
        scale: (x, y) scaling factors for output image
        padding: Padding in pixels around image elements
        contig_height: Height of contig ribbon in pixels
        threads: Number of threads for bowtie2
        work_dir: Working directory for intermediate files
        
    Returns:
        Numpy array of shape (height, width, 4) with RGBA values
    """
    # Setup colors
    colors = DEFAULT_COLORS.copy()
    if colormap:
        for k, v in colormap.items():
            colors[k] = np.array(v, dtype=np.uint8)
    
    # Setup working directory
    if work_dir is None:
        work_dir = os.path.dirname(contigs_fasta) or "."
    
    # Step 1: Run alignment
    print(f"Step 1: Aligning reads to contigs...")
    bam_path = run_bowtie2(fwd_fastq, rev_fastq, contigs_fasta, work_dir, threads)
    
    # Step 2: Load contig sequence
    print(f"Step 2: Loading contig '{contig_name}'...")
    contig_seq = load_contig(contigs_fasta, contig_name)
    print(f"  Contig length: {len(contig_seq)} bp")
    
    # Step 3: Find seed positions
    print(f"Step 3: Finding seed positions...")
    seed_positions = find_seed_positions(contig_seq, seed_sequence)
    print(f"  Found {len(seed_positions)} positions covered by seed")
    
    # Step 4: Parse alignments
    print(f"Step 4: Parsing alignments...")
    rows = get_read_pair_rows(bam_path, contig_name, contig_seq, seed_positions, colors)
    print(f"  Found {len(rows)} read pairs")
    
    # Step 5: Rendering
    print(f"Step 5: Rendering...")
    renderer = PileupRenderer(colors)
    image = renderer.render(
        contig_len=len(contig_seq),
        seed_positions=seed_positions,
        rows=rows,
        scale=scale,
        padding=padding,
        contig_height=contig_height
    )
    print(f"  Image size: {image.shape[1]} x {image.shape[0]} pixels")
    
    # Step 6: Save with legend
    if output_path:
        image_with_legend = renderer.add_legend(image)
        renderer.save(image_with_legend, output_path)

    return image


# =============================================================================
# CLI
# =============================================================================

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Create pileup visualization of read alignments to a contig"
    )
    parser.add_argument("fwd_fastq", help="Forward reads FASTQ file")
    parser.add_argument("rev_fastq", help="Reverse reads FASTQ file")
    parser.add_argument("contigs_fasta", help="Contigs FASTA file")
    parser.add_argument("contig_name", help="Name of contig to visualize")
    parser.add_argument("seed_sequence", help="Seed sequence to highlight")
    parser.add_argument("output", help="Output PNG path")
    parser.add_argument("--scale-x", type=int, default=2, help="Horizontal scale factor")
    parser.add_argument("--scale-y", type=int, default=2, help="Vertical scale factor")
    parser.add_argument("--padding", type=int, default=10, help="Padding in pixels")
    parser.add_argument("--contig-height", type=int, default=10, help="Contig ribbon height")
    parser.add_argument("--threads", type=int, default=4, help="Threads for bowtie2")
    
    args = parser.parse_args()
    
    create_pileup_visualization(
        fwd_fastq=args.fwd_fastq,
        rev_fastq=args.rev_fastq,
        contigs_fasta=args.contigs_fasta,
        contig_name=args.contig_name,
        seed_sequence=args.seed_sequence,
        output_path=args.output,
        scale=(args.scale_x, args.scale_y),
        padding=args.padding,
        contig_height=args.contig_height,
        threads=args.threads
    )


if __name__ == "__main__":
    main()