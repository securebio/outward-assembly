#!/usr/bin/env python3
"""
Generate synthetic test data and visualize it.

This creates fake alignments that exercise different CIGAR operations
and visualizes them without needing real sequencing data.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Set
from pileup_visualization import (
    PileupReadPairRow,
    PileupRenderer,
    find_seed_positions,
    DEFAULT_COLORS,
)


@dataclass
class MockRead:
    """Simulates a pysam.AlignedSegment for testing."""
    query_name: str
    reference_start: int
    query_sequence: str
    cigartuples: List[Tuple[int, int]]
    is_read1: bool = True
    is_unmapped: bool = False
    is_secondary: bool = False
    is_supplementary: bool = False


def process_mock_mate(
    read: MockRead,
    contig_seq: str,
    seed_positions: Set[int],
    colors: dict
) -> Tuple[List[np.ndarray], int, int, int, bool]:
    """
    Process a mock read - same logic as the real function.
    Returns (color_list, start_coord, end_coord, insertions, hit_seed)
    """
    color_list = []
    insertions = 0
    hit_seed = False

    query_pos = 0
    ref_pos = read.reference_start
    start_coord = ref_pos

    if read.cigartuples is None:
        return [], ref_pos, ref_pos, 0, False

    for op_idx, (op, length) in enumerate(read.cigartuples):

        if op == 0:  # M - alignment match
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

        elif op == 7:  # = - sequence match
            for _ in range(length):
                if ref_pos in seed_positions:
                    color_list.append(colors["seed"])
                    hit_seed = True
                else:
                    color_list.append(colors["match"])
                query_pos += 1
                ref_pos += 1

        elif op == 8:  # X - sequence mismatch
            for _ in range(length):
                color_list.append(colors["mismatch"])
                query_pos += 1
                ref_pos += 1

        elif op == 1:  # I - insertion
            insertions += length
            query_pos += length
            if color_list:
                color_list[-1] = colors["insertion_marker"]

        elif op == 2:  # D - deletion
            for _ in range(length):
                color_list.append(colors["deletion"])
                ref_pos += 1

        elif op == 4:  # S - soft clip
            if op_idx == 0:
                start_coord = ref_pos - length
                soft_colors = [colors["soft_clip"]] * length
                color_list = soft_colors + color_list
            else:
                for _ in range(length):
                    color_list.append(colors["soft_clip"])
            query_pos += length

    end_coord = start_coord + len(color_list)
    return color_list, start_coord, end_coord, insertions, hit_seed


def create_test_scenarios():
    """
    Create various alignment scenarios to test visualization.
    
    Returns contig sequence, seed, and list of mock read pairs.
    """
    # Contig: 100bp with seed in the middle (positions 40-59)
    contig_seq = "A" * 40 + "TCGATCGATCGATCGATCGA" + "A" * 40  # 100bp total
    seed_seq = "TCGATCGATCGATCGATCGA"  # 20bp seed
    
    scenarios = []
    
    # Scenario 1: Perfect alignment pair with gap
    # Read1: positions 10-59, Read2: positions 70-99
    read1 = MockRead(
        query_name="perfect_pair",
        reference_start=10,
        query_sequence="A" * 30 + "TCGATCGATCGATCGATCGA",  # 50bp
        cigartuples=[(0, 50)],  # 50M
        is_read1=True
    )
    read2 = MockRead(
        query_name="perfect_pair",
        reference_start=70,
        query_sequence="A" * 30,  # 30bp
        cigartuples=[(0, 30)],  # 30M
        is_read1=False
    )
    scenarios.append(("Perfect pair with gap", read1, read2))
    
    # Scenario 2: Pair with mismatches
    read1 = MockRead(
        query_name="mismatch_pair",
        reference_start=5,
        query_sequence="A" * 10 + "GGG" + "A" * 27 + "TCGATCGATCG",  # Mismatches at 15-17
        cigartuples=[(0, 51)],
        is_read1=True
    )
    read2 = MockRead(
        query_name="mismatch_pair",
        reference_start=65,
        query_sequence="CCC" + "A" * 32,  # Mismatches at start
        cigartuples=[(0, 35)],
        is_read1=False
    )
    scenarios.append(("Pair with mismatches", read1, read2))
    
    # Scenario 3: Read with insertion
    read1 = MockRead(
        query_name="insertion_pair",
        reference_start=20,
        query_sequence="A" * 20 + "NNNNN" + "TCGATCGATCGATCGATCGA",  # 5bp insertion
        cigartuples=[(0, 20), (1, 5), (0, 20)],  # 20M5I20M
        is_read1=True
    )
    read2 = MockRead(
        query_name="insertion_pair",
        reference_start=75,
        query_sequence="A" * 25,
        cigartuples=[(0, 25)],
        is_read1=False
    )
    scenarios.append(("Pair with insertion", read1, read2))
    
    # Scenario 4: Read with deletion
    read1 = MockRead(
        query_name="deletion_pair",
        reference_start=15,
        query_sequence="A" * 20 + "TCGATCGATCGATCGATCGA",  # 40bp but covers 43bp of ref
        cigartuples=[(0, 20), (2, 3), (0, 20)],  # 20M3D20M
        is_read1=True
    )
    read2 = MockRead(
        query_name="deletion_pair",
        reference_start=70,
        query_sequence="A" * 30,
        cigartuples=[(0, 30)],
        is_read1=False
    )
    scenarios.append(("Pair with deletion", read1, read2))
    
    # Scenario 5: Soft-clipped read (extends before contig start)
    read1 = MockRead(
        query_name="softclip_start",
        reference_start=0,  # Aligns at start, but has soft clip before
        query_sequence="GGGGG" + "A" * 45,  # 5bp soft clip + 45bp aligned
        cigartuples=[(4, 5), (0, 45)],  # 5S45M
        is_read1=True
    )
    read2 = MockRead(
        query_name="softclip_start",
        reference_start=60,
        query_sequence="A" * 35,
        cigartuples=[(0, 35)],
        is_read1=False
    )
    scenarios.append(("Soft-clip at start (negative coords)", read1, read2))
    
    # Scenario 6: Soft-clipped read (extends past contig end)
    read1 = MockRead(
        query_name="softclip_end",
        reference_start=30,
        query_sequence="A" * 10 + "TCGATCGATCGATCGATCGA" + "A" * 20,
        cigartuples=[(0, 50)],
        is_read1=True
    )
    read2 = MockRead(
        query_name="softclip_end",
        reference_start=90,  # Near end
        query_sequence="A" * 10 + "GGGGG",  # 10bp aligned + 5bp soft clip past end
        cigartuples=[(0, 10), (4, 5)],  # 10M5S
        is_read1=False
    )
    scenarios.append(("Soft-clip past contig end", read1, read2))
    
    # Scenario 7: Overlapping mates (short insert)
    read1 = MockRead(
        query_name="overlap_pair",
        reference_start=35,
        query_sequence="A" * 5 + "TCGATCGATCGATCGATCGA" + "A" * 15,  # 40bp
        cigartuples=[(0, 40)],
        is_read1=True
    )
    read2 = MockRead(
        query_name="overlap_pair",
        reference_start=55,  # Overlaps with read1 (read1 ends at 75)
        query_sequence="ATCGATCGA" + "A" * 31,  # 40bp
        cigartuples=[(0, 40)],
        is_read1=False
    )
    scenarios.append(("Overlapping mates", read1, read2))
    
    # Scenario 8: Complex CIGAR (soft clip + match + insertion + match + deletion + match)
    read1 = MockRead(
        query_name="complex_cigar",
        reference_start=25,
        query_sequence="GG" + "A" * 15 + "NNN" + "TCGATCGATCGA" + "A" * 8,
        cigartuples=[(4, 2), (0, 15), (1, 3), (0, 12), (2, 2), (0, 8)],  # 2S15M3I12M2D8M
        is_read1=True
    )
    read2 = MockRead(
        query_name="complex_cigar",
        reference_start=75,
        query_sequence="A" * 25,
        cigartuples=[(0, 25)],
        is_read1=False
    )
    scenarios.append(("Complex CIGAR", read1, read2))

    # Scenario 9: Reads that DON'T contain seed (should be muted blue)
    # Positions 0-35 (before seed at 40-59)
    read1 = MockRead(
        query_name="no_seed_pair",
        reference_start=0,
        query_sequence="A" * 20,
        cigartuples=[(0, 20)],
        is_read1=True
    )
    read2 = MockRead(
        query_name="no_seed_pair",
        reference_start=25,
        query_sequence="A" * 10,
        cigartuples=[(0, 10)],
        is_read1=False
    )
    scenarios.append(("Non-seed pair (muted)", read1, read2))

    return contig_seq, seed_seq, scenarios


def process_pair_to_row(
    read1: MockRead,
    read2: MockRead,
    contig_seq: str,
    seed_positions: Set[int],
    colors: dict
) -> PileupReadPairRow:
    """Process a read pair into a PileupReadPairRow."""

    # Sort by start position
    if read1.reference_start <= read2.reference_start:
        upstream, downstream = read1, read2
    else:
        upstream, downstream = read2, read1

    up_colors, up_start, up_end, up_ins, up_hit = process_mock_mate(
        upstream, contig_seq, seed_positions, colors
    )
    down_colors, down_start, down_end, down_ins, down_hit = process_mock_mate(
        downstream, contig_seq, seed_positions, colors
    )

    contains_seed = up_hit or down_hit

    # Swap match -> match_noseed for non-seed-containing reads
    def adjust_colors(color_list):
        if contains_seed:
            return color_list
        match_color = colors["match"]
        noseed_color = colors["match_noseed"]
        return [noseed_color if np.array_equal(c, match_color) else c for c in color_list]

    up_colors = adjust_colors(up_colors)
    down_colors = adjust_colors(down_colors)

    # Handle gap or overlap
    gap = down_start - up_end

    if gap > 0:
        gap_colors = [colors["insert_gap"]] * gap
        combined = up_colors + gap_colors + down_colors
    elif gap < 0:
        overlap = -gap
        if overlap < len(down_colors):
            combined = up_colors + down_colors[overlap:]
        else:
            combined = up_colors
    else:
        combined = up_colors + down_colors

    return PileupReadPairRow(
        colors=np.array(combined, dtype=np.uint8),
        start_coord=up_start,
        ins_upstream=up_ins,
        ins_downstream=down_ins,
        contains_seed=contains_seed
    )


def create_labeled_visualization():
    """Create visualization with labels for each scenario."""
    from PIL import Image, ImageDraw

    contig_seq, seed_seq, scenarios = create_test_scenarios()
    seed_positions = find_seed_positions(contig_seq, seed_seq)

    print(f"Contig length: {len(contig_seq)}")
    print(f"Seed positions: {min(seed_positions)}-{max(seed_positions)}")
    print(f"Number of scenarios: {len(scenarios)}")
    print()

    # Process each scenario
    rows = []
    labels = []
    for name, read1, read2 in scenarios:
        print(f"Processing: {name}")
        print(f"  Read1: start={read1.reference_start}, CIGAR={read1.cigartuples}")
        print(f"  Read2: start={read2.reference_start}, CIGAR={read2.cigartuples}")

        row = process_pair_to_row(read1, read2, contig_seq, seed_positions, DEFAULT_COLORS)
        rows.append(row)
        labels.append(name)

        print(f"  Result: start_coord={row.start_coord}, width={len(row.colors)}, ins_up={row.ins_upstream}, ins_down={row.ins_downstream}")
        print()

    # Render with uniform scaling
    scale = 4
    padding = 12
    contig_height = 12
    row_height = 6  # Height per read pair row

    renderer = PileupRenderer()
    image = renderer.render(
        contig_len=len(contig_seq),
        seed_positions=seed_positions,
        rows=rows,
        scale=(scale, scale),  # Uniform scaling
        padding=padding,
        contig_height=contig_height
    )

    # Layout constants
    label_margin = 20
    label_area_width = 220
    axis_height = 30
    legend_row_height = 25
    legend_margin = 15

    # Calculate final dimensions
    img_width = image.shape[1]
    img_height = image.shape[0]

    legend_items = [
        ("Match", DEFAULT_COLORS["match"]),
        ("Seed", DEFAULT_COLORS["seed"]),
        ("Mismatch", DEFAULT_COLORS["mismatch"]),
        ("Deletion", DEFAULT_COLORS["deletion"]),
        ("Soft Clip", DEFAULT_COLORS["soft_clip"]),
        ("Insertion", DEFAULT_COLORS["insertion_marker"]),
        ("Unsequenced", DEFAULT_COLORS["insert_gap"]),
    ]

    # Legend layout: 4 items per row
    items_per_row = 4
    legend_rows = (len(legend_items) + items_per_row - 1) // items_per_row
    legend_height = legend_rows * legend_row_height + legend_margin * 2

    final_width = img_width + label_area_width
    final_height = axis_height + img_height + legend_height

    # Create final canvas
    final_image = Image.new('RGBA', (final_width, final_height), (255, 255, 255, 255))

    # Paste the rendered image
    pil_image = Image.fromarray(image)
    final_image.paste(pil_image, (0, axis_height))

    draw = ImageDraw.Draw(final_image)

    # --- Draw coordinate axis ---
    # Calculate genome offset (same logic as renderer)
    min_x = min(0, min(r.start_coord for r in rows))
    max_ins_up = max(r.ins_upstream for r in rows)
    ins_up_width = max_ins_up * 2
    pad_left = padding if max_ins_up > 0 else 0
    genome_x_offset = ins_up_width + pad_left

    def coord_to_x(coord):
        canvas_x = genome_x_offset + (coord - min_x)
        return canvas_x * scale

    axis_y = axis_height - 8

    # Draw axis line
    start_x = coord_to_x(min_x)
    end_x = coord_to_x(len(contig_seq))
    draw.line([(start_x, axis_y), (end_x, axis_y)], fill=(100, 100, 100, 255), width=1)

    # Draw ticks every 25bp
    for coord in range(0, len(contig_seq) + 1, 25):
        x = coord_to_x(coord)
        draw.line([(x, axis_y), (x, axis_y - 4)], fill=(100, 100, 100, 255), width=1)
        text = str(coord)
        text_width = len(text) * 6
        draw.text((x - text_width // 2, axis_y - 18), text, fill=(80, 80, 80, 255))

    # --- Draw row labels ---
    row_y_start = axis_height + (padding + contig_height + padding) * scale

    for i, label in enumerate(labels):
        y = row_y_start + i * scale + scale // 2 - 4
        draw.text((img_width + label_margin, y), label, fill=(60, 60, 60, 255))

    # --- Draw legend ---
    legend_y_start = axis_height + img_height + legend_margin
    box_size = 14
    item_width = final_width // items_per_row

    for idx, (label, color) in enumerate(legend_items):
        row = idx // items_per_row
        col = idx % items_per_row

        x = col * item_width + 20
        y = legend_y_start + row * legend_row_height

        # Draw color box with border
        draw.rectangle([x, y, x + box_size, y + box_size], fill=tuple(color), outline=(180, 180, 180, 255))
        draw.text((x + box_size + 6, y), label, fill=(60, 60, 60, 255))

    # Save
    output_path = "synthetic_test_visualization.png"
    final_image.save(output_path)
    print(f"Saved: {output_path}")

    return final_image


def create_simple_visualization():
    """Create a clean visualization with legend but no row labels."""
    contig_seq, seed_seq, scenarios = create_test_scenarios()
    seed_positions = find_seed_positions(contig_seq, seed_seq)

    rows = []
    for name, read1, read2 in scenarios:
        row = process_pair_to_row(read1, read2, contig_seq, seed_positions, DEFAULT_COLORS)
        rows.append(row)

    scale = 4
    padding = 5

    renderer = PileupRenderer()
    image = renderer.render(
        contig_len=len(contig_seq),
        seed_positions=seed_positions,
        rows=rows,
        scale=(scale, scale),
        padding=padding
    )

    # Add legend and save
    image_with_legend = renderer.add_legend(image)
    output_path = "synthetic_test_simple.png"
    renderer.save(image_with_legend, output_path)

    return image


def create_large_realistic_dataset(
    contig_len: int = 500,
    seed_start: int = 200,
    seed_len: int = 30,
    num_reads: int = 200,
    read_len: int = 100,
    insert_size_mean: int = 250,
    insert_size_std: int = 50,
    mismatch_rate: float = 0.01,
    seed_coverage_boost: float = 3.0
):
    """
    Generate a large synthetic dataset similar to real sequencing data.

    Creates a contig with a seed region, and generates read pairs with:
    - Realistic insert size distribution
    - Higher coverage around the seed region
    - Occasional mismatches
    - Some reads that don't overlap the seed
    """
    import random

    # Create contig sequence
    bases = ['A', 'C', 'G', 'T']
    contig_seq = ''.join(random.choice(bases) for _ in range(contig_len))

    # Insert a specific seed sequence
    seed_seq = ''.join(random.choice(bases) for _ in range(seed_len))
    contig_seq = contig_seq[:seed_start] + seed_seq + contig_seq[seed_start + seed_len:]

    seed_positions = set(range(seed_start, seed_start + seed_len))

    scenarios = []

    for i in range(num_reads):
        # Determine read start position
        # Boost coverage near seed region
        if random.random() < 0.6:  # 60% of reads near seed
            # Center around seed with some spread
            center = seed_start + seed_len // 2
            start = int(random.gauss(center - insert_size_mean // 2, insert_size_std))
        else:
            # Random position across contig
            start = random.randint(-read_len // 2, contig_len - read_len // 2)

        start = max(-10, min(start, contig_len - 50))  # Clamp to reasonable range

        # Generate insert size
        insert_size = int(random.gauss(insert_size_mean, insert_size_std))
        insert_size = max(read_len + 20, min(insert_size, contig_len))

        # Read 1 sequence (may have mismatches)
        read1_seq = ""
        for j in range(read_len):
            pos = start + j
            if 0 <= pos < contig_len:
                base = contig_seq[pos]
                if random.random() < mismatch_rate:
                    base = random.choice([b for b in bases if b != base])
                read1_seq += base
            else:
                read1_seq += random.choice(bases)

        # Read 2 starts at start + insert_size - read_len
        read2_start = start + insert_size - read_len
        read2_seq = ""
        for j in range(read_len):
            pos = read2_start + j
            if 0 <= pos < contig_len:
                base = contig_seq[pos]
                if random.random() < mismatch_rate:
                    base = random.choice([b for b in bases if b != base])
                read2_seq += base
            else:
                read2_seq += random.choice(bases)

        # Create mock reads with simple CIGAR (all matches)
        # Add occasional soft clips for reads near edges
        cigar1 = [(0, read_len)]
        cigar2 = [(0, read_len)]

        # Add soft clips for reads that extend past contig
        if start < 0:
            clip_len = -start
            cigar1 = [(4, clip_len), (0, read_len - clip_len)]
            start = 0

        if read2_start + read_len > contig_len:
            clip_len = read2_start + read_len - contig_len
            cigar2 = [(0, read_len - clip_len), (4, clip_len)]

        read1 = MockRead(
            query_name=f"read_{i}",
            reference_start=max(0, start),
            query_sequence=read1_seq,
            cigartuples=cigar1,
            is_read1=True
        )

        read2 = MockRead(
            query_name=f"read_{i}",
            reference_start=max(0, read2_start),
            query_sequence=read2_seq,
            cigartuples=cigar2,
            is_read1=False
        )

        scenarios.append((f"read_{i}", read1, read2))

    return contig_seq, seed_seq, scenarios, seed_positions


def create_realistic_pileup():
    """Create a realistic-looking pileup visualization with many reads."""
    print("Generating large realistic dataset...")

    contig_seq, seed_seq, scenarios, seed_positions = create_large_realistic_dataset(
        contig_len=500,
        seed_start=200,
        seed_len=30,
        num_reads=150,
        read_len=80,
        insert_size_mean=200,
        insert_size_std=40
    )

    print(f"  Contig length: {len(contig_seq)}")
    print(f"  Seed positions: {min(seed_positions)}-{max(seed_positions)}")
    print(f"  Number of read pairs: {len(scenarios)}")

    rows = []
    for name, read1, read2 in scenarios:
        row = process_pair_to_row(read1, read2, contig_seq, seed_positions, DEFAULT_COLORS)
        rows.append(row)

    print(f"  Rows with seed: {sum(1 for r in rows if r.contains_seed)}")
    print(f"  Rows without seed: {sum(1 for r in rows if not r.contains_seed)}")

    renderer = PileupRenderer()
    image = renderer.render(
        contig_len=len(contig_seq),
        seed_positions=seed_positions,
        rows=rows,
        scale=(2, 2),
        padding=5
    )

    image_with_legend = renderer.add_legend(image)
    output_path = "synthetic_realistic_pileup.png"
    renderer.save(image_with_legend, output_path)

    return image


if __name__ == "__main__":
    print("=" * 60)
    print("Creating synthetic test visualization")
    print("=" * 60)
    print()

    create_labeled_visualization()
    print()
    create_simple_visualization()
    print()
    create_realistic_pileup()