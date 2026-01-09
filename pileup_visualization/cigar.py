import pysam
import numpy as np
from collections import defaultdict
from typing import List, Tuple, Dict, Optional

def get_read_pair_rows(
    bam_path: str, 
    contig_name: str,
    contig_seq: str,  # Need this to check for mismatches
    seed_positions: set,  # Set of positions, more flexible than range
    colormap: Dict
) -> List['PileupReadPairRow']:
    """Parse BAM file to generate PileupReadPairRow objects for a contig."""
    samfile = pysam.AlignmentFile(bam_path, "rb")
    read_pairs = defaultdict(list)

    for read in samfile.fetch(contig_name):
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            read_pairs[read.query_name].append(read)

    rows = []
    for qname, mates in read_pairs.items():
        if len(mates) != 2:
            continue
        
        mates.sort(key=lambda x: x.reference_start)
        upstream, downstream = mates[0], mates[1]

        row = process_pair(upstream, downstream, contig_seq, seed_positions, colormap)
        if row is not None:
            rows.append(row)

    samfile.close()
    return rows


def process_mate(
    read: pysam.AlignedSegment,
    contig_seq: str,
    seed_positions: set,
    colormap: Dict
) -> Tuple[List, int, int, int]:
    """
    Convert one mate's alignment to colors.
    
    Returns:
        (colors, start_coord, end_coord, insertion_count)
        
    start_coord/end_coord are in contig coordinates (can be negative for leading soft clips)
    """
    colors = []
    insertions = 0
    
    query_pos = 0  # Position in read sequence
    ref_pos = read.reference_start  # Position on contig
    
    # Track actual start (adjusted for leading soft clip)
    start_coord = ref_pos
    
    for op_idx, (op, length) in enumerate(read.cigartuples):
        
        if op == 0:  # M - alignment match (could be match or mismatch)
            for _ in range(length):
                read_base = read.query_sequence[query_pos]
                ref_base = contig_seq[ref_pos] if 0 <= ref_pos < len(contig_seq) else 'N'
                
                if ref_pos in seed_positions:
                    colors.append(colormap["seed"])
                elif read_base.upper() == ref_base.upper():
                    colors.append(colormap["match"])
                else:
                    colors.append(colormap["mismatch"])
                
                query_pos += 1
                ref_pos += 1
        
        elif op == 7:  # = - sequence match
            for _ in range(length):
                if ref_pos in seed_positions:
                    colors.append(colormap["seed"])
                else:
                    colors.append(colormap["match"])
                query_pos += 1
                ref_pos += 1
        
        elif op == 8:  # X - sequence mismatch
            for _ in range(length):
                colors.append(colormap["mismatch"])
                query_pos += 1
                ref_pos += 1
        
        elif op == 1:  # I - insertion (bases in read, not in reference)
            insertions += length
            query_pos += length
            # No color added, no ref_pos change
        
        elif op == 2:  # D - deletion (bases in reference, not in read)
            for _ in range(length):
                colors.append(colormap["mismatch"])
                ref_pos += 1
            # No query_pos change
        
        elif op == 4:  # S - soft clip
            if op_idx == 0:
                # Leading soft clip: bases are BEFORE reference_start
                start_coord = ref_pos - length
                for _ in range(length):
                    colors.append(colormap["mismatch"])
                # Note: we prepend conceptually, but since we process left-to-right
                # and adjusted start_coord, the colors list order is correct
            else:
                # Trailing soft clip: bases are AFTER aligned region
                for _ in range(length):
                    colors.append(colormap["mismatch"])
            query_pos += length
            # ref_pos doesn't change for soft clips
        
        elif op == 5:  # H - hard clip
            pass  # Bases not in query_sequence, nothing to do
    
    end_coord = start_coord + len(colors)
    return colors, start_coord, end_coord, insertions


def process_pair(
    upstream: pysam.AlignedSegment,
    downstream: pysam.AlignedSegment,
    contig_seq: str,
    seed_positions: set,
    colormap: Dict
) -> Optional['PileupReadPairRow']:
    """Combine two mates into a single PileupReadPairRow."""
    
    up_colors, up_start, up_end, up_ins = process_mate(
        upstream, contig_seq, seed_positions, colormap
    )
    down_colors, down_start, down_end, down_ins = process_mate(
        downstream, contig_seq, seed_positions, colormap
    )
    
    # Handle gap or overlap between mates
    gap = down_start - up_end
    
    if gap > 0:
        # There's a gap - fill with insert_gap color
        gap_colors = [colormap["insert_gap"]] * gap
        combined = up_colors + gap_colors + down_colors
    elif gap < 0:
        # Mates overlap - skip redundant portion of downstream
        overlap = -gap
        combined = up_colors + down_colors[overlap:]
    else:
        # Exactly adjacent
        combined = up_colors + down_colors
    
    return PileupReadPairRow(
        colors=np.array(combined, dtype=np.uint8),
        start_coord=up_start,
        ins_upstream=up_ins,
        ins_downstream=down_ins
    )