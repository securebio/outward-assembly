import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional
from PIL import Image

@dataclass
class PileupReadPairRow:
    """
    Represent a single row in the pileup.
    colors: (N, 4) RGBA array representing the visual span of the read pair.
    start_coord: The contig-relative index where this row begins.
    ins_upstream: Number of insertion events to render on the left.
    ins_downstream: Number of insertion events to render on the right.
    """
    colors: np.ndarray  
    start_coord: int
    ins_upstream: int = 0
    ins_downstream: int = 0

class PileupRenderer:
    def __init__(self, colormap: Optional[Dict] = None):
        # Default palette from your specifications
        self.colors = {
            "match": [74, 144, 226, 255],        # Blue
            "seed": [255, 140, 0, 255],         # Orange
            "mismatch": [80, 80, 80, 255],     # Dark Grey
            "insert_gap": [224, 224, 224, 255], # Light Grey
            "bg": [255, 255, 255, 255],         # White
            "ins_marker_a": [0, 0, 0, 255],     # Black for indel pattern
            "ins_marker_b": [255, 255, 255, 255] # White for indel pattern
        }
        if colormap:
            self.colors.update(colormap)

    def _create_insertion_pattern(self, count: int, height: int) -> np.ndarray:
        """Creates the <black pixel><white pixel> vertical indicator."""
        pattern = np.zeros((height, 2 * count, 4), dtype=np.uint8)
        for i in range(count):
            pattern[:, 2*i] = self.colors["ins_marker_a"]
            pattern[:, 2*i + 1] = self.colors["ins_marker_b"]
        return pattern

    def render_pileup(
        self, 
        contig_len: int, 
        seed_range: tuple, 
        rows: List[PileupReadPairRow],
        scale: tuple = (2, 5),
        padding: int = 10,
        contig_height: int = 10
    ) -> np.ndarray:
        """
        Assembles the rows into a final image.
        scale: (x_scale, y_scale) for the 'intentionally dumb' scaling.
        """
        seed_start, seed_end = seed_range
        
        # 1. Calculate horizontal bounds
        # We need to account for reads that might start before 0 or end after contig_len
        min_x = min(0, min((r.start_coord for r in rows), default=0))
        max_x = max(contig_len, max((r.start_coord + len(r.colors) for r in rows), default=contig_len))
        genome_width = max_x - min_x
        
        # Calculate max insertions for sidebars
        max_ins_up = max((r.ins_upstream for r in rows), default=0)
        max_ins_down = max((r.ins_downstream for r in rows), default=0)
        
        # 2. Setup Image Canvas
        # Width = LeftIns + Pad + Genome + Pad + RightIns
        total_width = (max_ins_up * 2) + (padding if max_ins_up > 0 else 0) + \
                      genome_width + \
                      (padding if max_ins_down > 0 else 0) + (max_ins_down * 2)
        
        total_height = padding + contig_height + padding + len(rows) + padding
        
        canvas = np.full((total_height, total_width, 4), self.colors["bg"], dtype=np.uint8)
        
        # Offset for the genome section
        genome_offset_x = (max_ins_up * 2) + (padding if max_ins_up > 0 else 0)
        
        # 3. Draw Contig Ribbon
        contig_y_start = padding
        contig_y_end = padding + contig_height
        # Calculate where contig 0 is relative to min_x
        contig_x_start = genome_offset_x + (0 - min_x)
        
        canvas[contig_y_start:contig_y_end, contig_x_start:contig_x_start + contig_len] = self.colors["match"]
        canvas[contig_y_start:contig_y_end, contig_x_start + seed_start:contig_x_start + seed_end] = self.colors["seed"]

        # 4. Draw Read Rows
        for i, row in enumerate(rows):
            y = padding + contig_height + padding + i
            
            # Position the read colors relative to the genome min_x
            x_start = genome_offset_x + (row.start_coord - min_x)
            x_end = x_start + len(row.colors)
            canvas[y, x_start:x_end] = row.colors
            
            # Draw insertions if present
            if row.ins_upstream > 0:
                pattern = self._create_insertion_pattern(row.ins_upstream, 1)
                canvas[y, 0:row.ins_upstream*2] = pattern
                
            if row.ins_downstream > 0:
                pattern = self._create_insertion_pattern(row.ins_downstream, 1)
                canvas[y, -row.ins_downstream*2:] = pattern

        # 5. "Intentionally Dumb" Scaling
        # Blow up each pixel to a patch of (scale_y, scale_x)
        scaled_img = np.repeat(canvas, scale[0], axis=1) # Scale width
        scaled_img = np.repeat(scaled_img, scale[1], axis=0) # Scale height
        
        return scaled_img

    def save(self, array: np.ndarray, filename: str):
        img = Image.fromarray(array)
        img.save(filename)