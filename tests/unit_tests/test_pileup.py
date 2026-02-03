""" pileup smoke tests """

import numpy as np
import pytest

from outward_assembly.pileup import (
    DEFAULT_COLORS,
    PileupRenderer,
    create_pileup_visualization,
    find_seed_positions,
)


class TestPileupImports:
    """Test that the pileup module imports correctly."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_import_main_function(self):
        assert callable(create_pileup_visualization)

    @pytest.mark.fast
    @pytest.mark.unit
    def test_import_renderer(self):
        assert PileupRenderer is not None

    @pytest.mark.fast
    @pytest.mark.unit
    def test_import_colors(self):
        expected_colors = [
            "match_and_contains_seed",
            "match_no_seed",
            "seed",
            "seed_mismatch",
            "mismatch",
            "deletion",
            "soft_clip",
            "insertion_marker",
            "insertion_stripe",
            "unsequenced",
            "background",
        ]
        for color_name in expected_colors:
            assert color_name in DEFAULT_COLORS, f"Missing color: {color_name}"
            assert DEFAULT_COLORS[color_name].shape == (
                4,
            ), f"Color {color_name} wrong shape"
            assert DEFAULT_COLORS[color_name].dtype == np.uint8


class TestFindSeedPositions:
    """Test seed position finding."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_find_seed_forward(self):
        contig = "AAAAAATCGATCGAAAAAA"
        seed = "ATCGATCG"
        positions = find_seed_positions(contig, seed)

        # Seed "ATCGATCG" starts at position 5 (the 6th A begins the match)
        # Spans positions 5-12 (8 bases)
        assert positions == set(range(5, 13))

    @pytest.mark.fast
    @pytest.mark.unit
    def test_find_seed_reverse_complement(self):
        contig = "AAAAAACGATCGATAAAAAA"  # Contains CGATCGAT, RC of ATCGATCG
        seed = "ATCGATCG"
        positions = find_seed_positions(contig, seed)

        # RC is at positions 6-13
        assert positions == set(range(6, 14))

    @pytest.mark.fast
    @pytest.mark.unit
    def test_seed_not_found(self, capsys):
        contig = "AAAAAAAAAAAAAAAA"
        seed = "GCGCGCGC"
        positions = find_seed_positions(contig, seed)

        assert positions == set()
        captured = capsys.readouterr()
        assert "Warning" in captured.out


class TestPileupRenderer:
    """Test the renderer class."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_renderer_initialization(self):
        renderer = PileupRenderer()
        assert renderer.colors is not None
        assert "background" in renderer.colors

    @pytest.mark.fast
    @pytest.mark.unit
    def test_renderer_empty_rows(self):
        renderer = PileupRenderer()
        image = renderer.render(
            contig_len=100,
            seed_positions=set(range(40, 60)),
            rows=[],
            scale=(1, 1),
            padding=5,
            contig_height=2,
        )

        assert image.ndim == 3
        assert image.shape[2] == 4  # RGBA
        assert image.dtype == np.uint8
