import csv
import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from outward_assembly.pipeline_steps import _subset_contigs


@pytest.fixture
def temp_workdir_with_contigs():
    """Create a temporary directory with mock contigs for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)

        # Create megahit output directory
        megahit_dir = workdir / "megahit_out_iter1-1"
        megahit_dir.mkdir(parents=True)

        # Create mock contigs
        contigs = [
            SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="contig1"),
            SeqRecord(Seq("TTTTTTTTTTTTTTTTT"), id="contig2"),
            SeqRecord(Seq("GGGGGGAAAAAACCCCCC"), id="contig3"),
            SeqRecord(Seq("TACGCGTACGCGTACGTCGTA"), id="contig4"),
            SeqRecord(Seq("ATAGACGATGGCGAGGCGCGTA"), id="contig5"),
        ]

        # Write contigs to file
        contigs_path = megahit_dir / "final.contigs.fa"
        SeqIO.write(contigs, contigs_path, "fasta")

        yield workdir


@pytest.fixture
def overlapping_contig(temp_workdir_with_contigs):
    """
    Adds an additional contig to temp_workdir_with_contigs. This new contig overlaps with
    contig 5
    """
    megahit_dir = temp_workdir_with_contigs / "megahit_out_iter1-1"
    contigs_path = megahit_dir / "final.contigs.fa"
    # Create an additional contig that overlaps with contig 5
    with open(contigs_path, "a") as f:
        f.write(">contig6\nGCGCGTAATATAAAGGCC\n")


@pytest.fixture
def overlapping_contig_reverse(temp_workdir_with_contigs):
    """
    Adds an additional contig to temp_workdir_with_contigs. This new contig overlaps with
    the reverse complement of contig 5
    """
    megahit_dir = temp_workdir_with_contigs / "megahit_out_iter1-1"
    contigs_path = megahit_dir / "final.contigs.fa"
    # Create an additional contig that overlaps with RC of contig 5
    with open(contigs_path, "a") as f:
        f.write(">contig7\nCGTCTATATATATATAT\n")


@pytest.mark.fast
@pytest.mark.unit
def test_subset_contigs_multiple_seeds(temp_workdir_with_contigs):
    """Test _subset_contigs with multiple seed sequences."""
    workdir = temp_workdir_with_contigs

    # Define multiple seeds
    seed_seqs = [
        Seq("ACGTAC"),  # seed1: in contig1, RC in contig4
        Seq("AAAAAA"),  # seed2: in contig3, RC in contig2
    ]

    # Run subset_contigs
    _subset_contigs(workdir, iter=1, seed_seqs=seed_seqs)

    # Check output
    filtered_path = workdir / "megahit_out_iter1-1" / "contigs_filtered.fasta"
    filtered_contigs = list(SeqIO.parse(filtered_path, "fasta"))

    # Should find 4 contigs (1, 2, 3, and 4)
    assert len(filtered_contigs) == 4

    # Check contig IDs
    contig_ids = [rec.id for rec in filtered_contigs]
    assert "contig1" in contig_ids
    assert "contig2" in contig_ids
    assert "contig3" in contig_ids
    assert "contig4" in contig_ids
    assert "contig5" not in contig_ids


@pytest.mark.fast
@pytest.mark.unit
def test_subset_contigs_with_overlaps(temp_workdir_with_contigs, overlapping_contig):
    """Test _subset_contigs with include_overlaps parameter."""
    workdir = temp_workdir_with_contigs
    megahit_dir = temp_workdir_with_contigs / "megahit_out_iter1-1"

    # Define seeds
    seed_seqs = [Seq("ATATAAAGGCC")]

    # Run subset_contigs with include_overlaps=True
    _subset_contigs(workdir, iter=1, seed_seqs=seed_seqs, include_overlaps=True)

    # Check output
    filtered_path = megahit_dir / "contigs_filtered.fasta"
    filtered_contigs = list(SeqIO.parse(filtered_path, "fasta"))

    # Should include contig5 due to overlap
    contig_ids = [rec.id for rec in filtered_contigs]
    assert len(contig_ids) == 2
    assert "contig5" in contig_ids

    graph_path = megahit_dir / "contig_overlap_graph.tsv"
    with open(graph_path) as f:
        overlap_edges = list(csv.DictReader(f, delimiter="\t"))

    assert overlap_edges == [
        {
            "source_contig": "contig5",
            "target_contig": "contig6",
            "source_oriented": "contig5",
            "target_oriented": "contig6",
            "edge_orientation": "forward",
            "source_contains_seed": "false",
            "target_contains_seed": "true",
        }
    ]


@pytest.mark.fast
@pytest.mark.unit
def test_subset_contigs_preserves_seed_orientation(temp_workdir_with_contigs):
    """
    Test that _subset_contigs returns contig sequences in the forward orientation with respect
    to the seed sequence
    """
    workdir = temp_workdir_with_contigs

    seed_seqs = [
        Seq("AAAAAA"),  # in contig3, RC in contig2
    ]

    # Run subset_contigs
    _subset_contigs(workdir, iter=1, seed_seqs=seed_seqs, include_overlaps=False)

    # Check output
    filtered_path = workdir / "megahit_out_iter1-1" / "contigs_filtered.fasta"
    filtered_contigs = list(SeqIO.parse(filtered_path, "fasta"))

    # Should find 2 contigs (forward in contig3 and reverse in contig2)
    assert len(filtered_contigs) == 2

    # Check that contigs were returned in a forward orientation with respect to the seed
    contig_sequences = {str(rec.seq) for rec in filtered_contigs}
    assert {
        "AAAAAAAAAAAAAAAAA",  # reverse complement of contig2
        "GGGGGGAAAAAACCCCCC",  # contig3
    } == contig_sequences

    graph_path = workdir / "megahit_out_iter1-1" / "contig_overlap_graph.tsv"
    assert not graph_path.exists()


@pytest.mark.fast
@pytest.mark.unit
def test_subset_contigs_preserves_seed_orientation_for_overlaps(
    temp_workdir_with_contigs, overlapping_contig, overlapping_contig_reverse
):
    """
    Test that when _subset_contigs returns overlapping contigs (that don't contain a seed
    themselves but overlap with a contig that does), the overlapping contigs are also in the
    forward direction with respect to the seed
    """
    workdir = temp_workdir_with_contigs

    # This seed is present in the RC of contig5. Contig5 overlaps with contig6, and the RC
    # of contig5 overlaps with contig6. So, once we orient all three contigs with respect
    # to the seed, we expect to see: contig5 in RC, contig6 in RC, and contig7 in forward.
    seed_seqs = [
        Seq("GGCCTTTATAT"),
    ]

    # Run subset_contigs
    _subset_contigs(workdir, iter=1, seed_seqs=seed_seqs, include_overlaps=True)

    # Check output
    filtered_path = workdir / "megahit_out_iter1-1" / "contigs_filtered.fasta"
    filtered_contigs = list(SeqIO.parse(filtered_path, "fasta"))

    # Should find 2 contigs (RC in contig5, contig5 overlaps in contig6)
    assert len(filtered_contigs) == 3

    # Check that contigs were returned in a forward orientation with respect to the seed
    contig_sequences = {str(rec.seq) for rec in filtered_contigs}
    assert {
        "TACGCGCCTCGCCATCGTCTAT",  # reverse compliment of contig5
        "GGCCTTTATATTACGCGC",  # reverse compliment of contig6
        "CGTCTATATATATATAT",  # contig7 (forward)
    } == contig_sequences
