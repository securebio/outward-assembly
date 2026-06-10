import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from outward_assembly.basic_seq_operations import SeqOrientation
from outward_assembly.pipeline_steps import _contig_ids_by_seed_ahocorasick


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_basic():
    """
    Tests that _contig_ids_by_seed_ahocorasick returns the subset of contigs that
    contain seeds, along with each contig's orientation with respect to the seed.
    """
    seed = Seq("AAAT")
    contigs = [
        SeqRecord(seq=Seq("AAATCGCGCGCG")),  # contains seed in FWD orientation
        SeqRecord(seq=Seq("ATTTGCGCGCGC")),  # contains seed in RC orientation
        SeqRecord(seq=Seq("CCCCCGGGGGG")),  # does not contain seed
        SeqRecord(seq=Seq("AAATATTTCCG")),  # contains seed in FWD and then RC
    ]
    result = _contig_ids_by_seed_ahocorasick(records=contigs, seed_seqs=[seed])
    assert result == {
        0: SeqOrientation.FORWARD,
        1: SeqOrientation.REVERSE,
        3: SeqOrientation.FORWARD,
    }


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_empty_inputs():
    """Test behavior with empty inputs."""
    assert _contig_ids_by_seed_ahocorasick(records=[], seed_seqs=[]) == {}
    assert _contig_ids_by_seed_ahocorasick(records=[], seed_seqs=[Seq("ACGT")]) == {}
    assert (
        _contig_ids_by_seed_ahocorasick(records=[SeqRecord(seq=Seq("ACGT"))], seed_seqs=[])
        == {}
    )


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_multiple_seeds():
    """Test with multiple seed sequences."""
    seeds = [Seq("GATTACA"), Seq("ACGTACGT")]
    contigs = [
        SeqRecord(seq=Seq("AAAGATTACAAAA")),
        SeqRecord(seq=Seq("AAATGTAATCAAA")),  # (TGTAATC == RC of GATTACA)
        SeqRecord(seq=Seq("AAAACGTACGTAA")),
        SeqRecord(seq=Seq("AAAAAAAAAAAAA")),
    ]
    result = _contig_ids_by_seed_ahocorasick(records=contigs, seed_seqs=seeds)
    assert set(result.keys()) == {0, 1, 2}
    assert result[0] == SeqOrientation.FORWARD
    assert result[1] == SeqOrientation.REVERSE
    assert result[2] == SeqOrientation.FORWARD


@pytest.mark.fast
@pytest.mark.unit
@pytest.mark.parametrize(
    ("contig_seq", "expected_orientation"),
    [
        # Two reverse-complement seed hits should orient the contig in reverse even when a
        # forward hit appears first.
        ("AAATCCCATTTCCTATTTGGG", SeqOrientation.REVERSE),
        # Two forward seed hits should orient the contig forward even when a reverse hit
        # appears first.
        ("ATTTCCTAAATCCCAAATGGG", SeqOrientation.FORWARD),
        # Equal forward/reverse counts keep the previous leftmost-match tie behavior.
        ("ATTTCCTAAATGGG", SeqOrientation.REVERSE),
    ],
)
def test_contig_ids_by_seed_uses_majority_orientation(contig_seq, expected_orientation):
    """Orient contigs by majority seed direction, with leftmost-match tie breaking."""
    result = _contig_ids_by_seed_ahocorasick(
        records=[SeqRecord(seq=Seq(contig_seq))],
        seed_seqs=[Seq("AAAT")],
    )

    assert result == {0: expected_orientation}


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_palindrome():
    """Test with palindromic seed (RC equals forward)."""
    seed = Seq("ACGT")  # RC is also ACGT
    contigs = [
        SeqRecord(seq=Seq("AAACGTAAA")),  # Contains palindrome
        SeqRecord(seq=Seq("AAAAAAA")),  # No match
    ]
    result = _contig_ids_by_seed_ahocorasick(records=contigs, seed_seqs=[seed])
    assert 0 in result
    assert 1 not in result


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_variable_length():
    """Test with variable length seeds."""
    seeds = [
        Seq("GATTACA"),  # 7bp
        Seq("ACGTACGTACGT"),  # 12bp
        Seq("AT"),  # 2bp
    ]
    contigs = [
        SeqRecord(seq=Seq("AAAGATTACAAAA")),  # Has seed 0 (7bp)
        SeqRecord(seq=Seq("ACGTACGTACGTAA")),  # Has seed 1 (12bp)
        SeqRecord(seq=Seq("CCCCCCCCCCCCCC")),  # No match (no AT either)
        SeqRecord(seq=Seq("GGGGATGGGG")),  # Has seed 2 (AT)
    ]
    result = _contig_ids_by_seed_ahocorasick(records=contigs, seed_seqs=seeds)
    assert set(result.keys()) == {0, 1, 3}
    assert result[0] == SeqOrientation.FORWARD
    assert result[1] == SeqOrientation.FORWARD
    assert result[3] == SeqOrientation.FORWARD


def _naive_contig_ids_by_seed(records, seed_seqs):
    """
    Naive O(seeds × contigs) implementation for comparison.
    Returns the majority seed orientation for each contig, breaking ties by leftmost match.
    """
    filtered_records = {}
    for i, rec in enumerate(records):
        contig_sequence = str(rec.seq)
        orientation_counts = {
            SeqOrientation.FORWARD: 0,
            SeqOrientation.REVERSE: 0,
        }
        first_match_by_orientation = {}

        for seed in seed_seqs:
            seed_str = str(seed)
            seed_rc = str(seed.reverse_complement())
            for orientation, seed_to_find in (
                (SeqOrientation.FORWARD, seed_str),
                (SeqOrientation.REVERSE, seed_rc),
            ):
                if seed_str == seed_rc and orientation == SeqOrientation.REVERSE:
                    continue

                start = contig_sequence.find(seed_to_find)
                while start != -1:
                    orientation_counts[orientation] += 1
                    first_match_by_orientation.setdefault(orientation, start)
                    start = contig_sequence.find(seed_to_find, start + 1)

        if not first_match_by_orientation:
            continue

        forward_count = orientation_counts[SeqOrientation.FORWARD]
        reverse_count = orientation_counts[SeqOrientation.REVERSE]
        if forward_count > reverse_count:
            filtered_records[i] = SeqOrientation.FORWARD
        elif reverse_count > forward_count:
            filtered_records[i] = SeqOrientation.REVERSE
        else:
            filtered_records[i] = min(
                first_match_by_orientation,
                key=lambda orientation: first_match_by_orientation[orientation],
            )

    return filtered_records


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_random_data_verification():
    """
    Generate random data and verify that Aho-Corasick matches the naive implementation.
    """
    import random

    random.seed(42)

    bases = "ACGT"
    num_seeds = 50
    num_contigs = 200

    # generate some seeds
    seeds = []
    for _ in range(num_seeds):
        length = random.randint(15, 40)
        seq = "".join(random.choices(bases, k=length))
        seeds.append(Seq(seq))

    # generate some contigs
    contigs = []
    for i in range(num_contigs):
        length = random.randint(200, 600)
        seq = "".join(random.choices(bases, k=length))
        contigs.append(SeqRecord(Seq(seq), id=f"contig_{i}"))

    # guarantee that some contigs contain seeds
    for i in range(0, min(30, num_contigs)):
        seed = str(random.choice(seeds))
        contig_seq = str(contigs[i].seq)
        if len(contig_seq) >= len(seed):
            pos = random.randint(0, len(contig_seq) - len(seed))
            new_seq = contig_seq[:pos] + seed + contig_seq[pos + len(seed) :]
            contigs[i].seq = Seq(new_seq)

    naive_result = _naive_contig_ids_by_seed(contigs, seeds)
    ac_result = _contig_ids_by_seed_ahocorasick(contigs, seeds)

    # check indices
    assert set(naive_result.keys()) == set(ac_result.keys()), (
        f"Different contig indices: naive={set(naive_result.keys())}, "
        f"ac={set(ac_result.keys())}"
    )

    # check orientations
    for key in naive_result.keys():
        assert naive_result[key] == ac_result[key], (
            f"Orientation mismatch for contig {key}: "
            f"naive={naive_result[key]}, ac={ac_result[key]}"
        )
