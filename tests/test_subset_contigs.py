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
        SeqRecord(
            seq=Seq("AAATGTAATCAAA")
        ),  # (TGTAATC == RC of GATTACA)
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


@pytest.mark.fast
@pytest.mark.unit
def test_contig_ids_by_seed_random_data_verification():
    """
    Generate random data and verify that Aho-Corasick matches the naive implementation.
    """
    import random
    import sys
    from pathlib import Path

    # grab the naive implementation
    benchmarks_dir = Path(__file__).parent.parent / "benchmarks"
    sys.path.insert(0, str(benchmarks_dir))

    from naive import contig_ids_by_seed as naive_contig_ids_by_seed

    sys.path.pop(0)
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


    naive_result = naive_contig_ids_by_seed(contigs, seeds)
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
