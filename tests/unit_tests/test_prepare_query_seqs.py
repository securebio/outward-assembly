import pytest

from outward_assembly.pipeline_steps import _prepare_query_seqs


@pytest.mark.fast
@pytest.mark.unit
def test_prepare_query_seqs_combines_seed_and_warm_start(tmp_path):
    seed = tmp_path / "current_contigs.fasta"
    seed.write_text(">seed\nACGT\n")
    warm = tmp_path / "warm.fasta"
    warm.write_text(">warm\nTGCA\n")

    result = _prepare_query_seqs(tmp_path, warm, adapters_path=None, read_subset_k=27)
    content = result.read_text()

    assert ">seed" in content
    assert ">warm" in content
    assert "ACGT" in content
    assert "TGCA" in content


@pytest.mark.fast
@pytest.mark.unit
def test_prepare_query_seqs_no_warm_start(tmp_path):
    seed = tmp_path / "current_contigs.fasta"
    seed.write_text(">seed\nACGT\n")

    result = _prepare_query_seqs(tmp_path, None, adapters_path=None, read_subset_k=27)
    assert result == seed
