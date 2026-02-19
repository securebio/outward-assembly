import pytest
from Bio import SeqIO

from outward_assembly.io_helpers import process_s3_paths
from outward_assembly.kmer_freq_filter import high_freq_kmers_split_files


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.requires_tools
@pytest.mark.requires_aws
def test_high_freq_kmer(temp_workdir):
    """Integration test of high_freq_kmers_split_files.

    Reads a tiny file in the testing s3 bucket. This file has 8 FASTQ records (for this
    test the paired end nature isn't relevant), and each sequence ends in AAACCCGGG. So
    with a min_kmer_freq of 8 we should see exactly one high freq kmer: AAACCCGGG."""
    s3_records = process_s3_paths(
        ["s3://nao-testing/outward-assembly-test-data/8_seqs_ending_aaacccggg.fastq.zst"]
    )
    hfk_path = high_freq_kmers_split_files(s3_records, temp_workdir, k=9, min_kmer_freq=8)
    records = list(SeqIO.parse(hfk_path, "fasta"))
    assert len(records) == 1
    assert str(records[0].seq) == "AAACCCGGG"
    assert records[0].id == "kmer_max_8_obs"
