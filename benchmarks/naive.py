"""
Naive O(seeds × contigs) implementation for seed-to-contig matching.

This is the original implementation, kept here for benchmarking comparison.
"""

from typing import Dict, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from outward_assembly.basic_seq_operations import SeqOrientation


def contig_ids_by_seed(
    records: List[SeqRecord], seed_seqs: List[Seq]
) -> Dict[int, SeqOrientation]:
    """
    Given a list of contigs and a list of seed sequences, returns the indices of each contig
    that contains at least one seed, along with its orientation with respect to the seed
    (forward or reverse complement). (If a contig has multiple seeds, we return its
    orientation with respect to the first seed in the contig.)

    Args:
        records: List of SeqRecord objects representing contigs
        seed_seqs: List of seed sequences to search for
    Returns:
        Dict whose keys correspond to the indices of records that contain seed sequences,
        and whose values correspond to the contig orientation with respect to the seed
    """
    filtered_records = {}
    for i, rec in enumerate(records):
        contig_sequence = str(rec.seq)
        for seed in seed_seqs:
            seed_str = str(seed)
            seed_rc = str(seed.reverse_complement())
            if seed_str in contig_sequence:
                filtered_records[i] = SeqOrientation.FORWARD
                break
            elif seed_rc in contig_sequence:
                filtered_records[i] = SeqOrientation.REVERSE
                break

    return filtered_records
