import warnings
from collections import deque
from typing import Dict, NamedTuple, Sequence

import networkx as nx
from Bio.Seq import Seq

from outward_assembly.basic_seq_operations import SeqOrientation


class OverlapGraphDiagnostics(NamedTuple):
    graph: nx.Graph
    connected_orientations: Dict[int, SeqOrientation]


def _seqs_overlap_single(
    back_seq: str | Seq,
    front_seq: str | Seq,
    n_0_error: int,
    n_1_error: int,
    max_overlap: int | None = None,
) -> bool:
    """Return whether two sequences overlap by at least n_0_error bases with no errors
    or n_1_error bases with at most 1 error. Only considers overlaps where the end
    of the back sequence overlaps with the start of the forward sequence, like:

    back:    GATTACA
    forward:     ACATTG

    Note that no checks for reverse complements are done, nor does this function
    consider beginning of back_seq overlapping end of front_seq.

    Args:
        back_seq: The sequence that may overlap at its end
        front_seq: The sequence that may overlap at its start
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch
        max_overlap: Maximum overlap length to check (defaults to min length of sequences)

    Returns:
        True if sequences overlap with given criteria, False otherwise
    """
    back_seq = Seq(back_seq) if isinstance(back_seq, str) else back_seq
    front_seq = Seq(front_seq) if isinstance(front_seq, str) else front_seq

    if max_overlap is None:
        max_overlap = min(len(back_seq), len(front_seq))

    for overlap_len in range(n_0_error, max_overlap + 1):
        back_suffix = back_seq[-overlap_len:]
        front_prefix = front_seq[:overlap_len]
        disagreements = sum(b != f for b, f in zip(back_suffix, front_prefix))

        if disagreements == 0 or (overlap_len >= n_1_error and disagreements == 1):
            return True

    return False


def get_seq_overlap_orientation(
    seq1: str | Seq,
    seq2: str | Seq,
    n_0_error: int,
    n_1_error: int,
) -> SeqOrientation | None:
    """
    Checks whether two sequences overlap, and returns their relative orientation if so
    (forward or reverse complement). The sequences are considered to overlap if they have
    an overlapping region of n_0_error bases with no errors or n_1_error bases with at most
    1 error.

    Args:
        seq1: First sequence to check for overlap
        seq2: Second sequence to check for overlap
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch
        relative_orientation: The relative orientations of the sequences to check. If
            Forward, checks whether the sequences overlap in forward orientation relative
            to each other. If Reverse, checks whether the sequences overlap in the reverse
            compliment relative to each other.

    Returns:
        True if sequences overlap with given criteria, False otherwise

    Raises:
        ValueError: If n_1_error < n_0_error
    """
    if n_1_error < n_0_error:
        raise ValueError("n_1_error must be at least n_0_error")

    seq1 = Seq(seq1) if isinstance(seq1, str) else seq1
    seq2 = Seq(seq2) if isinstance(seq2, str) else seq2

    seq2_rc = seq2.reverse_complement()
    seq1_rc = seq1.reverse_complement()

    if _seqs_overlap_single(seq1, seq2, n_0_error, n_1_error) or _seqs_overlap_single(
        seq1_rc, seq2_rc, n_0_error, n_1_error
    ):
        return SeqOrientation.FORWARD
    elif _seqs_overlap_single(seq1, seq2_rc, n_0_error, n_1_error) or _seqs_overlap_single(
        seq1_rc, seq2, n_0_error, n_1_error
    ):
        return SeqOrientation.REVERSE
    else:
        return None


def _overlap_graph(seqs: Sequence[str | Seq], n_0_error: int, n_1_error: int) -> nx.Graph:
    """Create simple graph with seqs at vertices and edges when sequences overlap.

    Two sequences are considered overlapping if they share n_0_error bases with
    no errors or n_1_error bases with at most 1 error.

    Args:
        seqs: Sequence of sequences to build graph from
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        NetworkX Graph where vertex i corresponds to seqs[i], edges correspond to
        overlapping seqs, and each edge is annotated with the attribute "orientation",
        whose value is the relative SeqOrientation (forward or reverse complement) of the
        two overlapping seqs
    """
    seqs = [Seq(s) if isinstance(s, str) else s for s in seqs]

    if len(set(str(seq) for seq in seqs)) != len(seqs):
        warnings.warn("Input sequences to overlap_graph are not unique", UserWarning)

    g = nx.Graph()
    for i in range(len(seqs)):
        g.add_node(i)
    for i, seq1 in enumerate(seqs[:-1]):
        for j, seq2 in enumerate(seqs[i + 1 :], start=i + 1):
            # Annotate overlap orientation as an attribute on the edge
            overlap_orientation = get_seq_overlap_orientation(
                seq1, seq2, n_0_error, n_1_error
            )
            if overlap_orientation is not None:
                # These sequences do overlap
                g.add_edge(i, j, orientation=overlap_orientation)
    return g


def _traverse_subgraph_and_orient(
    g: nx.Graph, reference_seq_idx: int, reference_seq_orientation: SeqOrientation
) -> Dict[int, SeqOrientation]:
    """
    Given a graph of overlapping sequences, and the index of a reference sequence, finds
    all sequences that are connected to the reference, and returns them along with their
    orientations with respect to the reference. The graph is traversed via BFS to find
    all nodes connected to the reference.

    Args:
        g: A networkx graph whose nodes represent sequences. For every pair of nodes A and
           B, there is an edge between A and B iff the sequences of A and B overlap, and
           the edge is annotated with the metadata key "orientation", whose value is the
           relative orientation those sequences (forward or reverse compliment)
        reference_seq_idx: The sequence for which to find all connected sequences
        reference_seq_orientation: The initial orientation of the reference sequence

    Returns:
        A dict whose keys are the indices of sequences that are connected to the reference
        sequence, and whose values are the the orientation of each sequence with respect
        to the reference orientation
    """
    orientations: Dict[int, SeqOrientation] = {reference_seq_idx: reference_seq_orientation}
    # Queue of nodes whose neighbors we want to traverse next in the breadth-first search
    bfs_queue = deque([reference_seq_idx])
    visited = {reference_seq_idx}

    while bfs_queue:
        current = bfs_queue.popleft()
        current_orientation = orientations[current]
        for neighbor in g.neighbors(current):
            if neighbor not in visited:
                # The relative orientation between current and neighbor is stored as an
                # attribute on the edge between them
                neighbor_relative_orientation: SeqOrientation = g[current][neighbor][
                    "orientation"
                ]
                neighbor_canonical_orientation = (
                    current_orientation * neighbor_relative_orientation
                )
                orientations[neighbor] = neighbor_canonical_orientation
                visited.add(neighbor)
                bfs_queue.append(neighbor)

    return orientations


def get_overlap_graph_diagnostics(
    seqs: Sequence[str | Seq],
    seqs_containing_seed: Dict[int, SeqOrientation],
    n_0_error: int,
    n_1_error: int,
) -> OverlapGraphDiagnostics:
    """
    Build the overlap graph and find sequences connected to any seed-containing sequence.

    Args:
        seqs: Sequences to analyze
        seqs_containing_seed: Dict to use to subset the above sequences. The keys of this
            dict are indices of seqs, corresponding to sequences that contain a seed, and
            the values are the orientation of each sequence relative to the seed it contains
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        OverlapGraphDiagnostics containing the full overlap graph and a dict whose keys are
        the indices of any sequence that is connected via overlap to a seed-containing
        sequence, and whose values are the relative orientation of the sequence with
        respect to the seed.
    """
    seqs = [Seq(s) if isinstance(s, str) else s for s in seqs]

    # Construct a graph where nodes are sequences and edges represent overlaps between
    # sequences
    g: nx.Graph = _overlap_graph(seqs, n_0_error, n_1_error)

    all_connected_contigs: Dict[int, SeqOrientation] = {}

    # Search for sequences that are connected to seed-containing sequences
    for seed_contig_idx, seed_contig_orientation in seqs_containing_seed.items():
        if seed_contig_idx in all_connected_contigs:
            # We've already seen this sequence before in a previous connected component, so
            # don't process it again
            continue

        # Get all sequences that are connected to this seed-containing sequence, along with
        # their orientation relative to this seed-containing sequence
        oriented_connected_contigs = _traverse_subgraph_and_orient(
            g,
            reference_seq_idx=seed_contig_idx,
            reference_seq_orientation=seed_contig_orientation,
        )
        all_connected_contigs.update(oriented_connected_contigs)

    return OverlapGraphDiagnostics(
        graph=g,
        connected_orientations=all_connected_contigs,
    )


def get_overlapping_sequence_ids(
    seqs: Sequence[str | Seq],
    seqs_containing_seed: Dict[int, SeqOrientation],
    n_0_error: int,
    n_1_error: int,
) -> Dict[int, SeqOrientation]:
    """
    Get indices of sequences that are connected to any seed-containing sequence, where
    connectedness is determined by the sequences having an overlapping region. Returns
    0-based indices into the provided seqs list, along with the relative orientation of
    each sequence with respect to the seed.

    Args:
        seqs: Sequences to analyze
        seqs_containing_seed: Dict to use to subset the above sequences. The keys of this
            dict are indices of seqs, corresponding to sequences that contain a seed, and
            the values are the orientation of each sequence relative to the seed it contains
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        Dict whose keys are the indices of any sequence that is connected via overlap to a
        seed-containing sequence, and whose values are the relative orientation of the
        sequence with respect to the seed
    """
    return get_overlap_graph_diagnostics(
        seqs,
        seqs_containing_seed,
        n_0_error,
        n_1_error,
    ).connected_orientations
