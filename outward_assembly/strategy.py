from typing import Callable, List

from .actions import decrease_k, increase_k, next_priority
from .strategy_helper import (
    CONTIG_COUNT,
    LONGEST_CONTIG,
    READ_PAIR_COUNT,
    OutwardAssemblyMetrics,
)

"""
Strategies: 
User can use constants from strategy_helper.py to define their own strategies, 
and the predefined actions to define the actions to take. 
The user can provide a name for their strategy, or it will default to the function name. 
They will then specify this name in the config file.

Condition variables:
- **CONTIG_COUNT**: The number of contigs produced by Outward Assembly.
- **LONGEST_CONTIG**: The length of the longest contig produced by Outward Assembly.
- **READ_PAIR_COUNT**: The number of reads processed by Outward Assembly. 
  If the user has specified an adapter path or high frequency kmers, 
  then this count will reflect the number of read pairs after each of those.

Actions:
- **next_priority**: Advances the current dataset selection to the next priority.
- **decrease_k**: Decreases the k-mer size used for read filtering for the current dataset.
- **increase_k**: Increases the k-mer size used for read filtering for the current dataset.
"""


def example_strategy(
    outward_assembly_metrics: OutwardAssemblyMetrics,
) -> List[Callable]:
    """
    Example strategy function that demonstrates how to define conditions and actions.
    Returns a list of actions to take based on the metrics.
    """
    actions = []

    if (
        outward_assembly_metrics[READ_PAIR_COUNT] < 100
        and outward_assembly_metrics[LONGEST_CONTIG] < 200
    ):
        actions.append(next_priority())
    if outward_assembly_metrics[CONTIG_COUNT] < 2:
        actions.append(decrease_k(3))
    if outward_assembly_metrics[READ_PAIR_COUNT] > 1000:
        actions.append(increase_k(5))
    return actions


"""
User can define their own strategies below by following the format above.
"""
