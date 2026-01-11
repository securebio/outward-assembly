from enum import Enum


class SeqOrientation(int, Enum):
    """
    Orientation of a sequence (forward or reverse compliment) relative to some reference
    sequence.
    """

    FORWARD = 1
    REVERSE = -1

    def __mul__(self, other: "SeqOrientation | int") -> "SeqOrientation":
        if not isinstance(other, SeqOrientation):
            other = SeqOrientation(other)
        return SeqOrientation(self.value * other.value)

    def __rmul__(self, other: "SeqOrientation | int") -> "SeqOrientation":
        if not isinstance(other, SeqOrientation):
            other = SeqOrientation(other)
        return SeqOrientation(self.value * other.value)
