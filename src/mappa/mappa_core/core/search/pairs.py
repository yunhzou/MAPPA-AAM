"""Immutable atom-pair sequences sharing unchanged blocks between search states."""
from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import chain


@dataclass(frozen=True, eq=False, slots=True)
class SharedPairs(Sequence):
    blocks: tuple
    size: int
    _hash: object = field(default=None, init=False, repr=False)

    def __iter__(self):
        return chain.from_iterable(self.blocks)

    def __len__(self):
        return self.size

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(self)[index]
        if index < 0:
            index += self.size
        if index < 0 or index >= self.size:
            raise IndexError(index)
        for block in self.blocks:
            if index < len(block):
                return block[index]
            index -= len(block)
        raise IndexError(index)

    def __eq__(self, other):
        if isinstance(other, SharedPairs):
            if self.blocks == other.blocks:
                return True
            return len(self) == len(other) and all(a == b for a, b in zip(self, other))
        if isinstance(other, tuple):
            return tuple(self) == other
        return NotImplemented

    def __hash__(self):
        if self._hash is None:
            object.__setattr__(self, '_hash', hash(tuple(self)))
        return self._hash

    def __deepcopy__(self, memo):
        return self

    def __reduce__(self):
        # Persist shared blocks directly; the cached hash is derived local data.
        return type(self), (self.blocks, self.size)
