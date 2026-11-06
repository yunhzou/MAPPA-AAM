"""Unique completed fragments shared by extension queries across cuts."""

from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
import hashlib
from mappa.mappa_core.core.matching.growth import native


@dataclass(frozen=True)
class FragmentPoolView:
    source: object
    target: object
    source_view: object
    target_view: object
    engine: object
    cuts: tuple


class FragmentPool:
    """Use known fragment witnesses to certify temporary growth fragments."""
    def __init__(self, source, target, *, cache_bytes=64*1024*1024,
                 bond_tolerance=None):
        self.source = source
        self.target = target
        self.source_view = native.source_graph(source,
            None if bond_tolerance is None else bond_tolerance.resolve(source))
        self.target_view = native.target_graph(target)
        if self.source_view is None or self.target_view is None:
            raise ValueError('fragment pool requires native WBO graphs')
        self.engine = native._engine.FragmentPool(
            self.source_view.graph, self.target_view.graph, cache_bytes)

    def for_cut(self, cut):
        index = self.source_view.index
        source = self.source.copy()
        source.remove_edges_from(cut)
        return FragmentPoolView(source, self.target, self.source_view, self.target_view, self.engine,
            tuple((index[a], index[b]) for a, b in cut))

    def stats(self):
        return dict(self.engine.stats())

    @cached_property
    def snapshot_identity(self):
        return hashlib.sha256(self.engine.snapshot_context()).hexdigest()

    def dumps(self):
        return self.engine.dump(self.snapshot_identity)

    def loads(self, data):
        return self.engine.load(data, self.snapshot_identity)
