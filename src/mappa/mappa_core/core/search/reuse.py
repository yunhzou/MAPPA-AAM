"""Dependency-checked reuse across cut searches."""

from __future__ import annotations
from dataclasses import dataclass
from mappa.mappa_core.core.matching import native


@dataclass(frozen=True)
class CutReplayView:
    source: object
    target: object
    source_view: object
    target_view: object
    engine: object
    cuts: tuple


class _CutSession:
    def __init__(self, source, target, target_orbits, bond_tolerance=None):
        self.source = source
        self.target = target
        self.source_view = native.source_graph(source,
            None if bond_tolerance is None else bond_tolerance.resolve(source))
        self.target_view = native.target_graph(target, target_orbits)
        if self.source_view is None or self.target_view is None:
            raise ValueError("cut replay requires native WBO graphs and exact target orbits")

    def for_cut(self, cut):
        index = self.source_view.index
        source = self.source.copy()
        source.remove_edges_from(cut)
        return CutReplayView(source, self.target, self.source_view, self.target_view, self.engine,
            tuple((index[a], index[b]) for a, b in cut))

    def stats(self):
        return dict(self.engine.stats())


class CutReplay(_CutSession):
    def __init__(self, source, target, target_orbits, *, cache_bytes=64*1024*1024,
                 checkpoints=True, reference_only=False, bond_tolerance=None):
        super().__init__(source, target, target_orbits, bond_tolerance)
        self.engine = native._engine.GrowthReplay(
            self.source_view.graph, self.target_view.graph, cache_bytes, checkpoints, reference_only)


class FragmentRepair(_CutSession):
    """Reuse a later fragment only when its full conditional result is unchanged.

    Target occupancy and observed source/boundary inputs are dependencies.
    Changed earlier mapping histories do not themselves invalidate an otherwise
    independent fragment. Search ordering and admission remain the caller's.
    """
    def __init__(self, source, target, target_orbits, *, cache_bytes=64*1024*1024,
                 extension_cache_bytes=0, bond_tolerance=None):
        super().__init__(source, target, target_orbits, bond_tolerance)
        self.engine = native._engine.FragmentRepair(
            self.source_view.graph, self.target_view.graph, cache_bytes, extension_cache_bytes)
