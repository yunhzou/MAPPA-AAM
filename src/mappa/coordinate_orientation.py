"""Coordinate orientation measurements; no mapping search or witness selection."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import numpy as np

ORIENTATION_DEGENERACY_TOL = 0.1


class OrientationError(ValueError):
    """Invalid coordinate-orientation input."""


@dataclass(frozen=True)
class _OrientationMeasure:
    normalized: float
    determinant: float
    determinant_error_bound: float
    sign: int
    zero_length: bool


def _orientation_measure(coords, origin, other_points, *,
                         degeneracy_tol=ORIENTATION_DEGENERACY_TOL):
    xyz = np.asarray(coords, dtype=np.longdouble)
    vectors = np.stack(
        [xyz[int(point)] - xyz[int(origin)] for point in other_points],
        axis=0,
    )
    squared_lengths = np.sum(vectors * vectors, axis=1)
    zero_length = bool(np.any(squared_lengths == 0))
    denominator = np.sqrt(np.prod(squared_lengths))
    a, b, c = vectors
    positive = (
        a[0] * b[1] * c[2],
        a[1] * b[2] * c[0],
        a[2] * b[0] * c[1],
    )
    negative = (
        a[2] * b[1] * c[0],
        a[1] * b[0] * c[2],
        a[0] * b[2] * c[1],
    )
    determinant = sum(positive) - sum(negative)
    permanent = sum(abs(term) for term in (*positive, *negative))
    eps = np.longdouble(np.finfo(np.longdouble).eps)
    gamma = 16 * eps / (1 - 16 * eps)
    error_bound = gamma * max(permanent, denominator)
    normalized = 0 if zero_length else determinant / denominator
    if (zero_length or abs(determinant) <= error_bound
            or abs(normalized) <= float(degeneracy_tol)):
        sign = 0
    else:
        sign = 1 if determinant > 0 else -1
    return _OrientationMeasure(
        normalized=float(normalized),
        determinant=float(determinant),
        determinant_error_bound=float(error_bound),
        sign=sign,
        zero_length=zero_length,
    )


def _simplex_measure(coords, center, neighbors, degeneracy_tol):
    """Affine orientation of one local substituent simplex."""
    neighbors = tuple(int(atom) for atom in neighbors)
    if len(neighbors) == 3:
        origin, points = int(center), neighbors
        return _orientation_measure(
            coords, origin, points, degeneracy_tol=degeneracy_tol)
    if len(neighbors) != 4:
        raise OrientationError(
            "an oriented substituent simplex needs three or four neighbors")

    # A four-ligand simplex is an affine tetrahedron.  Normalizing its volume
    # by the three edges incident to ``neighbors[0]`` makes the degeneracy
    # decision depend on which ligand happens to be listed first.  Use one
    # canonical determinant and the permutation-invariant geometric mean of
    # all six tetrahedron edges instead; recover ordered orientation solely
    # from permutation parity.
    canonical = tuple(sorted(neighbors))
    base = _orientation_measure(
        coords, canonical[0], canonical[1:], degeneracy_tol=0.0)
    xyz = np.asarray(coords, dtype=np.longdouble)
    squared_distances = []
    for left, right in combinations(canonical, 2):
        delta = xyz[left] - xyz[right]
        squared_distances.append(np.sum(delta * delta))
    zero_length = bool(any(value == 0 for value in squared_distances))
    denominator = (
        np.longdouble(0.0) if zero_length else
        np.prod(np.asarray(squared_distances, dtype=np.longdouble))
        ** np.longdouble(0.25))
    positions = {atom: index for index, atom in enumerate(canonical)}
    order = tuple(positions[atom] for atom in neighbors)
    inversions = sum(
        order[left] > order[right]
        for left in range(4) for right in range(left + 1, 4))
    parity = -1 if inversions % 2 else 1
    determinant = np.longdouble(parity) * np.longdouble(base.determinant)
    normalized = (
        np.longdouble(0.0) if zero_length else determinant / denominator)
    if (zero_length
            or abs(determinant) <= base.determinant_error_bound
            or abs(normalized) <= float(degeneracy_tol)):
        sign = 0
    else:
        sign = 1 if determinant > 0 else -1
    return _OrientationMeasure(
        normalized=float(normalized),
        determinant=float(determinant),
        determinant_error_bound=float(base.determinant_error_bound),
        sign=sign,
        zero_length=zero_length,
    )
