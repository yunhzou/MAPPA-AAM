"""Mode parsing and coordinate reindexing; scoring lives only in mappa.scoring."""
from __future__ import annotations
import re
from pathlib import Path

import numpy as np


def parse_g98_modes(path):
    """Parse vibration modes from xtb's g98.out (Gaussian 98 format).

    Returns:
      freqs:   (n_modes,) float array, cm^-1, negative = imaginary
      modes:   (n_modes, n_atoms, 3) float array of Cartesian displacements
    """
    text = Path(path).read_text()
    lines = text.splitlines()
    freqs_all = []
    blocks = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith('Frequencies --'):
            tokens = line.split('--', 1)[1].split()
            block_freqs = [float(t) for t in tokens]
            j = i + 1
            while j < len(lines) and not lines[j].lstrip().startswith('Atom'):
                j += 1
            j += 1  # skip header
            atoms_disp = []  # rows of [d1x,d1y,d1z,d2x,d2y,d2z,d3x,d3y,d3z]
            while j < len(lines):
                ln = lines[j].strip()
                if not ln or ln.startswith('Frequencies'):
                    break
                parts = ln.split()
                if len(parts) < 2 + 3 * len(block_freqs):
                    break
                vals = [float(x) for x in parts[2:2 + 3 * len(block_freqs)]]
                atoms_disp.append(vals)
                j += 1
            arr = np.asarray(atoms_disp, dtype=float)
            n_atoms = arr.shape[0]
            for k in range(len(block_freqs)):
                disp = arr[:, 3 * k:3 * (k + 1)]  # (n_atoms, 3)
                freqs_all.append(block_freqs[k])
                blocks.append(disp)
            i = j
        else:
            i += 1
    if not blocks:
        return np.zeros(0), np.zeros((0, 0, 3))
    freqs = np.asarray(freqs_all, dtype=float)
    modes = np.stack(blocks, axis=0)  # (n_modes, n_atoms, 3)
    return freqs, modes


def bond_reaction_vector(xyz_TS_in_R, broken_bonds_R, formed_bonds_R):
    """Per-atom reaction direction at TS coordinates.

    For each broken bond (i, j): atoms should be moving APART --
    V[i] -= unit, V[j] += unit  where unit = (xyz[j]-xyz[i])/||.||.
    For each formed bond (i, j): atoms should be moving TOGETHER --
    V[i] += unit, V[j] -= unit.

    Sign convention is chosen so that a true concerted reaction mode
    accumulates contributions coherently. Returns (n_atoms, 3) ndarray,
    NOT normalized."""
    xyz = np.asarray(xyz_TS_in_R, dtype=float)
    V = np.zeros_like(xyz)
    for i, j in broken_bonds_R:
        v = xyz[j] - xyz[i]
        n = float(np.linalg.norm(v))
        if n < 1e-9: continue
        u = v / n
        V[i] -= u
        V[j] += u
    for i, j in formed_bonds_R:
        v = xyz[j] - xyz[i]
        n = float(np.linalg.norm(v))
        if n < 1e-9: continue
        u = v / n
        V[i] += u
        V[j] -= u
    return V


def core_atoms_in_R_frame(mapping_R_to_P, broken, formed):
    """Atoms touching broken (R indices) or formed (P, mapped back via inv).
    Returns sorted list of R-frame indices."""
    inv = {v: k for k, v in mapping_R_to_P.items()}
    core = set()
    for (i, j, _, _) in broken:
        core.add(int(i)); core.add(int(j))
    for (i, j, _, _) in formed:
        if i in inv: core.add(int(inv[i]))
        if j in inv: core.add(int(inv[j]))
    return sorted(core)


def reindex_modes_to_R(modes_TS, mapping_R_to_TS, n_R):
    """Reindex mode displacements from TS-frame to R-frame indexing,
    so modes_R[k][i_R] == modes_TS[k][mapping_R_to_TS[i_R]]. Returns
    (n_modes, n_R, 3)."""
    n_modes = modes_TS.shape[0]
    out = np.zeros((n_modes, n_R, 3))
    for r, t in mapping_R_to_TS.items():
        out[:, r] = modes_TS[:, t]
    return out
