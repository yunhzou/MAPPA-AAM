"""Gaussian-format normal-mode input adapter."""
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
