"""Shared element classification and conditional bond thresholds."""
METAL_ELEMENTS = frozenset({
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs",
    "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv",
})


def is_metal_element(element):
    """Return True for elements treated with the metal WBO-event cutoff."""
    if element is None:
        return False
    s = str(element).strip()
    if not s:
        return False
    return s[0].upper() + s[1:].lower() in METAL_ELEMENTS


def bond_event_threshold(elements, i, j, *,
                         default_threshold=0.5,
                         metal_threshold=None):
    """Pair-specific delta-WBO threshold for mechanism events.

    The default threshold is used for ordinary covalent/organic pairs.  If a
    metal threshold is supplied and either endpoint is a metal, the lower
    metal-aware cutoff is used.
    """
    if elements is None or metal_threshold is None:
        return float(default_threshold)
    if is_metal_element(elements[i]) or is_metal_element(elements[j]):
        return float(metal_threshold)
    return float(default_threshold)
