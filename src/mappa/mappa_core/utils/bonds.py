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


# Main-group centres whose bonds are not limited to two-centre two-electron
# pairs. In a three-centre bond the Wiberg order of each link is about one
# half, so one link forming or breaking changes the order by at most ~0.5 and
# falls under the ordinary cutoff.
#   three-centre four-electron (hypervalent and secondary/sigma-hole bonds):
#     p-block elements from the third period on, which can hold more than
#     four electron pairs around the centre;
#   three-centre two-electron (electron-deficient): boron. Its heavier
#     congeners are already listed as metals.
# H, C, N, O and F keep the ordinary cutoff: their stable endpoints are
# described by two-centre bonds with near-integer orders.
MULTICENTER_ELEMENTS = frozenset({
    "B",
    "Si", "P", "S", "Cl",
    "Ge", "As", "Se", "Br", "Kr",
    "Sb", "Te", "I", "Xe",
    "At", "Rn",
})
FRACTIONAL_BOND_ELEMENTS = METAL_ELEMENTS | MULTICENTER_ELEMENTS


def _symbol(element):
    if element is None:
        return ""
    s = str(element).strip()
    return s[0].upper() + s[1:].lower() if s else ""


def is_metal_element(element):
    """Return True for metals."""
    return _symbol(element) in METAL_ELEMENTS


def is_fractional_bond_element(element):
    """Return True for elements treated with the lower WBO-event cutoff.

    These are metals and the main-group elements that form three-centre
    bonds, where a single link carries a bond order near one half.
    """
    return _symbol(element) in FRACTIONAL_BOND_ELEMENTS


def bond_event_threshold(elements, i, j, *,
                         default_threshold=0.5,
                         fractional_bond_threshold=None):
    """Pair-specific delta-WBO threshold for mechanism events.

    The default threshold is used for ordinary two-centre covalent pairs. If a
    fractional-bond threshold is supplied and either endpoint is a metal or a main-group
    element that forms three-centre bonds, the lower cutoff is used.
    """
    if elements is None or fractional_bond_threshold is None:
        return float(default_threshold)
    if is_fractional_bond_element(elements[i]) or is_fractional_bond_element(elements[j]):
        return float(fractional_bond_threshold)
    return float(default_threshold)
