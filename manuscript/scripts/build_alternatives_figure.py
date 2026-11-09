"""Draw saved EnzymeMap alternatives and check their hydrogen-completed counts."""
from pathlib import Path
from itertools import permutations
import json,os
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[1]/'build/mpl-cache'))
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Circle
from rdkit import Chem
from build_molecule_figure import molecule,bond_event_marker,INK,MUTED,TEAL
R='[C:1]([O:2][H:6])([H:3])([H:4])[H:5]'
A='[C:1](=[O:2])([H:3])[H:4]'
B='[C:1](=[O:2])([H:3])[H:6]'

def check_example(man):
 d=json.loads((man/'evidence/multicandidate_example.json').read_text())
 assert d['complete_saved_families'] and [c['total'] for c in d['candidates']]==[4,6]
 assert d['queries']=={'same_class_H3_H4_swap':'allowed','other_class_preserve_events':'forbidden','other_class_any_saved_events':'allowed'}
 wr=np.array(d['wbo_R']);wp=np.array(d['wbo_P'])
 expected=[[(0,1,1),(0,4,-1),(1,5,-1),(4,5,1)],[(0,1,1),(0,3,-1),(0,4,-1),(0,5,1),(1,5,-1),(3,4,1)]]
 from rdkit import Chem
 def bonds(s):
  p=Chem.SmilesParserParams();p.removeHs=False;m=Chem.MolFromSmiles(s,p);assert m
  return {tuple(sorted((b.GetBeginAtom().GetAtomMapNum()-1,b.GetEndAtom().GetAtomMapNum()-1))):b.GetBondTypeAsDouble() for b in m.GetBonds()}
 r=bonds(R)
 assert len(r)==5 and all(wr[i,j]==v for (i,j),v in r.items())
 for m,expected_events,smiles in zip(d['figure_witnesses_zero_based'],expected,[A+'.[H:5][H:6]',B+'.[H:4][H:5]']):
  m={int(k):v for k,v in m.items()}
  delta=wp[np.ix_([m[i] for i in range(6)],[m[i] for i in range(6)])]-wr
  events=[(i,j,int(np.sign(delta[i,j]))) for i in range(6) for j in range(i+1,6) if abs(delta[i,j])>.5]
  assert events==expected_events
  b=bonds(smiles)
  assert len(b)==4
  assert all(b.get((i,j),0)==wp[m[i],m[j]] for i in range(6) for j in range(i+1,6))
 counts=json.loads((man/'evidence/output_multiplicity.json').read_text())
 assert all(m['cases']==1851 for m in counts['methods'])
 assert all(m['multiple_record_cases']==0 for m in counts['methods'] if not m['method'].startswith('slap'))
 assert all(m['multiple_record_cases']>0 for m in counts['methods'] if m['method'].startswith('slap'))
 return d


def check_case(man):
    """Recount the illustrated saved mappings directly from their endpoint matrices."""
    from collections import Counter

    data = json.loads((man/'evidence/enzymemap_case13175.json').read_text())
    assert data['dataset'] == 'EnzymeMap' and data['case'] == 13175
    source, target = (data['input'][side] for side in ('reactant', 'product'))
    assert Counter(source['elements']) == Counter(target['elements']) == dict(C=4, H=6, O=5)
    weights = [np.asarray(endpoint['wbo']) for endpoint in (source, target)]
    heavy = [i for i, element in enumerate(source['elements']) if element != 'H']
    assert heavy == list(range(9))
    reference = dict(data['reference_mapping'])
    assert set(reference) == set(heavy) and set(reference.values()) == set(heavy)
    assert all(dict(data['patterns'][0]['mapping'])[a] == b for a, b in reference.items())
    params = Chem.SmilesParserParams()
    params.removeHs = False
    for endpoint, expected in zip((source, target), weights):
        mol = Chem.MolFromSmiles(endpoint['smiles'], params)
        actual = np.zeros_like(expected)
        for atom in mol.GetAtoms():
            assert atom.GetSymbol() == endpoint['elements'][atom.GetAtomMapNum()-1]
        for bond in mol.GetBonds():
            a, b = (atom.GetAtomMapNum()-1 for atom in (bond.GetBeginAtom(), bond.GetEndAtom()))
            actual[a, b] = actual[b, a] = bond.GetBondTypeAsDouble()
        assert np.array_equal(actual, expected)
    for pattern in data['patterns']:
        mapping = dict(pattern['mapping'])
        assert sorted(mapping) == sorted(mapping.values()) == list(range(15))
        assert all(source['elements'][a] == target['elements'][b] for a, b in mapping.items())
        assert sorted(a for fragment in pattern['source_fragments'] for a in fragment) == list(range(15))
        order = [mapping[a] for a in range(15)]
        delta = weights[1][np.ix_(order, order)] - weights[0]
        pairs = [(a, b) for a in range(15) for b in range(a+1, 15) if abs(delta[a, b]) > .5]
        heavy_pairs = [(a, b) for a, b in pairs if a in heavy and b in heavy]
        assert len(pairs) == pattern['total']
        assert len(heavy_pairs) == pattern['heavy_total']
        assert sorted(tuple(event['r']) for event in pattern['heavy_events']) == heavy_pairs
        for event in pattern['heavy_events']:
            a, b = event['r']
            before, after = weights[0][a, b], weights[1][mapping[a], mapping[b]]
            assert event['p'] == [mapping[a], mapping[b]]
            assert event['wbo'] == [before, after]
            expected = ('broken' if after == 0 else 'weakened') if before > after else ('formed' if before == 0 else 'strengthened')
            assert event['kind'] == expected
    assert [p['heavy_total'] for p in data['patterns']] == [4, 2]
    assert [p['total'] for p in data['patterns']] == [8, 6]
    assert data['provenance']['minimum_events'] == data['patterns'][1]['total']
    completion = data['hydrogen_completion_check']
    assert completion['scope'] == 'all hydrogen assignments with each heavy-atom mapping fixed'
    hydrogen = [i for i, element in enumerate(source['elements']) if element == 'H']
    target_hydrogen = [i for i, element in enumerate(target['elements']) if element == 'H']
    pairs = np.triu_indices(len(source['elements']), 1)
    for pattern, row in zip(data['patterns'], completion['results'], strict=True):
        fixed = {a: b for a, b in pattern['mapping'] if a in heavy}
        totals = []
        for assigned in permutations(target_hydrogen):
            mapping = fixed | dict(zip(hydrogen, assigned))
            order = [mapping[a] for a in range(len(source['elements']))]
            delta = weights[1][np.ix_(order, order)] - weights[0]
            totals.append(int(np.count_nonzero(np.abs(delta[pairs]) > .5)))
        assert len(totals) == completion['assignments_per_mapping'] == 720
        assert row['name'] == pattern['name']
        assert min(totals) == row['minimum_total'] == 6
        assert row['minimum_hydrogen_events'] == 6-pattern['heavy_total']
        mapping = dict(row['mapping'])
        assert len(mapping) == len(set(mapping.values())) == 15
        assert all(mapping[a] == b for a, b in fixed.items())
        assert all(source['elements'][a] == target['elements'][b] for a, b in mapping.items())
        order = [mapping[a] for a in range(15)]
        delta = weights[1][np.ix_(order, order)] - weights[0]
        assert np.count_nonzero(np.abs(delta[pairs]) > .5) == 6
    return data


def build(man):
    data = check_case(man)
    palette = [(.70, .87, .81), (.81, .68, .92), (.98, .81, .49)]
    params = Chem.SmilesParserParams()
    params.removeHs = False
    endpoints = [Chem.RemoveHs(Chem.MolFromSmiles(data['input'][side]['smiles'], params))
                 for side in ('reactant', 'product')]
    # Fixed atom-origin colors expose the assignment difference across both rows.
    origin_groups = [[1, 2], [3, 4, 5], [6, 7, 8, 9]]
    owners = {atom: color for atoms, color in zip(origin_groups, palette) for atom in atoms}
    carbon_labels = {i+1: str(i+1) for i, element in enumerate(data['input']['reactant']['elements'])
                     if element == 'C'}
    source_positions = {1: (0, 0), 2: (1.5, 0),
                        3: (0, 1.5), 4: (0, 0), 5: (-1.3, -.75),
                        6: (1.3, -.75), 7: (1.3, -2.25), 8: (2.6, 0), 9: (3.9, -.75)}
    target_positions = {1: (0, -1.5), 2: (0, 0), 3: (-1.3, .75),
                        4: (-2.6, 0), 5: (1.3, .75), 6: (2.6, 0),
                        7: (-1.5, 0), 8: (0, 0), 9: (1.5, 0)}

    fig = plt.figure(figsize=(11.8, 6.9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1416), ylim=(828, 0))
    ax.axis('off')

    def text(x, y, label, size=11, color=INK, weight='normal', align='left'):
        ax.text(x, y, label, fontsize=size, color=color, weight=weight, ha=align,
                va='center', linespacing=1.3, zorder=9)

    def box(x, y, width, height, color='#F5F7FB', edge=None):
        ax.add_patch(FancyBboxPatch((x, y), width, height,
            boxstyle='round,pad=0,rounding_size=12', facecolor=color,
            edgecolor=edge or color, lw=.8, zorder=0))

    box(12, 12, 1392, 804)
    text(32, 42, 'Same endpoints, different carbon origins', 16, weight='bold')
    text(32, 74, 'EnzymeMap 13175  ·  complete composition C₄H₆O₅  ·  all 9 heavy atoms shown', 10.5, MUTED)
    for x, color, label in zip((44, 308, 564), palette,
                               ('Formaldehyde', 'Carboxyl group', 'Remaining carbon chain')):
        ax.add_patch(Circle((x, 111), 9, facecolor=color, edgecolor='none'))
        text(x+18, 111, label, 10, MUTED)
    for x, sign, label in ((914, -1, 'Break / weaken'), (1153, 1, 'Form / strengthen')):
        ax.plot([x, x+37], [111, 111], color='black', lw=1.1)
        bond_event_marker(ax, (x, 111), (x+37, 111), sign, size=8)
        text(x+48, 111, label, 9.5)

    titles = ['A   Recovered reference mapping', 'B   MAPPA alternative mapping']
    descriptions = ['Carbon 1 joins the three-carbon molecule; carbon 4 becomes CO₂.',
                    'Carbon 1 becomes CO₂; carbon 4 stays in the three-carbon molecule.']
    event_summaries = ['1 broken + 1 formed + 1 weakened + 1 strengthened',
                       '1 broken + 1 formed']
    for number, pattern in enumerate(data['patterns']):
        y = 144 + number*295
        box(28, y, 1359, 281, 'white', '#C7DFD5' if number == 0 else '#DCE3E9')
        text(47, y+25, titles[number], 12, TEAL if number == 0 else INK, 'bold')
        text(1358, y+25, f"{pattern['heavy_total']} heavy-atom changes", 11, TEAL, 'bold', 'right')
        mapping = dict(pattern['mapping'])
        inverse = {b: a for a, b in mapping.items()}
        layouts = [[(40, 160), (254, 373)], [(755, 363), (1173, 180)]]
        for side in (0, 1):
            components = list(Chem.GetMolFrags(endpoints[side], asMols=True))
            for component, (x, width) in zip(components, layouts[side]):
                positions = {}
                for atom in component.GetAtoms():
                    original = atom.GetAtomMapNum()
                    base = original-1 if side == 0 else inverse[original-1]
                    atom.SetAtomMapNum(base+1)
                    positions[base+1] = (source_positions if side == 0 else target_positions)[original]
                sign = -1 if side == 0 else 1
                kinds = ('broken', 'weakened') if side == 0 else ('formed', 'strengthened')
                changes = {tuple(a+1 for a in event['r']): sign for event in pattern['heavy_events']
                           if event['kind'] in kinds and all(a+1 in positions for a in event['r'])}
                molecule(ax, Chem.MolToSmiles(component), x, y+55, width, 173,
                         owners=owners, notes=carbon_labels, font=24, align_acetyl=False,
                         bond_changes=changes, coordinates=positions,
                         atom_labels={1: 'H<sub>2</sub>C'} if side == 0 else None)
            text(225 if side == 0 else 1142, y+141, '+', 18, MUTED, align='center')
        ax.add_patch(FancyArrowPatch((657, y+141), (722, y+141),
            arrowstyle='-|>', mutation_scale=13, color=INK, lw=1.8, zorder=6))
        text(47, y+245, descriptions[number], 10.5, TEAL, 'bold')
        text(47, y+267, event_summaries[number], 9.5, MUTED)

    text(40, 760, 'Colors track atom origins; carbon numbers keep the same identity across both endpoints.', 10.3, MUTED)
    text(40, 790, 'H atoms are implicit. With their assignments optimized, either heavy-atom mapping has 6 total events.', 10.3, MUTED)
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(man/'figs'/f'fig4_alternatives.{extension}', facecolor='white',
                    bbox_inches='tight', pad_inches=.025, dpi=220)
    plt.close(fig)
    print('EnzymeMap alternatives figure generated')


if __name__ == '__main__':
    build(Path(__file__).resolve().parents[1])
