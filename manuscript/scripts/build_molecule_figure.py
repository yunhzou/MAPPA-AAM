"""Vector chemical worked example; formal bond weights, no benchmark search.

RDKit generates chemical depictions. Its vector paths are transferred to
Matplotlib without rasterization. Analytical checks validate the illustrated
partial mappings, intrinsic response symmetry, and heavy-atom event lists.
"""
from pathlib import Path
import os,json,re,xml.etree.ElementTree as ET
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parents[1]/"build/mpl-cache"))
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from matplotlib import pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch,Ellipse,FancyArrowPatch,FancyBboxPatch,Circle
from matplotlib.transforms import Affine2D

INK='#223340'; MUTED='#586974'; TEAL='#087F8C'; ORANGE='#CB702A'; PURPLE='#77559B'; RED='#B63F52'; BLUE='#3366A5'; LINE='#D5DEE2'
PALE_A=(.70,.88,.86); PALE_B=(1.,.82,.59); PALE_O=(.84,.76,.94)


def check_design_examples():
    """Direct pair counting for the molecular heavy-atom examples, not AAM runs."""
    from itertools import combinations, permutations

    def count(source, target, mapping, internal):
        events = []
        for a, b in combinations(sorted(mapping), 2):
            change = target.get(tuple(sorted((mapping[a], mapping[b]))), 0)-source.get((a, b), 0)
            if abs(change) > .5:
                events.append(dict(pair=[a, b], change=change,
                                   location='inner' if {a, b} <= internal else 'outer'))
        return events

    ring = {(0, 1): 1, (1, 2): 1, (1, 3): 1, (2, 3): 1}
    chain = {(0, 2): 1, (1, 2): 1, (2, 3): 1}
    ring_maps = {'center': {0: 0, 1: 2, 2: 1, 3: 3},
                 'end': {0: 0, 1: 1, 2: 2, 3: 3}}
    ring_events = {key: count(ring, chain, mapping, {1, 2, 3})
                   for key, mapping in ring_maps.items()}
    assert [len(ring_events[key]) for key in ('center', 'end')] == [1, 3]
    assert all(sum(e['location'] == 'inner' for e in events) == 1
               for events in ring_events.values())
    # This accepted path forces C1 onto a chain end, preventing growth to A.
    carbon_maps = [{0: 0, **dict(zip((1, 2, 3), perm))} for perm in permutations((1, 2, 3))]
    def present(mapping, edge):
        return chain.get(tuple(sorted(mapping[i] for i in edge)), 0) > .2
    path_matches = [mapping for mapping in carbon_maps
                    if all(present(mapping, edge) for edge in ((1, 2), (2, 3)))]
    assert path_matches and not any(present(mapping, (0, 1)) for mapping in path_matches)
    nitrile = {(1, 2): 1, (2, 3): 3}
    single = {(1, 2): 1, (2, 3): 1}
    nitrile_maps = {'retain': {1: 1, 2: 2, 3: 3}, 'exchange': {1: 2, 2: 1, 3: 3}}
    nitrile_events = {key: count(nitrile, single, mapping, {1, 2})
                      for key, mapping in nitrile_maps.items()}
    assert [len(nitrile_events[key]) for key in ('retain', 'exchange')] == [1, 2]
    assert all(e['location'] == 'outer' for events in nitrile_events.values() for e in events)
    assert abs(nitrile[(2, 3)]-single[(2, 3)]) > 1
    assert abs(nitrile[(2, 3)]-single[(2, 3)]) <= 2
    for smiles, expected in ((RING_ALCOHOL, ring), (CHAIN_END, chain),
                             (NITRILE, nitrile), (AMINE, single)):
        assert weighted(smiles)[1] == expected
    assert len(count(ring, weighted(CHAIN_CENTER)[1], {i: i for i in range(4)}, {1, 2, 3})) == 1
    assert len(count(nitrile, weighted(AMINE_EXCHANGED)[1], {i: i for i in (1, 2, 3)}, {1, 2})) == 2
    source_weights, target_weights = (weighted(smiles)[1] for smiles in (SWEEP_BASE, SWEEP_TARGET))
    assert source_weights[(3, 4)] == target_weights[(3, 4)] == 1
    assert all(source_weights[edge] == target_weights[edge] for edge in source_weights if min(edge) >= 3)
    return dict(scope='Analytical heavy-atom graph illustrations; no H constraints or benchmark reruns.',
                ring_events=ring_events, nitrile_events=nitrile_events,
                triple_to_single_weight_difference=2, completion=check_completion_example(),
                sweep=dict(base_smiles=SWEEP_BASE, target_smiles=SWEEP_TARGET,
                           masked_edge=[3, 4], highlighted_atoms=list(range(3, 10))))


RING_FRAGMENT = '[CH2:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1'
RING_SITES = ['[CH:1]1([F:7])[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1',
              '[CH:11]1([Cl:17])[CH2:12][CH2:13][CH2:14][CH2:15][CH2:16]1']
SWEEP_BASE = '[CH3:1][NH:2][CH2:3][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1'
SWEEP_TARGET = '[CH3:1][NH2:2].[CH3:3][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1'
RING_ALCOHOL = '[OH:0][CH:1]1[CH2:2][CH2:3]1'
CHAIN_CENTER = '[OH:0][CH:1]([CH3:2])[CH3:3]'
CHAIN_END = '[OH:0][CH:2]([CH3:1])[CH3:3]'
NITRILE = '[CH3:1][C:2]#[N:3]'
AMINE = '[CH3:1][CH2:2][NH2:3]'
AMINE_EXCHANGED = '[CH3:2][CH2:1][NH2:3]'


def check_completion_example():
    """Independent small-graph oracle for the two illustrated carbon placements."""
    from collections import Counter
    target = Chem.MolFromSmiles('.'.join(RING_SITES))
    query = Chem.MolFromSmarts('[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1')
    matches = target.GetSubstructMatches(query, uniquify=False, maxMatches=1000)
    sets = Counter(tuple(sorted(target.GetAtomWithIdx(i).GetAtomMapNum() for i in match))
                   for match in matches)
    assert sets == {tuple(range(1, 7)): 12, tuple(range(11, 17)): 12}
    # Different F/Cl attachments prevent a whole-target operation exchanging rings.
    carbon_sets = [{atom.GetIdx() for atom in target.GetAtoms() if atom.GetAtomMapNum() in labels}
                   for labels in (range(1, 7), range(11, 17))]
    for action in target.GetSubstructMatches(target, uniquify=False, maxMatches=1000):
        assert all({action[i] for i in group} == group for group in carbon_sets)
    return dict(scope='Fixed six-carbon fragment; hydrogen constraints omitted.',
                base_smiles=RING_FRAGMENT, target_smiles='.'.join(RING_SITES),
                placements=[dict(atoms=list(atoms), internal_maps=count) for atoms, count in sorted(sets.items())],
                placement_count=2, whole_target_ring_exchange=False)


def _chemical_canvas(height):
    fig = plt.figure(figsize=(11.8, height/120))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1416), ylim=(height, 0))
    ax.axis('off')
    return fig, ax


def _figure_text(ax, x, y, value, size=11, color=INK, weight='normal', align='left'):
    ax.text(x, y, value, fontsize=size, color=color, weight=weight,
            ha=align, va='center', linespacing=1.3, zorder=9)


def _figure_panel(ax, x, y, width, height, color='#F5F7FB', edge=None):
    ax.add_patch(FancyBboxPatch((x, y), width, height,
        boxstyle='round,pad=0,rounding_size=12', facecolor=color,
        edgecolor=edge or color, lw=.8, zorder=0))


def _figure_arrow(ax, start, end, color=INK):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>',
        mutation_scale=12, color=color, lw=1.6, zorder=6))


def _save_chemical_figure(fig, man, stem):
    for extension in ('pdf', 'png', 'svg'):
        fig.savefig(man/'figs'/f'{stem}.{extension}', dpi=240, facecolor='white',
                    bbox_inches='tight', pad_inches=.025)
    plt.close(fig)


def build_completion_example(man):
    check_completion_example()
    fig, ax = _chemical_canvas(620)
    _figure_panel(ax, 12, 12, 1392, 596)
    _figure_text(ax, 34, 44, 'MAPPA saturation: one ring, two placements', 15, weight='bold')
    _figure_text(ax, 145, 145, 'Base fragment', 11, weight='bold', align='center')
    _figure_text(ax, 545, 145, 'Target: two separate rings', 11, weight='bold', align='center')
    molecule(ax, RING_FRAGMENT, 30, 210, 235, 220,
             owners={i: PALE_A for i in range(1, 7)}, notes=False, align_acetyl=False)
    for x, smiles, label in zip((335, 560), RING_SITES, ('Ring A', 'Ring B')):
        molecule(ax, smiles, x, 215, 200, 220, notes=False, align_acetyl=False, font=26)
        _figure_text(ax, x+100, 458, label, 10.5, MUTED, align='center')
    _figure_text(ax, 550, 325, '+', 18, MUTED, align='center')
    _figure_arrow(ax, (265, 325), (322, 325))
    ax.plot([775, 850, 850], [325, 325, 208], color=INK, lw=1.6)
    ax.plot([850, 850], [325, 453], color=INK, lw=1.6)
    for number, y in enumerate((95, 340)):
        _figure_panel(ax, 910, y, 478, 225, 'white', LINE)
        _figure_arrow(ax, (850, y+113), (902, y+113))
        _figure_text(ax, 930, y+25, f'Branch {number+1}: select ring {"AB"[number]}', 11.5, TEAL, 'bold')
        for site, (x, smiles) in enumerate(zip((925, 1160), RING_SITES)):
            atoms = range(1, 7) if site == 0 else range(11, 17)
            owners = {i: PALE_A for i in atoms} if site == number else None
            molecule(ax, smiles, x, y+49, 205, 151, owners=owners,
                     notes=False, align_acetyl=False, font=24)
        _figure_text(ax, 1145, y+124, '+', 16, MUTED, align='center')
    _figure_text(ax, 42, 522, 'Same base atoms', 10.5, MUTED)
    _figure_text(ax, 363, 522, 'Different target atom sets', 10.5, MUTED)
    _figure_text(ax, 34, 583, 'Rotations and reflections stay inside each placement; they do not create extra branches.', 10.5, TEAL)
    _save_chemical_figure(fig, man, 'fig11_saturation')


def build_design_examples(man):
    """Molecular depictions of the checked heavy-atom design examples."""
    proof = check_design_examples()
    (man/'evidence/design_examples.json').write_text(json.dumps(proof, indent=2)+'\n')
    fig, ax = _chemical_canvas(1160)
    text = lambda x, y, value, **kw: _figure_text(ax, x, y, value, **kw)
    for y, height in ((12, 278), (302, 450), (764, 384)):
        _figure_panel(ax, 12, y, 1392, height)

    text(34, 42, 'a   Greedy constraints and the single-edge sweep', size=14, weight='bold')
    for x, label in ((222, 'Base'), (705, 'Target'), (1185, 'One sweep condition')):
        text(x, 84, label, size=11, weight='bold', align='center')
    owners = {i: PALE_B for i in range(3, 10)}
    molecule(ax, SWEEP_BASE, 40, 99, 370, 144, owners=owners, notes=False,
             align_acetyl=False, font=24)
    _figure_arrow(ax, (424, 173), (474, 173))
    molecule(ax, '[CH3:1][NH2:2]', 492, 118, 123, 105, notes=False,
             align_acetyl=False, font=24)
    text(633, 173, '+', size=17, align='center')
    molecule(ax, SWEEP_TARGET.split('.')[1], 650, 99, 250, 144,
             owners=owners, notes=False, align_acetyl=False, font=24)
    positions = molecule(ax, SWEEP_BASE, 977, 99, 370, 144, owners=owners,
                         notes=False, align_acetyl=False, font=24)
    a, b = positions[3], positions[4]
    ax.plot(*zip(a, b), color='white', lw=4, zorder=6)
    ax.plot(*zip(a, b), color=PURPLE, lw=1.8, ls=(0, (3, 2)), zorder=7)
    middle = ((a[0]+b[0])/2, (a[1]+b[1])/2)
    ax.plot([middle[0]-5, middle[0]+5], [middle[1]+10, middle[1]-10], color=PURPLE, lw=2, zorder=8)
    text(222, 263, 'Large benzyl region highlighted', size=10, color=MUTED, align='center')
    text(705, 263, 'A large region can be preserved', size=10, color=MUTED, align='center')
    text(1185, 263, 'Mask one edge during growth', size=10, color=PURPLE, align='center')

    text(34, 332, 'b   Ring opening can couple an inner shuffle to an outer bond', size=14, weight='bold')
    for x, label in ((222, 'Base: cyclopropanol'), (715, 'Target: propan-2-ol'), (1180, 'Heavy-atom events')):
        text(x, 375, label, size=11, weight='bold', align='center')
    ring_positions = {0: (0, 2.5), 1: (0, 1), 2: (-.866, -.5), 3: (.866, -.5)}
    owners = {0: PALE_B, 1: PALE_A, 2: PALE_A, 3: PALE_A}
    for y, target, choices, summary in ((391, CHAIN_CENTER, 'C₁ stays next to O', '1 inner + 0 outer = 1 event'),
                                       (558, CHAIN_END, 'C₁ moves to a chain end', '1 inner + 2 outer = 3 events')):
        center = target == CHAIN_CENTER
        cuts = {(2, 3): -1} if center else {(1, 3): -1, (0, 1): -1}
        molecule(ax, RING_ALCOHOL, 50, y, 345, 151, owners=owners,
                 notes={i: str(i) for i in (1, 2, 3)}, font=24, align_acetyl=False,
                 coordinates=ring_positions, bond_changes=cuts)
        _figure_arrow(ax, (425, y+77), (500, y+77))
        coords = ({0: (0, 1.5), 1: (0, 0), 2: (-1.3, -.75), 3: (1.3, -.75)} if center
                  else {0: (0, 1.5), 2: (0, 0), 1: (-1.3, -.75), 3: (1.3, -.75)})
        molecule(ax, target, 525, y, 380, 151, owners=owners,
                 notes={i: str(i) for i in (1, 2, 3)}, font=24, align_acetyl=False,
                 coordinates=coords, bond_changes={} if center else {(0, 2): 1})
        text(1170, y+60, choices, size=11, weight='bold', align='center')
        text(1170, y+96, summary, size=10.5, color=TEAL if center else RED, align='center')
    text(34, 733, 'The carbon fragment is the same size; the oxygen attachment can change with its assignment.', size=10, color=MUTED)

    text(34, 794, 'c   A permissive tolerance can keep the changing bond inside a fragment', size=14, weight='bold')
    for x, label in ((222, 'Base'), (705, 'Target: retain carbons'), (1185, 'Target: exchange carbons')):
        text(x, 837, label, size=11, weight='bold', align='center')
    owners = {1: PALE_A, 2: PALE_A, 3: PALE_B}
    coords = {1: (-1.5, 0), 2: (0, 0), 3: (1.5, 0)}
    molecule(ax, NITRILE, 45, 850, 350, 121, owners=owners,
             notes={1: '1', 2: '2'}, font=24, align_acetyl=False, coordinates=coords,
             bond_changes={(2, 3): -1})
    _figure_arrow(ax, (422, 910), (480, 910))
    for x, smiles, exchanged in ((530, AMINE, False), (1005, AMINE_EXCHANGED, True)):
        coords = {1: (-1.3, -.75), 2: (0, 0), 3: (1.3, -.75)}
        if exchanged: coords[1], coords[2] = coords[2], coords[1]
        molecule(ax, smiles, x, 850, 350, 121, owners=owners,
                 notes={1: '1', 2: '2'}, font=24, align_acetyl=False, coordinates=coords,
                 bond_changes={(1, 3): 1} if exchanged else {})
    text(222, 993, 'τ = 1: C–C and N separate', size=10, color=MUTED, align='center')
    text(705, 993, '1 outer weakening', size=10.5, color=TEAL, align='center')
    text(1185, 993, '2 outer events', size=10.5, color=RED, align='center')
    molecule(ax, NITRILE, 45, 1023, 190, 88, owners={i: PALE_A for i in (1, 2, 3)},
             notes=False, font=20, align_acetyl=False)
    _figure_arrow(ax, (250, 1065), (301, 1065))
    molecule(ax, AMINE, 320, 1023, 190, 88, owners={i: PALE_A for i in (1, 2, 3)},
             notes=False, font=20, align_acetyl=False)
    text(565, 1063, 'τ = 2: C–N can stay inside one fragment', size=11, color=TEAL, weight='bold')
    text(34, 1126, 'Heavy-atom illustrations; H constraints omitted. Red: break / weaken. Green: form / strengthen.', size=9.5, color=MUTED)
    _save_chemical_figure(fig, man, 'fig8_design_examples')

ESTER='[CH3:1][C:2](=[O:3])[O:4][CH3:5]'
HYDROXIDE='[OH-:6]'
ACETATE='[CH3:1][C:2](=[O:3])[O-:4]'
ALCOHOL='[CH3:5][OH:6]'
R=ESTER+'.'+HYDROXIDE; P=ACETATE+'.'+ALCOHOL
ALCOHOL_FRAGMENT='[*:14][CH:12]([OH:13])[CH3:11]'
ALCOHOL_SITES='[N:50]([CH:22]([OH:23])[CH3:21])([CH2:51][CH:32]([OH:33])[CH3:31])[CH2:52][CH2:53][CH:42]([OH:43])[CH3:41]'

def weighted(smiles):
 m=Chem.MolFromSmiles(smiles); assert m is not None
 z={a.GetAtomMapNum():a.GetAtomicNum() for a in m.GetAtoms()}
 w={tuple(sorted((b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()))):b.GetBondTypeAsDouble() for b in m.GetBonds()}
 return z,w

def check_example():
 zr,wr=weighted(R);zp,wp=weighted(P)
 identity={i:i for i in zr};exchange={**identity,3:4,4:3}
 def events(m):
  assert sorted(m.values())==sorted(zp) and all(zr[i]==zp[j] for i,j in m.items())
  return [[i,j,1 if d>0 else -1] for i in zr for j in zr if i<j and abs(d:=wp.get(tuple(sorted((m[i],m[j]))),0)-wr.get((i,j),0))>.5]
 def preserves(m,fragments):
  return all(wp.get(tuple(sorted((m[i],m[j]))),0)>=.2 and abs(wp.get(tuple(sorted((m[i],m[j]))),0)-v)<=1 for f in fragments for (i,j),v in wr.items() if i in f and j in f)
 assert preserves(identity,[[1,2,3,4],[5],[6]])
 assert preserves(exchange,[[1,2,3,4],[5],[6]])
 # The accepted acetoxy placement cannot grow through O4-C5 in either case.
 assert all(wp.get(tuple(sorted((m[4],m[5]))),0)<.2 for m in [identity,exchange])
 # The O3/O4 exchange is an exact symmetry of the matching-response graph,
 # although it is not a symmetry of the formal bond-order graph.
 response=lambda w:[w>=.2,*[w>=.2 and abs(w-v)<=1 for v in [1,2]]]
 assert response(1)==response(2)==[True,True,True]
 assert events(identity)==[[4,5,-1],[5,6,1]]
 assert events(exchange)==[[2,3,-1],[2,4,1],[4,5,-1],[5,6,1]]
 # Validate complete local alcohol fragments, including every attached H.
 source=Chem.AddHs(Chem.MolFromSmiles(ALCOHOL_FRAGMENT))
 target=Chem.AddHs(Chem.MolFromSmiles(ALCOHOL_SITES))
 si={a.GetAtomMapNum():a.GetIdx() for a in source.GetAtoms() if a.GetAtomMapNum()}
 ti={a.GetAtomMapNum():a.GetIdx() for a in target.GetAtoms() if a.GetAtomMapNum()}
 placements={}
 for site,base in [('a',20),('b',30),('c',40)]:
  mapping={si[i]:ti[base+i-10] for i in [11,12,13]}
  for i in [11,12,13]:
   sh=[a.GetIdx() for a in source.GetAtomWithIdx(si[i]).GetNeighbors() if a.GetAtomicNum()==1]
   th=[a.GetIdx() for a in target.GetAtomWithIdx(mapping[si[i]]).GetNeighbors() if a.GetAtomicNum()==1]
   assert len(sh)==len(th)
   mapping.update(zip(sh,th))
  assert len(set(mapping.values()))==len(mapping)
  for bond in source.GetBonds():
   i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
   if i in mapping and j in mapping:
    other=target.GetBondBetweenAtoms(mapping[i],mapping[j]);assert other and other.GetBondTypeAsDouble()==bond.GetBondTypeAsDouble()
  placements[site]={'heavy_atoms':[base+1,base+2,base+3],'hydrogen_inclusive_mapping':sorted(mapping.items())}
 # Sites have different whole-molecule environments, not interchangeable copies.
 plain=Chem.MolFromSmiles(ALCOHOL_SITES)
 for atom in plain.GetAtoms():atom.SetAtomMapNum(0)
 ranks=Chem.CanonicalRankAtoms(plain,breakTies=False)
 assert len({ranks[ti[k]] for k in [22,32,42]})==3
 tree=[]
 for a in placements:
  for b in placements:
   if a==b:continue
   c=next(s for s in placements if s not in [a,b]);tree.append({'A':a,'B':b,'C':c})
 assert len(tree)==6 and all(len(set(row.values()))==3 for row in tree)
 # Extend the two illustrated family witnesses to all explicit hydrogens.
 reactant=Chem.AddHs(Chem.MolFromSmiles(R));product=Chem.AddHs(Chem.MolFromSmiles(P))
 ri={a.GetAtomMapNum():a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomMapNum()}
 pi={a.GetAtomMapNum():a.GetIdx() for a in product.GetAtoms() if a.GetAtomMapNum()}
 for heavy in [identity,exchange]:
  mapping={ri[i]:pi[j] for i,j in heavy.items()}
  for i,j in heavy.items():
   hs=[a.GetIdx() for a in reactant.GetAtomWithIdx(ri[i]).GetNeighbors() if a.GetAtomicNum()==1]
   ht=[a.GetIdx() for a in product.GetAtomWithIdx(pi[j]).GetNeighbors() if a.GetAtomicNum()==1]
   assert len(hs)==len(ht);mapping.update(zip(hs,ht))
  assert len(mapping)==reactant.GetNumAtoms()==len(set(mapping.values()))==product.GetNumAtoms()
  for bond in reactant.GetBonds():
   i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
   if reactant.GetAtomWithIdx(i).GetAtomicNum()==1 or reactant.GetAtomWithIdx(j).GetAtomicNum()==1:
    target_bond=product.GetBondBetweenAtoms(mapping[i],mapping[j]);assert target_bond and target_bond.GetBondTypeAsDouble()==bond.GetBondTypeAsDouble()
 # In panel c both acetate oxygens have zero attached H: unlike an OH/ether
 # exchange, the displayed oxygen action does not hide a lost O-H bond.
 acetate=Chem.MolFromSmiles(ACETATE)
 assert all(a.GetTotalNumHs()==0 for a in acetate.GetAtoms() if a.GetAtomicNum()==8)
 return {'placement_example':{'source_fragment_smiles':ALCOHOL_FRAGMENT,'target_smiles':ALCOHOL_SITES,'admissible_local_images':placements,'conditional_tree':tree,'scope':'Three distinct local alcohol sites; all attached H checked. R is outside the matched fragment. Constructed conditional placement illustration, not a complete reaction-search trajectory.'},'status':'passed','scope':'Constructed molecular illustration with formal bond orders. Panel a verifies attached H explicitly; panels b/c display heavy-atom events of ester saponification. Saved constraints govern full-family actions.',
 'reactants_smiles':R,'products_smiles':P,'iso_tol':1.,'edge_floor':.2,'event_threshold':.5,
 'mappings':{'retain_A':identity,'oxygen_exchange':exchange},
 'heavy_atom_events':{'retain_A':events(identity),'oxygen_exchange':events(exchange)},
 'response_labels':{'single':response(1),'double':response(2)},'saturated_fragment':[1,2,3,4]}

def svg_path(d):
 tokens=re.findall(r'[MLQCZmlqcz]|[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?',d)
 vertices=[];codes=[];i=0;start=(0,0)
 while i<len(tokens):
  c=tokens[i];i+=1
  assert c in {'M','L','Q','C','Z'},c
  if c=='Z':vertices.append(start);codes.append(MPath.CLOSEPOLY);continue
  n={'M':2,'L':2,'Q':4,'C':6}[c]
  a=list(map(float,tokens[i:i+n]));i+=n
  pts=list(zip(a[::2],a[1::2]));vertices.extend(pts)
  codes.extend([{'M':MPath.MOVETO,'L':MPath.LINETO,'Q':MPath.CURVE3,'C':MPath.CURVE4}[c]]*len(pts))
  if c=='M':start=pts[0]
 return MPath(vertices,codes)

EVENT_RED = '#EF1026'
EVENT_GREEN = '#00B83F'

def bond_event_marker(ax, a, b, sign, size=None, y_scale=1):
 """Small, vivid marks centered on a black bond; no filled badge or leader.

 Red crossing strokes indicate loss; two green inward arrows indicate gain.
 The bond-axis coordinates keep the glyph attached to slanted bonds as well.
 """
 import math
 from matplotlib import patheffects
 a,b=(a[0],a[1]/y_scale),(b[0],b[1]/y_scale)
 length = math.dist(a, b)
 if not length: return
 ux, uy = (b[0]-a[0])/length, (b[1]-a[1])/length
 mid = ((a[0]+b[0])/2, (a[1]+b[1])/2)
 size = size if size is not None else min(8.5, max(5.0, length*.15))
 def point(x,y): return (mid[0]+ux*x-uy*y, mid[1]+uy*x+ux*y)
 def stroke(points):
  q=[(px,py*y_scale) for px,py in [point(x,y) for x,y in points]]
  ax.plot(*zip(*q),color=EVENT_RED if sign<0 else EVENT_GREEN,lw=1.8,
          solid_capstyle='round',solid_joinstyle='round',zorder=8,
          path_effects=[patheffects.Stroke(linewidth=3.1,foreground='white'),patheffects.Normal()])
 if sign<0:
  ux,uy=1.,0.
  stroke([(-size*.63,-size*.63),(size*.63,size*.63)])
  stroke([(-size*.63,size*.63),(size*.63,-size*.63)])
 else:
  # Opposed arrowheads point toward the new connection, with a visible center gap.
  for direction in [-1,1]:
   stroke([(direction*size*1.3,0),(direction*size*.20,0)])
   stroke([(direction*size*.76,-size*.54),(direction*size*.20,0),(direction*size*.76,size*.54)])

def molecule(ax,smiles,x,y,w,h,owners=None,notes=None,font=24,explicit_hydrogens=False,align_acetyl=True,rotation=0,bond_changes=None,preserve_bond_orders=False,coordinates=None,atom_labels=None):
 params=Chem.SmilesParserParams();params.removeHs=not explicit_hydrogens;params.sanitize=not preserve_bond_orders
 m=Chem.MolFromSmiles(smiles,params);assert m is not None
 if preserve_bond_orders:Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
 ids={a.GetAtomMapNum():a.GetIdx() for a in m.GetAtoms()}
 for a in m.GetAtoms():
  k=a.GetAtomMapNum()
  if notes is None:a.SetProp('atomNote',str(k))
  elif notes is not False and k in notes:a.SetProp('atomNote',str(notes[k]))
  a.SetAtomMapNum(0)
 rdDepictor.Compute2DCoords(m)
 if coordinates:
  for k,xy in coordinates.items():
   if k in ids:m.GetConformer().SetAtomPosition(ids[k],(*xy,0))
 # Shared acetyl coordinates make the ester and acid directly comparable.
 if align_acetyl and {1,2,3,4}.issubset(ids):
  points={1:(-1.30,.75),2:(0,0),3:(0,-1.5),4:(1.30,.75),5:(2.60,0)}
  for k,i in ids.items():
   if k in points:m.GetConformer().SetAtomPosition(i,(*points[k],0))
 d=rdMolDraw2D.MolDraw2DSVG(int(w),int(h));o=d.drawOptions();o.clearBackground=False;o.rotate=rotation;o.useBWAtomPalette();o.fixedFontSize=font;o.annotationFontScale=.62;o.padding=.11;o.bondLineWidth=1.8
 colors={ids[k]:v for k,v in (owners or {}).items() if k in ids}
 bonds={b.GetIdx():colors[b.GetBeginAtomIdx()] for b in m.GetBonds() if b.GetBeginAtomIdx() in colors and colors.get(b.GetBeginAtomIdx())==colors.get(b.GetEndAtomIdx())}
 for atom in m.GetAtoms():
  if atom.GetAtomicNum()==0:o.atomLabels[atom.GetIdx()]='R'
 for label,value in (atom_labels or {}).items():
  if label in ids:o.atomLabels[ids[label]]=value
 for pair,sign in (bond_changes or {}).items():
  bond=m.GetBondBetweenAtoms(ids[pair[0]],ids[pair[1]])
  assert bond is not None, ('Changed bond is absent from depiction',pair)
 d.DrawMolecule(m,highlightAtoms=list(colors),highlightBonds=list(bonds),highlightAtomColors=colors,highlightBondColors=bonds)
 d.FinishDrawing();root=ET.fromstring(d.GetDrawingText())
 transform=Affine2D().translate(x,y)+ax.transData
 for e in root:
  kind=e.tag.split('}')[-1]
  if kind not in {'path','ellipse','rect'}:continue
  style=dict(item.split(':',1) for item in e.get('style','').split(';') if ':' in item)
  fill=style.get('fill',e.get('fill','none'));stroke=style.get('stroke',e.get('stroke','none'))
  kw=dict(facecolor=fill,edgecolor=stroke,linewidth=float(style.get('stroke-width','0').removesuffix('px'))*.60,transform=transform,zorder=3)
  if fill=='none' or fill in {'#000000','black'}:
   kw['zorder']=5
  if kind=='path':patch=PathPatch(svg_path(e.attrib['d']),**kw)
  elif kind=='ellipse':patch=Ellipse((float(e.attrib['cx']),float(e.attrib['cy'])),2*float(e.attrib['rx']),2*float(e.attrib['ry']),**kw)
  else:raise AssertionError('Unexpected background rectangle in molecule')
  ax.add_patch(patch)
 positions={k:(x+d.GetDrawCoords(i).x,y+d.GetDrawCoords(i).y) for k,i in ids.items()}
 # Preserve the black endpoint bond and overlay only a compact event glyph.
 for pair,sign in (bond_changes or {}).items():
  bond_event_marker(ax,positions[pair[0]],positions[pair[1]],sign)
 return positions


def build(man):
 proof=check_example();(man/'evidence/molecule_example.json').write_text(json.dumps(proof,indent=2)+'\n')
 fig=plt.figure(figsize=(10.8,8.15));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1440),ylim=(1087,0));ax.axis('off')
 navy='#22324D';muted='#64748B';faint='#BCC7D5';bg='#F5F7FB';violet='#7957BA';orange='#CC762F';green='#128477'
 def t(x,y,s,size=11,color=navy,weight='normal',ha='left'):
  return ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.3,zorder=6)
 def box(x,y,w,h,fc='white',ec=LINE,lw=.8,r=12):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',facecolor=fc,edgecolor=ec,lw=lw,zorder=1))
 def ar(a,b,color=navy,lw=1.9):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=11,lw=lw,color=color,zorder=2))
 def tree(points,color=navy,lw=1.7):
  xs,ys=zip(*points[:-1]);ax.plot(xs,ys,color=color,lw=lw,solid_joinstyle='round',solid_capstyle='round',zorder=2)
  ar(points[-2],points[-1],color,lw)
 def dot(x,y,color=navy,r=3):ax.add_patch(Circle((x,y),r,fc=color,ec='none',zorder=3))
 def label(x,y,letter,title):
  t(x,y,letter,14,navy,'bold');t(x+36,y,title,14,navy,'bold')

 # The local motif and open placement tree distinguish witness growth from
 # completion. The tree's six paths are the checked conditional assignments.
 box(12,12,1416,545,bg,bg,r=16)
 label(34,46,'a','MAPPA matching: grow, complete and branch')
 t(34,79,'A, B and C are successive local alcohol fragments; R lies outside the matched motif.',10.5,muted)
 t(236,116,'Grow on base; keep one witness',11,green,'bold',ha='center')
 for x,atoms,name in [(24,[11],'Seed'),(170,[11,12],'Extend'),(316,[11,12,13],'Include O-H')]:
  molecule(ax,ALCOHOL_FRAGMENT,x,145,132,116,{i:PALE_A for i in atoms},notes=False,font=28)
  t(x+66,278,name,10,green,ha='center')
 ar((151,203),(167,203),green,1.4);ar((297,203),(313,203),green,1.4)
 t(236,320,'Target: three distinct sites',10.5,navy,'bold',ha='center')
 sites=molecule(ax,ALCOHOL_SITES,33,338,405,166,{i:PALE_A for i in [21,22,23,31,32,33,41,42,43]},notes=False,font=26)
 for atom,name,dx,dy in [(23,'a',30,-24),(33,'b',-30,24),(43,'c',-25,-25)]:
  x,y=sites[atom];t(x+dx,y+dy,name,11,green,'bold',ha='center')
 t(236,532,'C-C, C-O and every attached H are checked.',9.5,muted,ha='center')
 t(654,116,'At saturation',11,green,'bold',ha='center')
 t(654,145,'Complete A',10,green,ha='center')
 t(950,116,'Complete B',11,orange,'bold',ha='center')
 t(1263,116,'Complete C',11,violet,'bold',ha='center')
 t(950,145,'A remains fixed',9.5,muted,ha='center')
 t(1263,145,'A and B remain fixed',9.5,muted,ha='center')
 ar((446,203),(510,203),green,1.6)
 for site,y,remaining in [('a',203,['b','c']),('b',335,['a','c']),('c',467,['a','b'])]:
  tree([(510,203),(546,203),(546,y),(592,y)],green)
  t(654,y,'A → '+site,14,green,'bold',ha='center')
  for dest,yy in zip(remaining,[y-34,y+34]):
   last=next(k for k in remaining if k!=dest)
   tree([(716,y),(796,y),(796,yy),(870,yy)],navy,1.5)
   t(950,yy,'B → '+dest,13,orange,'bold',ha='center')
   ar((1023,yy),(1173,yy),navy,1.5)
   t(1263,yy,'C → '+last,13,violet,'bold',ha='center')
  dot(796,y)
 t(977,532,'Three placements → two unused sites → one remaining site',10,muted,ha='center')

 # A retained branch is a composition of paired endpoint atom sets. Show the
 # full endpoints so the inter-fragment changes have visible before/after bonds.
 box(12,576,708,454,bg,bg,r=16)
 label(34,612,'b','Matching output: paired fragments')
 t(34,646,'The same color pairs atoms across the two endpoints.',10.5,muted)
 owners={**{i:PALE_A for i in [1,2,3,4]},5:PALE_B,6:PALE_O}
 for y,title in [(688,'Base'),(830,'Target')]:
  t(40,y+53,title,11,navy,'bold')
  first=ESTER if title=='Base' else ACETATE
  second=HYDROXIDE if title=='Base' else ALCOHOL
  molecule(ax,first,143,y,290,111,owners,notes={2:'2',3:'3',4:'4',5:'5'},font=28,
           bond_changes={(4,5):-1} if title=='Base' else None)
  t(455,y+53,'+',13,muted,ha='center')
  molecule(ax,second,485,y+15,175,87,owners,notes={5:'5',6:'6'},font=28,
           bond_changes={(5,6):1} if title=='Target' else None)
 ar((361,802),(361,826),navy,1.6)
 t(40,966,'One branch:',10.5,navy,'bold')
 for x,name,color in [(200,'F',PALE_A),(320,'G',PALE_B),(440,'H',PALE_O)]:
  dot(x,966,color,8);t(x+17,966,name,11,navy,'bold')
  if name!='H':t(x+73,966,'+',11,muted)
 t(40,993,'2 boundary events: break 4-5; form 5-6.',10,navy)
 t(40,1015,'Share repeated fragment pairs across seeds and cuts.',9.5,muted)

 # Every local response shows both halves of F. The base ester retains its
 # unhighlighted methyl neighbor as context; that neighbor is not part of F.
 box(738,576,690,454,bg,bg,r=16)
 label(760,612,'c','Decoding: compare both sides of F')
 t(760,646,'Binary-connectivity automorphism permits the O3 / O4 swap.',10.5,muted)
 t(878,685,'Base side of F',10.5,green,'bold',ha='center')
 t(1145,685,'Target side of F',10.5,green,'bold',ha='center')
 t(1340,685,'Branch total',10.5,navy,'bold',ha='center')
 for y,exchange,count,color in [(724,False,2,green),(870,True,4,RED)]:
  t(1008,y,'Exchange O identities' if exchange else 'Keep O identities',10,violet,ha='center')
  molecule(ax,ESTER,754,y+15,247,112,{i:PALE_A for i in [1,2,3,4]},
           notes={2:'2',3:'3',4:'4',5:'5'},font=28,
           bond_changes={(2,3):-1} if exchange else None)
  ar((1003,y+71),(1040,y+71),navy,1.5)
  molecule(ax,ACETATE,1044,y+15,200,112,{i:PALE_A for i in [1,2,3,4]},
           notes={2:'2',3:'4' if exchange else '3',4:'3' if exchange else '4'},font=28,
           bond_changes={(2,3):1} if exchange else None)
  t(1340,y+48,str(count),21,color,'bold',ha='center')
  t(1340,y+78,'events',10,color,ha='center')
  t(1340,y+111,f'2 boundary\n+ {count-2} internal',9.5,muted,ha='center')
 t(760,1015,'Only highlighted atoms belong to F; C5 is outside F.',9.5,muted)

 # Preserve the original chemical/event legend; search cuts are not events.
 ar((35,1064),(84,1064),navy,1.7);t(95,1064,'Matching / correspondence',10,muted)
 t(490,1064,'Numbers follow base atom identities',10,muted)
 for x,sign,text in [(886,-1,'Break / weaken'),(1155,1,'Form / strengthen')]:
  ax.plot([x,x+50],[1064,1064],color='black',lw=1.1)
  bond_event_marker(ax,(x,1064),(x+50,1064),sign,size=8)
  t(x+61,1064,text,10,muted)
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 for text in ax.texts:
  extent=text.get_window_extent(renderer)
  assert extent.x0>=-1 and extent.x1<=fig.bbox.width+1 and extent.y0>=-1 and extent.y1<=fig.bbox.height+1,text.get_text()
 for ext in ['pdf','svg','png']:fig.savefig(man/'figs'/f'fig1_algorithm.{ext}',facecolor='white',dpi=260)
 plt.close(fig)
 print('Molecular overview: witness growth, conditional placements, paired fragments and local events.')

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
