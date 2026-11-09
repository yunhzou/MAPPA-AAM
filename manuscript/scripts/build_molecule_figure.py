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
    """Direct pair counting for the schematic heavy-atom examples, not AAM runs."""
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
    return dict(scope='Analytical heavy-atom graph illustrations; no H constraints or benchmark reruns.',
                ring_events=ring_events, nitrile_events=nitrile_events,
                triple_to_single_weight_difference=2)


def build_design_examples(man):
    """Use the manuscript's palette and bond-event glyphs for design examples."""
    import numpy as np
    proof = check_design_examples()
    (man/'evidence/design_examples.json').write_text(json.dumps(proof, indent=2)+'\n')
    fig, axes = plt.subplots(3, 1, figsize=(10.4, 8.4))
    fig.subplots_adjust(left=.035, right=.98, top=.98, bottom=.03, hspace=.22)
    for ax in axes:
        ax.set(xlim=(0, 112), ylim=(100, 0))
        ax.axis('off')

    def text(ax, x, y, label, size=9, color=INK, weight='normal', align='center'):
        ax.text(x, y, label, ha=align, va='center', fontsize=size, color=color,
                weight=weight, linespacing=1.25, zorder=10)

    def arrow(ax, a, b, dashed=False, curve=0):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle='->', mutation_scale=12,
                                    color=MUTED, lw=1.2, linestyle='--' if dashed else '-',
                                    connectionstyle=f'arc3,rad={curve}', zorder=2))

    def region(ax, x, y, width, height, color=PALE_A):
        ax.add_patch(Ellipse((x, y), width, height, color=color, alpha=.65, lw=0, zorder=0))

    def event(ax, a, b, sign):
        scale = (ax.bbox.width/112)/(ax.bbox.height/100)
        bond_event_marker(ax, a, b, sign, size=.95, y_scale=scale)

    def graph(ax, positions, edges, labels, colors=None, changes=None, masked=()):
        for edge, order in edges.items():
            a, b = (np.array(positions[i], dtype=float) for i in edge)
            normal = np.array([-(b-a)[1], (b-a)[0]])
            normal /= np.linalg.norm(normal)
            for offset in np.arange(order)-(order-1)/2:
                start, end = a+normal*offset*1.25, b+normal*offset*1.25
                ax.plot(*zip(start, end), color=INK, lw=1.5,
                        linestyle='--' if edge in masked else '-', zorder=2)
            if edge in (changes or {}):
                event(ax, a, b, changes[edge])
        for atom, (x, y) in positions.items():
            ax.scatter([x], [y], s=390, color=(colors or {}).get(atom, 'white'),
                       edgecolors='none', zorder=4)
            text(ax, x, y, labels[atom], 10)

    a = axes[0]
    text(a, 0, 5, 'a  Greedy constraints and the single-edge sweep', 11, weight='bold', align='left')
    text(a, 19, 23, 'Base', weight='bold')
    text(a, 59, 23, 'Target', weight='bold')
    source = {0: (10, 53), 1: (23, 39), 2: (23, 68), 3: (36, 68)}
    labels = {0: 'R', 1: r'B$_1$', 2: r'B$_2$', 3: 'A'}
    region(a, 29.5, 68, 27, 25, PALE_B)
    graph(a, source, {(0, 1): 1, (0, 2): 1, (2, 3): 1}, labels)
    arrow(a, (40, 53), (46, 53))
    target = {0: (50, 39), 1: (65, 39), 2: (50, 68), 3: (65, 68)}
    region(a, 57.5, 68, 27, 25, PALE_B)
    graph(a, target, {(0, 1): 1, (2, 3): 1}, labels)
    text(a, 25, 91, 'Large B–A region: many seed sites', 8.5, MUTED)
    text(a, 59, 91, 'A repeated preserved unit', 8.5, MUTED)
    text(a, 95, 23, 'One sweep condition', weight='bold')
    cut = {0: (80, 47), 1: (94, 35), 2: (94, 65), 3: (108, 65)}
    graph(a, cut, {(0, 1): 1, (0, 2): 1, (2, 3): 1}, labels, masked=((2, 3),))
    text(a, 101, 65, '/', 18, PURPLE, weight='bold')
    text(a, 95, 87, 'Mask a preservation edge;\nlet growth choose a new constraint', 8.3, MUTED)

    a = axes[1]
    text(a, 0, 4, 'b  Different connectivity can couple an outer bond to an inner shuffle', 11, weight='bold', align='left')
    text(a, 18, 23, 'A–C₃ ring', weight='bold')
    text(a, 57, 23, 'C₁ at the chain center', weight='bold')
    text(a, 94, 23, 'C₁ at a chain end', weight='bold')
    source = {0: (18, 39), 1: (18, 57), 2: (8, 76), 3: (28, 76)}
    labels = {0: 'A', 1: r'C$_1$', 2: r'C$_2$', 3: r'C$_3$'}
    region(a, 18, 69, 34, 36)
    graph(a, source, {(0, 1): 1, (1, 2): 1, (1, 3): 1, (2, 3): 1}, labels, {0: PALE_B})
    center = {0: (57, 39), 1: (57, 68), 2: (45, 68), 3: (69, 68)}
    region(a, 57, 68, 36, 24)
    graph(a, center, {(0, 1): 1, (1, 2): 1, (1, 3): 1}, labels, {0: PALE_B})
    end = {0: (94, 39), 1: (82, 68), 2: (94, 68), 3: (106, 68)}
    region(a, 94, 68, 36, 24)
    graph(a, end, {(0, 2): 1, (1, 2): 1, (2, 3): 1}, labels, {0: PALE_B}, {(0, 2): 1})
    a.plot([94, 82], [39, 68], ls='--', color=RED, lw=1.2, zorder=1)
    event(a, (94, 39), (82, 68), -1)
    arrow(a, (22, 57), (52, 65), dashed=True, curve=-.12)
    arrow(a, (59, 72), (79, 72), dashed=True, curve=.3)
    text(a, 18, 96, 'Final carbon fragment highlighted', 8.3, MUTED)
    text(a, 57, 96, '1 inner + 0 outer = 1 event', 9, TEAL, 'bold')
    text(a, 94, 96, '1 inner + 2 outer = 3 events', 9, RED, 'bold')

    a = axes[2]
    text(a, 0, 4, 'c  Permissive growth avoids a tolerance-created boundary', 11, weight='bold', align='left')
    text(a, 19, 23, 'Base', weight='bold')
    text(a, 59, 23, 'Target: retain carbons', weight='bold')
    text(a, 96, 23, 'Target: exchange carbons', weight='bold')
    labels = {1: r'C$_1$', 2: r'C$_2$', 3: 'N'}
    source = {1: (7, 45), 2: (20, 45), 3: (33, 45)}
    target = {1: (45, 45), 2: (58, 45), 3: (71, 45)}
    swapped = {1: (96, 45), 2: (83, 45), 3: (109, 45)}
    for x in (13.5, 51.5, 89.5):
        region(a, x, 45, 25, 24)
    graph(a, source, {(1, 2): 1, (2, 3): 3}, labels, {3: PALE_B})
    graph(a, target, {(1, 2): 1, (2, 3): 1}, labels, {3: PALE_B}, {(2, 3): -1})
    graph(a, swapped, {(1, 2): 1, (1, 3): 1}, labels, {3: PALE_B}, {(1, 3): 1})
    a.plot([83, 96, 109], [43, 32, 43], color=RED, ls='--', lw=1.1, zorder=1)
    event(a, (93, 32), (99, 32), -1)
    arrow(a, (37, 45), (41, 45))
    text(a, 19, 66, 'τ = 1: C–C and N separate', 8.5, MUTED)
    text(a, 58, 66, '1 outer weakening', 9, TEAL, 'bold')
    text(a, 96, 66, '2 outer events', 9, RED, 'bold')
    region(a, 26, 88, 40, 21)
    graph(a, {1: (13, 88), 2: (26, 88), 3: (39, 88)}, {(1, 2): 1, (2, 3): 3}, labels)
    arrow(a, (47, 88), (54, 88))
    text(a, 57, 88, 'τ = 2: absorb C–N; its change is internal to one fragment', 9,
         TEAL, align='left')
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(man/'figs'/f'fig8_design_examples.{ext}', dpi=240, facecolor='white',
                    bbox_inches='tight', pad_inches=.06)
    plt.close(fig)
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
 def ar(a,b,color=navy,lw=1.9,dashed=False):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=11,lw=lw,color=color,linestyle=(0,(3,3)) if dashed else '-',zorder=2))
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

 # Cut conditions converge into one catalogue. The union is schematic: it does
 # not claim each illustrated cut returns the same fragment or search history.
 box(12,576,708,454,bg,bg,r=16)
 label(34,612,'b','Explore cuts; collect the retained families')
 t(34,646,'Methyl acetate: three base conditions are shown.',10.5,muted)
 for y,title,cut in [(675,'Uncut',None),(780,'Cut O4-C5',(4,5)),(885,'Cut C2-O4',(2,4))]:
  t(40,y+39,title,10.5,green,'bold')
  coords=molecule(ax,ESTER,157,y-7,280,96,{i:PALE_A for i in [1,2,3,4,5]},notes={2:'2',4:'4',5:'5'},font=24)
  if cut:
   start,end=coords[cut[0]],coords[cut[1]]
   left=tuple(start[k]+.23*(end[k]-start[k]) for k in range(2))
   right=tuple(start[k]+.77*(end[k]-start[k]) for k in range(2))
   ax.plot([left[0],right[0]],[left[1],right[1]],color='white',lw=4.5,zorder=7)
   ax.plot([left[0],right[0]],[left[1],right[1]],color=violet,lw=1.5,ls=(0,(2,3)),zorder=8)
  tree([(446,y+39),(491,y+39),(491,819),(564,819)],navy,1.5)
 t(616,800,r'$\bigcup$',28,navy,ha='center')
 t(616,870,'Family\ncatalogue',11,navy,'bold',ha='center')
 t(616,940,'Across seeds\nand cuts',10,muted,ha='center')
 t(34,1006,'Merge repeated families; retain their constraints.',10,muted)

 # Local decoding reads shared fragments and their saved context. The two
 # pictured realizations share two boundary changes and differ internally.
 box(738,576,690,454,bg,bg,r=16)
 label(760,612,'c','Decode shared fragments and boundaries')
 t(760,646,'Use original bond weights and the saved joint choices.',10.5,muted)
 t(867,699,'Unique fragment',11,violet,'bold',ha='center')
 molecule(ax,ACETATE,761,728,210,159,{3:PALE_O,4:PALE_O},notes={3:'3',4:'4'},font=33)
 ar((655,819),(758,819),navy,1.7)
 t(866,907,r'Allowed O exchange: $3\leftrightarrow4$',10,violet,ha='center')
 t(866,958,'One compressed family\nshared local responses',10,violet,ha='center')
 for y,mapping,count,color in [(690,{3:3,4:4,5:5,6:6},2,green),(861,{3:4,4:3,5:5,6:6},4,RED)]:
  box(1056,y,347,135,'white','#D6E1E9',.9,10)
  molecule(ax,ACETATE,1060,y+2,146,104,{1:PALE_A,2:PALE_A,3:PALE_O,4:PALE_O},notes={i:mapping[i] for i in [3,4]},font=29,bond_changes={} if count==2 else {(2,3):1,(2,4):-1})
  t(1209,y+51,'+',11)
  molecule(ax,ALCOHOL,1220,y+27,95,72,{5:PALE_B},notes={6:'6'},font=28,bond_changes={(5,6):1})
  t(1359,y+43,str(count),20,color,'bold',ha='center')
  t(1359,y+74,'events',10,color,ha='center')
  t(1230,y+119,f'2 boundary + {count-2} internal',10,muted,ha='center')
 ar((974,799),(1053,751),violet,1.8,True)
 ar((974,825),(1053,922),violet,1.8,True)
 t(1230,844,'Exchange O3 / O4',9.5,violet,ha='center')

 # Preserve the original chemical/event legend; search cuts are not events.
 ar((35,1064),(84,1064),navy,1.7);t(95,1064,'Search / aggregation',10,muted)
 ar((414,1064),(463,1064),violet,1.7,True);t(474,1064,'Event projection',10,muted)
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
 print('Molecular overview: witness growth, conditional placements, cut aggregation and local events.')

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
