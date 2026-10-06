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
  return [[i,j,1 if d>0 else -1] for i in zr for j in zr if i<j and abs(d:=wp.get(tuple(sorted((m[i],m[j]))),0)-wr.get((i,j),0))>=.5]
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

def bond_event_marker(ax, a, b, sign, size=None):
 """Small, vivid marks centered on a black bond; no filled badge or leader.

 Red crossing strokes indicate loss; two green inward arrows indicate gain.
 The bond-axis coordinates keep the glyph attached to slanted bonds as well.
 """
 import math
 from matplotlib import patheffects
 length = math.dist(a, b)
 if not length: return
 ux, uy = (b[0]-a[0])/length, (b[1]-a[1])/length
 mid = ((a[0]+b[0])/2, (a[1]+b[1])/2)
 size = size if size is not None else min(8.5, max(5.0, length*.15))
 def point(x,y): return (mid[0]+ux*x-uy*y, mid[1]+uy*x+ux*y)
 def stroke(points):
  q=[point(x,y) for x,y in points]
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

def molecule(ax,smiles,x,y,w,h,owners=None,notes=None,font=24,explicit_hydrogens=False,align_acetyl=True,rotation=0,bond_changes=None,preserve_bond_orders=False,coordinates=None):
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
 # A connected branching overview, with chemical drawings as the tree states.
 fig=plt.figure(figsize=(10.8,8.15));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1440),ylim=(1087,0));ax.axis('off')
 navy='#22324D';muted='#64748B';faint='#BCC7D5';bg='#F5F7FB';violet='#7957BA';orange='#CC762F';green='#128477'
 def t(x,y,s,size=11,color=navy,weight='normal',ha='left'):
  return ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.3,zorder=6)
 def box(x,y,w,h,fc='white',ec=LINE,lw=.8,r=12):
  patch=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',facecolor=fc,edgecolor=ec,lw=lw,zorder=1)
  ax.add_patch(patch)
 def ar(a,b,color=navy,lw=1.9,dashed=False,rad=0,style='-|>'):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=11,lw=lw,color=color,linestyle=(0,(3,3)) if dashed else '-',connectionstyle=f'arc3,rad={rad}',zorder=2))
 def tree(points,color=navy,lw=2.0):
  xs,ys=zip(*points[:-1]);ax.plot(xs,ys,color=color,lw=lw,solid_joinstyle='round',solid_capstyle='round',zorder=2)
  ar(points[-2],points[-1],color,lw)
 def dot(x,y,color=navy,r=4.3):ax.add_patch(Circle((x,y),r,fc=color,ec='none',zorder=3))
 def label(x,y,letter,title):
  t(x,y,letter,14,navy,'bold');t(x+36,y,title,14,navy,'bold')
 def continuation(x,y):
  tree([(x,y),(x+35,y),(x+35,y-27),(x+96,y-27)],faint,1.3)
  tree([(x,y),(x+35,y),(x+35,y+27),(x+96,y+27)],faint,1.3)
  dot(x+35,y,faint,3)
  for yy in [y-27,y+27]:t(x+114,yy,'···',13,faint)

 # A: two genuine conditional forks, followed by the forced third placement.
 box(12,12,1416,557,bg,bg,r=16)
 label(34,47,'a','MAPPA growth: earlier placements condition every later choice')
 t(34,80,'Local alcohol fragments A, B and C share this motif; R lies outside the matched atoms.',10.5,muted)
 for x,atoms,name in [(24,[11],'Seed'),(163,[11,12],'Grow'),(302,[11,12,13],'Match O-H')]:
  t(x+64,119,name,10,green,'bold',ha='center')
  molecule(ax,ALCOHOL_FRAGMENT,x,137,130,116,{i:PALE_A for i in atoms},notes=False,font=28)
 ar((149,195),(160,195),green,1.4);ar((288,195),(299,195),green,1.4)
 t(236,280,'Three distinct target sites',10.5,navy,'bold',ha='center')
 molecule(ax,ALCOHOL_SITES,33,302,405,207,{i:PALE_A for i in [21,22,23,31,32,33,41,42,43]},notes={23:'a',33:'b',43:'c'},font=26)
 t(239,536,'Each site preserves C-C, C-O and O-H.',10,muted,ha='center')
 ax.plot([458,458],[108,545],color='#DDE4EE',lw=.8)
 t(680,120,'Place A',11,green,'bold',ha='center')
 t(956,120,'Place B, given A',11,orange,'bold',ha='center')
 t(1258,120,'Place C, given A and B',10.5,violet,'bold',ha='center')
 root=(485,330)
 def card(x,y,w,main,sub,color):
  box(x,y-28,w,56,'white','#D6E1E9',.8,8)
  t(x+w/2,y-8,main,12,color,'bold',ha='center')
  t(x+w/2,y+14,sub,9,muted,ha='center')
 for site,y,remaining in [('a',206,['b','c']),('b',338,['a','c']),('c',470,['a','b'])]:
  tree([root,(525,330),(525,y),(588,y)],navy,1.7)
  card(592,y,176,'A → '+site,'lock '+site,green)
  for dest,yy in zip(remaining,[y-36,y+36]):
   last=next(k for k in remaining if k!=dest)
   tree([(770,y),(807,y),(807,yy),(858,yy)],navy,1.6)
   card(862,yy,187,'B → '+dest,'A:'+site+' stays fixed',orange)
   ar((1052,yy),(1133,yy),navy,1.5)
   card(1137,yy,249,'C → '+last,'only unused site remains',violet)
  dot(807,y,navy,3)
 dot(525,330,navy,3)
 t(953,550,'Each path is a joint assignment; occupied target atoms cannot be reused.',10,muted,ha='center')

 # Lower panels use the same geometry, translated down for the deeper tree.
 original_t,original_box,original_ar,original_tree,original_dot,original_molecule=t,box,ar,tree,dot,molecule
 def t(x,y,*args,**kw):return original_t(x,y+120,*args,**kw)
 def box(x,y,*args,**kw):return original_box(x,y+120,*args,**kw)
 def ar(a,b,*args,**kw):return original_ar((a[0],a[1]+120),(b[0],b[1]+120),*args,**kw)
 def tree(points,*args,**kw):
  # original_tree resolves the wrapped ar; translate its line directly here.
  color=args[0] if args else navy;lw=args[1] if len(args)>1 else 2.0
  xs,ys=zip(*points[:-1]);ax.plot(xs,[y+120 for y in ys],color=color,lw=lw,zorder=2)
  ar(points[-2],points[-1],color,lw)
 def dot(x,y,*args,**kw):return original_dot(x,y+120,*args,**kw)
 def mol(smiles,x,y,*args,**kw):return original_molecule(ax,smiles,x,y+120,*args,**kw)
 # B: show source-edge conditions, not takeover or a fabricated trajectory.
 box(12,467,702,448,bg,bg,r=16)
 label(34,503,'b','Vary the constraints with a cut sweep')
 t(34,536,'Methyl acetate: three illustrated source conditions',10.5,muted)
 for y,title,cut in [(565,'Uncut source',None),(668,'Mask O4-C5',(4,5)),(771,'Mask C2-O4',(2,4))]:
  t(46,y+39,title,10.5,green,'bold')
  coords=mol(ESTER,214,y-7,306,96,{i:PALE_A for i in [1,2,3,4,5]},notes={2:'2',4:'4',5:'5'},font=24)
  if cut:
   a0,b0=coords[cut[0]],coords[cut[1]]
   a0,b0=tuple(a0[k]+.23*(b0[k]-a0[k]) for k in range(2)),tuple(a0[k]+.77*(b0[k]-a0[k]) for k in range(2))
   ax.plot([a0[0],b0[0]],[a0[1],b0[1]],color='white',lw=4.5,zorder=7)
   ax.plot([a0[0],b0[0]],[a0[1],b0[1]],color=violet,lw=1.5,ls=(0,(2,3)),zorder=8)
  ar((535,y+39),(578,y+39),navy,1.5)
  t(590,y+39,'Grow',10.5,navy)
 t(370,887,'Pool the retained branches; decode on original bond weights.',9.8,muted,ha='center')

 # C: output projection uses dashed arrows; it is not a search branch.
 box(732,467,696,448,bg,bg,r=16)
 label(754,503,'c','Decode each compressed family')
 t(754,536,'The oxygen exchange preserves matching at '+r'$\tau_{\mathrm{iso}}=1$',10.5,muted)
 box(757,621,216,189,'white','#CEC1E3',1.0,12)
 coords=mol(ACETATE,759,630,210,159,{3:PALE_O,4:PALE_O},notes={3:'c',4:'d'},font=33)
 t(865,594,'Matching symmetry',10.5,violet,'bold',ha='center')
 t(865,790,r'$g=(c\ d)$',12,violet,ha='center')
 t(865,851,'One branch\ncorrelated alternatives',10.5,violet,ha='center')
 # Witness cards show oxygen identities while leaving unchanged carbons unlabeled.
 for y,mapping,count,color in [(575,{3:3,4:4,5:5,6:6},2,green),(757,{3:4,4:3,5:5,6:6},4,RED)]:
  box(1056,y,347,144,'white','#D6E1E9',.9,10)
  mol(ACETATE,1060,y+10,146,118,{1:PALE_A,2:PALE_A,3:PALE_O,4:PALE_O},notes={i:mapping[i] for i in [3,4]},font=29,bond_changes={} if count==2 else {(2,3):1,(2,4):-1})
  t(1209,y+60,'+',11)
  mol(ALCOHOL,1220,y+39,95,72,{5:PALE_B},notes={6:'6'},font=28,bond_changes={(5,6):1})
  t(1359,y+51,str(count),20,color,'bold',ha='center')
  t(1359,y+83,'events',10,color,ha='center')
 ar((975,704),(1053,647),violet,1.8,True)
 ar((975,730),(1053,829),violet,1.8,True)
 t(1230,739,'Exchange O3 / O4',10,violet,ha='center')
 # Shared legend, kept outside the chemical trees.
 ar((35,949),(84,949),navy,1.7);t(95,949,'Search continuation',10,muted)
 ar((394,949),(443,949),violet,1.7,True);t(454,949,'Event projection',10,muted)
 for x,sign,label in [(886,-1,'Break / weaken'),(1155,1,'Form / strengthen')]:
  ax.plot([x,x+50],[1069,1069],color='black',lw=1.1)
  bond_event_marker(ax,(x,1069),(x+50,1069),sign,size=8)
  t(x+61,949,label,10,muted)
 fig.canvas.draw()
 renderer=fig.canvas.get_renderer()
 for label in ax.texts:
  extent=label.get_window_extent(renderer)
  assert extent.x0>=-1 and extent.x1<=fig.bbox.width+1 and extent.y0>=-1 and extent.y1<=fig.bbox.height+1,label.get_text()
 for ext in ['pdf','svg','png']:fig.savefig(man/'figs'/f'fig1_algorithm.{ext}',facecolor='white',dpi=260)
 plt.close(fig)
 print('Branching overview: RDKit vector molecules and checked illustrative mappings.')

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
