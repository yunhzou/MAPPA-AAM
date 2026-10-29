"""Reference-blind inputs and symmetry-aware scoring of saved AAM graphs.

Coverage queries use Schreier-Sims transversal choices, not group-element or
bijection enumeration. Reference labels are only consumed after search.
"""
import ast
from collections import defaultdict
import time

import numpy as np
import pynauty
from rdkit import Chem

from mappa.mappa_core.domain import AAMProblem, MolecularEndpoint
from mappa.mappa_workflow.mapping_analytics.family import query_path
from mappa.mappa_core.core.search.graph import frozen_value


def prepare(mapped_reaction):
    endpoints, features, labels = [], [], []
    for side, smiles in zip(('R', 'P'), mapped_reaction.split('>>')):
        mol = Chem.MolFromSmiles(smiles)
        maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.MolToSmiles(mol, canonical=True)
        order = ast.literal_eval(mol.GetProp('_smilesAtomOutputOrder'))
        mol = Chem.AddHs(Chem.RenumberAtoms(mol, order))
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
        hindex = {a:i for i,a in enumerate(heavy)}
        label = {maps[old]: new for new,old in enumerate(order) if maps[old]}
        labels.append(label)
        matrix = np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()))
        bonds = []
        for bond in mol.GetBonds():
            a,b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            matrix[a,b] = matrix[b,a] = bond.GetBondTypeAsDouble()
            if a in hindex and b in hindex:
                bonds.append((hindex[a],hindex[b],(bond.GetBondTypeAsDouble(),str(bond.GetStereo()))))
        colors = []
        for i in heavy:
            atom = mol.GetAtomWithIdx(i)
            colors.append((atom.GetSymbol(),atom.GetFormalCharge(),atom.GetIsotope(),
                           atom.GetTotalNumHs(includeNeighbors=True),
                           atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else ''))
        features.append(dict(heavy=heavy, colors=colors, bonds=bonds))
        endpoints.append(MolecularEndpoint(tuple(a.GetSymbol() for a in mol.GetAtoms()),
            np.zeros((mol.GetNumAtoms(),3)), matrix, side,
            metadata={'unmapped_smiles':Chem.MolToSmiles(mol), 'scoring_features':features[-1]}))
    reference = {labels[0][k]:labels[1][k] for k in labels[0].keys() & labels[1].keys()}
    return AAMProblem(*endpoints), features, reference


def colored_graph(features, mapping=None):
    adjacency, colors = defaultdict(set), defaultdict(set)
    offsets, count = [], 0
    for side, feature in enumerate(features):
        offsets.append(count)
        for color in feature['colors']:
            colors[('atom',side,tuple(color))].add(count); count += 1
    def link(a,b,color):
        nonlocal count
        colors[color].add(count)
        adjacency[count].update((a,b)); adjacency[a].add(count);adjacency[b].add(count)
        count += 1
    for side,feature in enumerate(features):
        for a,b,color in feature['bonds']:
            link(offsets[side]+a,offsets[side]+b,('bond',side,tuple(color)))
    if mapping is not None:
        for a,b in sorted(mapping.items()):
            link(a,offsets[1]+b,('mapping',))
    return pynauty.Graph(count, adjacency_dict={i:sorted(adjacency[i]) for i in range(count)},
                         vertex_coloring=[colors[k] for k in sorted(colors,key=repr)])


def project(mapping, features):
    source = {a:i for i,a in enumerate(features[0]['heavy'])}
    target = {a:i for i,a in enumerate(features[1]['heavy'])}
    return {source[a]:target[b] for a,b in dict(mapping).items() if a in source and b in target}


def endpoint_generators(feature):
    n = len(feature['heavy'])
    return tuple(tuple(g[:n]) for g in pynauty.autgrp(colored_graph([feature]))[0])


def exact_action(g, feature):
    colors = feature['colors']
    if any(colors[i] != colors[g[i]] for i in range(len(g))):
        return False
    bonds = {(min(a,b),max(a,b)):tuple(c) for a,b,c in feature['bonds']}
    return all(bonds.get(tuple(sorted((g[a],g[b])))) == c for (a,b),c in bonds.items())


def changes_chemical_correspondence(placement, features):
    """Only choices that change chemical correspondence need a reference query."""
    target_index={a:i for i,a in enumerate(features[1]['heavy'])}
    projected=tuple(dict.fromkeys(tuple(target_index[g.images[a]] for a in features[1]['heavy'])
                                   for g in placement.target_generators))
    projected=tuple(g for g in projected if g!=tuple(range(len(target_index))))
    changes=any(not exact_action(g,features[1]) for g in projected)
    for domain in placement.symmetry_domains:
        if domain.source=='exact_automorph_group':continue
        atoms=[target_index[a] for a in domain.p_atoms if a in target_index]
        for atom in atoms[1:]:
            swap=list(range(len(target_index)));swap[atoms[0]],swap[atom]=swap[atom],swap[atoms[0]]
            changes |= not exact_action(swap,features[1])
    return changes


def reference_event_index(features):
    """Count only annotated-endpoint heavy bonds, using the shared event rule.

    Chemical reference labels define this scope; full decoded all-H counts
    cannot be compared with an annotation lacking hydrogen correspondence.
    """
    from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
    endpoints=[]
    for feature in features:
        elements=tuple(color[0] for color in feature['colors'])
        weights=np.zeros((len(elements),len(elements)))
        for a,b,label in feature['bonds']:weights[a,b]=weights[b,a]=label[0]
        endpoints.append(MolecularEndpoint(elements,np.zeros((len(elements),3)),weights))
    return SignedEventIndex(AAMProblem(*endpoints))


def reference_membership(problem, representatives, paths, features, reference, *, seconds=120):
    """Membership in saved compressed choices, with the original Golden policy.

    ``paths`` is lazy: representative hits need no constraint construction.
    Unmatched heavy atoms remain part of the exact relation. Hydrogens are
    feasible but unscored; endpoint chemical equivalence preserves correlations.
    """
    start, cpu = time.perf_counter(), time.process_time()
    deadline = start + seconds
    reference = dict(reference)
    expected = pynauty.certificate(colored_graph(features, project(reference, features)))
    event_index=reference_event_index(features)
    expected_events=int(event_index.counts([event_index.vector(project(reference,features))])[0].sum())
    certificates = {}
    queries = []
    filtered=0
    def matches(mapping):
        nonlocal filtered
        heavy = tuple(sorted(project(mapping, features).items()))
        if heavy not in certificates:
            count=int(event_index.counts([event_index.vector(dict(heavy))])[0].sum())
            if count!=expected_events:
                filtered+=1
                certificates[heavy]=None
            else:
                certificates[heavy] = pynauty.certificate(colored_graph(features, dict(heavy)))
        return certificates[heavy] == expected
    def answer(status, witness=None, **details):
        return dict(status=status,witness=witness,queries=queries,**details,
            event_filter=dict(scope='annotated_endpoint_heavy_bonds',reference_count=expected_events,
                              rejected_representatives=filtered),
            timing=dict(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu))
    for identifier, mapping in representatives:
        if matches(mapping):
            return answer('recovered',dict(mapping=sorted(dict(mapping).items()),branch=identifier,
                                           method='representative'))
        if time.perf_counter() >= deadline:
            return answer('unknown',reason='reference-check budget')
    def lifted(feature, degree):
        result=[]
        for generator in endpoint_generators(feature):
            images=list(range(degree))
            for i, atom in enumerate(feature['heavy']):
                images[atom]=feature['heavy'][generator[i]]
            result.append(tuple(images))
        return tuple(result)
    source = lifted(features[0], problem.source_atom_count)
    target = lifted(features[1], problem.target_atom_count)
    seen=set(); unknown=False
    for identifier, path in paths():
        remaining=deadline-time.perf_counter()
        if remaining<=0:
            return answer('unknown',reason='reference-check budget')
        if len(project(path.mapping,features)) != len(project(reference,features)):
            continue
        if not any(changes_chemical_correspondence(placement,features)
                   for edge in path.transitions
                   if (placement:=path.graph.fragment_placement(edge)) is not None):
            continue
        identity=(tuple(sorted(path.mapping.items())),path.context.cuts,
                  path.context.graph_floor,path.context.iso_tolerance,frozen_value(path.fragments))
        if identity in seen:
            continue
        seen.add(identity)
        status, result=query_path(path,problem,reference,source_atoms=features[0]['heavy'],
            source_generators=source,target_generators=target,complete_reference=True,
            timeout_ms=min(10000,1000*remaining))
        queries.append(dict(branch=identifier,status=status,**result))
        if status=='recovered':
            if not matches(result['mapping']):
                raise ValueError('Reference witness fails independent chemical equivalence check')
            return answer(status,dict(result,branch=identifier,method='compressed_choices'))
        unknown |= status=='unknown'
    return answer('unknown' if unknown else 'not_recovered')


def rank_key(mapping, problem):
    """Frozen graph-only representative rank; no reference data consumed."""
    mapping = dict(mapping)
    inverse = {p:r for r,p in mapping.items()}
    r,p = problem.reactant.wbo,problem.product.wbo
    broken=formed=changed=0
    for a,b in zip(*np.where(np.triu(r,1)>.2)):
        if a not in mapping or b not in mapping:
            broken += int(a in mapping or b in mapping)
        elif p[mapping[a],mapping[b]] <= .2: broken += 1
        elif abs(r[a,b]-p[mapping[a],mapping[b]])>.5: changed += 1
    for a,b in zip(*np.where(np.triu(p,1)>.2)):
        if a not in inverse or b not in inverse or r[inverse.get(a,0),inverse.get(b,0)]<=.2:
            formed += 1
    heavy = sum(problem.product.elements[i]!='H' for i in mapping.values())
    return (-heavy,-len(mapping),broken+formed+changed,tuple(sorted(mapping.items()))),dict(
        broken=int(broken),formed=int(formed),bond_order_changed=int(changed))


def evaluate(aam, features, reference, seconds=120, reference_side='target', ranker=None):
    start=time.perf_counter(); deadline=start+seconds
    full_reference=dict(reference)
    reference=project(reference,features)
    side={'source':0,'target':1}[reference_side]
    complete_reference=len(reference)==len(features[side]['heavy'])
    ranker = ranker or (lambda mapping:rank_key(mapping,aam.problem))
    expected=pynauty.certificate(colored_graph(features,reference))
    terminals=list(aam.graph.terminals)
    ranked=sorted(terminals,key=lambda t:ranker(aam.graph.states[t].mapping)[0])
    candidates=[project(aam.graph.states[t].mapping,features) for t in ranked]
    # Explicit-H alternatives often have the identical heavy-atom relation.
    # Canonicalize that relation once, without discarding any AAM terminal.
    certificate_cache={}
    certificates=[]
    for mapping in candidates:
        key=tuple(sorted(mapping.items()))
        if key not in certificate_cache:
            certificate_cache[key]=pynauty.certificate(colored_graph(features,mapping))
        certificates.append(certificate_cache[key])
    hit=next((i for i,c in enumerate(certificates) if c==expected),None)
    result=dict(reference_pairs=len(reference),reference_annotation_complete=complete_reference,
        top1_correct=bool(certificates and certificates[0]==expected),
        representative_recovery=hit is not None,reference_recovery='recovered' if hit is not None else 'not_recovered',
        top_terminal=ranked[0] if ranked else None,
        witness_terminal=ranked[hit] if hit is not None else None,
        candidate_terminals=len(terminals),unique_representative_chemistries=len(set(certificates)),
        unique_heavy_representatives=len(certificate_cache),
        best_target_heavy_coverage=max((len(m)/max(1,len(features[side]['heavy'])) for m in candidates),default=0),
        best_target_all_atom_coverage=max((len(aam.graph.states[t].mapping)/(aam.problem.source_atom_count if side==0 else aam.problem.target_atom_count) for t in terminals),default=0),
        capped=aam.graph.capped, compatibility_queries=0, unknown_queries=0)
    if hit is None:
        answer=reference_membership(aam.problem,(),
            lambda:(((path.terminal,tuple(path.transitions)),path) for path in aam.graph.paths()),
            features,full_reference,seconds=max(0,deadline-time.perf_counter()))
        result.update(reference_recovery=answer['status'],
                      compatibility_queries=len(answer['queries']),
                      unknown_queries=sum(q['status']=='unknown' for q in answer['queries']),
                      reference_check=answer)
        if answer['witness'] is not None:
            witness=answer['witness']
            terminal,transitions=witness['branch']
            witness['heavy_mapping']=sorted(project(witness['mapping'],features).items())
            result.update(witness_terminal=terminal,witness_path=list(transitions),witness_actions=witness)
    result['evaluation_seconds']=time.perf_counter()-start
    result['verifier']='ground_truth_fragment_constraints_v1'
    result['reference_semantics']='exact heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry'
    if ranked:result['top_events']=ranker(aam.graph.states[ranked[0]].mapping)[1]
    return result


def evaluate_planned(aam, plan, features, reference, **options):
    """Score native compressed actions in search space; rank in input R/P space."""
    ordered_features = list(reversed(features)) if plan.reversed else features
    result = evaluate(aam,ordered_features,plan.to_search_mapping(reference),
        reference_side='source' if plan.reversed else 'target',
        ranker=lambda m:rank_key(plan.to_input_mapping(m),plan.input_problem),**options)
    result['search_direction']=plan.direction
    result['mapping_space']='search source/target; use orientation.json to recover input R/P indices'
    result['coverage_space']='original input product'
    if result['witness_terminal'] is not None and result['representative_recovery']:
        result['input_orientation_witness']=sorted(plan.to_input_mapping(
            aam.graph.states[result['witness_terminal']].mapping).items())
    elif result.get('witness_actions'):
        heavy=dict(result['witness_actions']['heavy_mapping'])
        concrete={ordered_features[0]['heavy'][r]:ordered_features[1]['heavy'][p]
                  for r,p in heavy.items()}
        result['input_orientation_heavy_witness']=sorted(plan.to_input_mapping(concrete).items())
        result['input_orientation_witness']=sorted(plan.to_input_mapping(result['witness_actions']['mapping']).items())
    return result
