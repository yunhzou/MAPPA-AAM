# MAPPA search, signed-event decoding, and explicit downstream consumers

The main computation is `search_aam -> decode_events`. Search retains fragment
choices and correlated symmetry actions. Event decoding finds one witness per
signed bond-event class within the saved families and requested event window.
It does not enumerate all atom bijections or select a geometry-ranked mapping.

```text
AAMProblem + AAMSearchConfig
          |
      search_aam
          |
       AAMResult -----------------> raw graph / growth viewer
          |
     decode_events
          |
     DecodedEvents ---------------> family-constrained symmetry queries
          |
          +-- explicitly chosen mapping/events: ReactionContext -> TS analysis
```

`group_mechanisms`, `compile_mapping_families`, and
`compile_mechanism_families` are separate family-analysis utilities. They do not
implicitly select a 3D orientation. Their exact relational code lives in
`alignment/mapping_families.py`.

The legacy R/P orientation-graph selector and its convenience wrappers have
been removed. There is no compatibility route that invokes them. The CLI
uses `events` by default and accepts only `aam`, `events`, or explicit `ts`.

## 1. Fragment Candidate State

### 1.1 `_SymBlock`

A symmetry block represents a set-valued assignment domain:

```text
SymBlock
|- r_atoms: R roles participating in this domain
|- p_atoms: interchangeable target pool
`- extendable: whether later growth may add another R role
```

It does not mean that every factorial permutation is valid. Correlations are
resolved by the exact candidate automorphism group.

### 1.2 `_SymCand`

One compressed candidate stores:

```text
_SymCand
|- mapping: one deterministic representative R -> P assignment
|- blocks: open/closed assignment domains
|- exact_fixed: roles individualized by prior constraints
|- multiplicity: number of encountered states represented
|- automorph_blocks: connected display domains
`- automorph_generators: authoritative exact target permutations
```

The representative is needed to continue deterministic code, but it is not
treated as the only mapping represented by the candidate.

### 1.3 Candidate canonicalization

After extension, `_CandidateAutomorphismCanonicalizer` constructs a colored
pynauty graph containing:

- target atom elements and WBO buckets;
- fixed target atoms;
- distinct R-role colors;
- set-valued symmetry pools;
- already locked assignments.

Two candidates are merged only when equal canonical certificates and an exact
transporter prove that one full state maps to the other. The transporter is
retained as an exact automorphism generator. Merely placing two atoms in the
same vertex orbit is not sufficient because independent orbit labels discard
correlations.

## 2. Weighted Island Growth

`grow_island` starts from one R seed and uses a priority queue of R edges,
normally strongest WBO first. Queue order affects traversal, not validity.

For a proposed new R atom `n` and candidate target `v`, all active edges from
`n` to the already grown R fragment are checked:

```text
element_R[n] == element_P[v]
v is unused

for every grown r with WBO_R[n,r] >= graph_floor:
    WBO_P[v,m(r)] >= graph_floor
    abs(WBO_R[n,r] - WBO_P[v,m(r)]) <= iso_tol
```

R non-edges are not local negative constraints. An extra P bond may be a
formed bond and is classified at mechanism level.

If a touched R role belongs to a symmetry block, the support question is
existential and correlated:

```text
Does one injective assignment inside all touched blocks make the complete
active-edge WBO vector valid?
```

This small internal matching is bounded by `SYM_SUPPORT_MAX_STATES`. It is not
a global witness enumeration.

The growth transition is explicitly one of:

```text
0 outputs       -> defer the observed boundary; keep the candidate alive
1 output        -> commit the unique compressed state
many outputs    -> merge automorphic outputs; branch only on distinct states
```

At saturation, the fragment record retains:

- its R atom set;
- its representative assignments;
- exact target automorphism generators;
- symmetry domains;
- deferred boundary edges.

## 3. Deferred Boundaries

A fragment that cannot absorb a frontier atom is not necessarily invalid.
The failed weighted relationship is stored as a deferred boundary. This is
essential when an internally symmetric fragment has two sides but only one
side has already encountered another island.

Candidate equality therefore includes:

```text
exact internal pynauty certificate
              +
deferred one-hop boundary state
```

Without the second term, future-distinguishable branches would be collapsed.

## 4. Multi-Fragment AAM Branches

`find_islands` grows fragments sequentially for each seed order.

### 4.1 Live branch object

```text
_Branch
|- mapping: current concrete representative
|- islands_R: R atom -> island ID
|- islands_P: P atom -> island ID
|- deferred_edges
`- symmetry_paths[]
    `- ordered committed fragment records
```

An AAM branch is therefore a combination of fragment candidates, not a single
fragment and not merely a final mapping.

Hard `anchor_map` pairs are preloaded into the branch as fixed mapping/island
state. An anchored atom may seed growth, but it can also remain outside every
grown fragment. During analytical compilation, such an uncovered anchor is
represented as an individually fixed singleton fragment. Missing non-anchor
atoms remain an error. Anchor colors are also carried into the relational
graph, so family dedupe, chirality, symmetry repair, and RMSD selection cannot
move an anchored pair.

### 4.2 Exact live dedupe

Live branches merge only when these concrete cumulative states are equal:

```text
(mapping, islands_R, islands_P, deferred_edges)
```

Their distinct fragment histories are retained in `symmetry_paths`.

We tested a more aggressive coupled R/P automorphic live-state quotient. It
reduced TS04 from four correct mechanisms to two because the live graph did
not encode the full accumulated hierarchy. That optimization was removed.
Semantic automorphic dedupe is delayed until the completed hierarchy exists.

### 4.3 Branch cap

The configured BGCP cap is normally `max_branches = 100`. It applies to the
post-dedupe live leaves of one parent subtree:

```text
if accepting this parent's descendants would make live leaves > 100:
    discard only this overflowing descendant subtree
    retain siblings and other seed paths
```

Exactly 100 is legal. This is an intentional computational completeness cap,
not a symmetry equivalence rule.

## 5. Completed AAM Data Hierarchy

The raw pool and typed post-AAM model have the following ownership:

```text
CutSweepPool
`- mechanism entry, keyed by exact event certificate
   |- representative_mapping
   |- cuts / has_no_cut / encounter count
   `- completed branches[]
      `- AAMBranch
         |- representative_mapping
         |- encounter_count
         |- cuts
         |- hierarchy: AAMHierarchy
         |  `- fragments[]: FragmentMatch
         |     |- R atom set
         |     |- island ID
         |     |- deferred edges
         |     |- symmetry domains
         |     `- exact target generators
         |- analytical mapping-family record
         `- path provenance[]

PostAAMMechanism
|- exact mechanism key
|- endpoint R automorphism group   [auxiliary graph information]
|- endpoint P automorphism group   [auxiliary graph information]
`- maximal AAMBranch families[]    [authoritative mapping choices]
```

Endpoint automorphism groups are not automatically free mapping candidates.
They describe endpoint graph symmetry. Allowed mapping changes come from the
selected branch's analytical relation/coset.

### 5.1 Meaning of a witness

A witness is one complete bijection inside a family. It is useful for:

- continuing deterministic computation;
- recording provenance;
- initializing an exact isomorphism.

It is not sampled uniformly, is not a generator, and must not be ranked as if
the encountered witness list were the solution space.

## 6. Exact Mechanism Identity

For a complete mapping, bond events are computed in R index order:

```text
broken if WBO_R[i,j] - WBO_P[m(i),m(j)] >= event_threshold(i,j)
formed if WBO_P[m(i),m(j)] - WBO_R[i,j] >= event_threshold(i,j)
```

The entire broken/formed edge set is then attached to the full WBO-colored R
graph using typed event vertices. A pynauty certificate of this decorated
graph is the mechanism key.

```text
R graph atom vertices
      |
      +-- WBO-colored graph edges
      |
      +-- broken-event vertices -- broken type marker
      |
      `-- formed-event vertices -- formed type marker
                         |
                         v
               pynauty canonical certificate
```

This fixes a subtle but important bug: equal endpoint vertex-orbit IDs do not
prove that two edges lie in the same edge orbit, and they certainly do not
prove that two multi-edge event sets are equivalent.

Mechanism selection keeps the classes with the minimum number of bond-breaking
plus bond-forming events. Different exact event certificates at that minimum
remain distinct mechanisms.

## 7. Analytical Mapping Families

Every completed branch is compiled into a colored relational isomorphism
between endpoint A and endpoint B.

Two nested relations are compiled from the same completed AAM branch:

- atom element and fragment-owner colors;
- selected fragment WBO relations at `iso_tol`;
- anchors;
- the structural relation stops here and represents the pre-event AAM group;
- the conservative event stabilizer additionally contains event-invariant
  pair colors;
- coordinate orientation is not part of analytical family compilation.

If `g` is one isomorphism and `G` is the target automorphism group of the
relation, the branch represents the coset:

```text
F = gG
```

`AnalyticalMappingFamily` stores:

```text
AnalyticalMappingFamily
|- source_mapping
|- representative_mapping g
|- target_generators of G
|- target_orbits
|- group_order
|- structural_target_generators
|- structural_group_order
`- colored relational records for exact membership
```

The event-invariant colors preserve threshold behavior against every endpoint
WBO value. That makes them a valid subgroup stabilizer, but occasionally
stricter than preserving the concrete event realized by the selected mapping.
The algorithm therefore traverses the finite quotient of the structural AAM
group by this conservative stabilizer. It retains one representative for each
quotient coset whose actual broken/formed event equals the mechanism event.
These event-family representatives remain available for explicit family inspection.

This is a Schreier-style traversal of event cosets, not enumeration of group
elements or atom bijections. TS01, for example, has two same-event quotient
cosets even though the conservative stabilizer contains only one of them.

### 7.1 Membership

For a proposed complete mapping `m`, `contains(m)` directly transports all A
atom colors and relation records through `m` and compares them with B. It does
not invoke a new geometry matching operation.

### 7.2 Family inclusion and dedupe

For finite cosets, `F1` is proven to be a subset of `F2` by checking:

1. `F1`'s representative belongs to `F2`;
2. applying every generator of `F1` to that representative remains in `F2`.

Identical mapping/fragment payloads are quotient-collapsed before family
compilation. Compiled families are then reduced in two exact phases:

1. canonical relation-certificate buckets prove and merge equal cosets;
2. families are visited in descending group order, using orbit compatibility
   to reject impossible containment before relation transport.

The result is a maximal antichain:

```text
equal family       -> merge provenance
strict subset      -> attach its provenance to the containing family
strict superset    -> replace contained families and inherit provenance
incomparable       -> retain both
```

Thus 82 TS01 growth paths become one maximal analytical family without
pretending the paths themselves are group elements. Large cases may retain
several incomparable families under the same mechanism.


## 8. Resource and completeness scope

`branch_limit` bounds live branch admissions within an individual seed/cut
search. It does not bound the accumulated graph history or process RAM.
`search_aam` retains the combined search graph; combining cuts can allocate new
state/transition wrappers while the individual graphs are still live.
`search_aam_checkpoints` persists cuts individually and avoids this combined
history, but an individual cut can still require substantial memory.

Event decoding has separate per-family time and pattern budgets. Its
`complete` flag covers retained families within the requested event window;
it is not an exhaustive-search claim after branch pruning. Resource limits,
incomplete decoding, and solver `unknown` outcomes are not negative proofs.

## 9. Coordinate orientation and transition-state analysis

`analyze_transition_state` consumes an explicitly constructed `ReactionContext`
and a `TransitionStateTarget`. It runs R-to-TS and P-to-TS partial core searches
and scores the resulting assignments. It does not choose or repair the supplied
R-to-P correspondence through another alignment pipeline.

## 10. Validation and viewer boundaries

Witnesses must satisfy their saved family and reported signed bond events.
Independent atom orbits are not interchangeable with correlated group actions.
The offline viewers display saved data and do not rerun matching, change atom
indices, or certify chemistry from an interpolation. Tests cover structural
family membership, current event decoding, coordinate measurements, explicit TS
contexts, and rejection of removed imports and CLI modes.
