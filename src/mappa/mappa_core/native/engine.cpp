// Native compressed matching and fragment growth. Python matching provides an
// independent scientific oracle. Element matching uses this backend by default;
// custom node policies and diagnostic traces use the Python implementation.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <map>
#include <list>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <chrono>

extern "C" {
#include "nauty.h"
}

namespace py = pybind11;

namespace {

using Pair = std::pair<int, int>;

thread_local std::vector<std::vector<int>>* generator_output = nullptr;

void collect_generator(int, int* permutation, int*, int, int, int size) {
    generator_output->emplace_back(permutation, permutation + size);
}

// ---------------------------------------------------------------------------
// nauty
// ---------------------------------------------------------------------------

struct DenseGraph {
    int n = 0;
    int m = 0;
    std::vector<setword> adjacency;
    std::vector<int> lab, ptn, orbits;
    std::vector<setword> workspace;
    optionblk options;
    statsblk stats;

    DenseGraph(int vertices, const std::vector<Pair>& edges)
        : n(vertices), m(SETWORDSNEEDED(vertices)) {
        nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
        adjacency.assign((size_t)n * m, 0);
        for (const auto& e : edges) ADDONEEDGE(adjacency.data(), e.first, e.second, m);
        lab.assign(n, 0);
        ptn.assign(n, 0);
        orbits.assign(n, 0);
        workspace.assign((size_t)66 * m, 0);
        static DEFAULTOPTIONS_GRAPH(defaults);
        options = defaults;
        options.getcanon = FALSE;
        options.defaultptn = FALSE;
        options.digraph = FALSE;
        options.writeautoms = FALSE;
    }

    long generator_calls = 0;

    void partition(const std::vector<std::vector<int>>& cells) {
        int k = 0;
        for (const auto& cell : cells) {
            for (int v : cell) {
                lab[k] = v;
                ptn[k] = 1;
                ++k;
            }
            if (!cell.empty()) ptn[k - 1] = 0;
        }
        if (k != n) throw std::runtime_error("colouring does not cover every vertex");
    }

    std::vector<std::vector<int>> generators(const std::vector<std::vector<int>>& cells) {
        partition(cells);
        ++generator_calls;
        auto settings = options;
        settings.getcanon = FALSE;
        settings.userautomproc = collect_generator;
        std::vector<std::vector<int>> result;
        generator_output = &result;
        nauty(adjacency.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &settings,
              &stats, workspace.data(), (int)workspace.size(), m, n, nullptr);
        generator_output = nullptr;
        return result;
    }
};

// ---------------------------------------------------------------------------
// graphs
// ---------------------------------------------------------------------------

struct Graph {
    int n = 0;
    std::vector<std::string> elem;
    std::vector<int> ecode;             // element code (shared table)
    std::vector<std::vector<int>> adj;  // ascending neighbour lists (active edges)
    std::vector<uint8_t> edge;          // n*n adjacency flags
    std::vector<double> w;              // n*n full WBO matrix
    double bond_cut = 0.2;
    std::vector<double> bond_tolerances; // source pair rules, resolved before growth
    bool tolerance_inclusive = true;

    double tolerance(int a, int b, double uniform) const {
        return bond_tolerances.empty() ? uniform : bond_tolerances[(size_t)a * n + b];
    }

    bool has_edge(int a, int b) const { return edge[(size_t)a * n + b] != 0; }
    // primitives._edge_wbo: 0.0 on the diagonal, else the matrix value
    double wbo(int a, int b) const { return a == b ? 0.0 : w[(size_t)a * n + b]; }
};


// primitives._growth_edge_supported
inline bool growth_edge_supported(double w_R, double w_P, double iso_tol, double floor, bool inclusive = true) {
    double difference = std::fabs(w_R - w_P);
    return w_P >= floor && (inclusive ? difference <= iso_tol : difference < iso_tol);
}

struct Target : Graph {
    std::vector<int> component;    // distinct target components stay separate
    std::vector<int> orbit;        // p_orbits[v]
    // Nauty base graph: atoms, then one vertex per nonzero or active pair.
    int n_vertices = 0;
    std::map<std::pair<int64_t, bool>, std::vector<int>> edge_cells;          // exact WBO colour -> edge vertices
    std::unique_ptr<DenseGraph> nauty;
    std::vector<std::vector<int>> same_element;          // ecode -> atoms ascending
    mutable bool group_prepared = false;
    mutable std::vector<std::vector<int>> global_group;

    const std::vector<std::vector<int>>& global_generators() const {
        if (!group_prepared) {
            std::map<std::pair<std::string,int>,std::vector<int>> atoms;
            for (int p=0;p<n;++p) atoms[{elem[p],component[p]}].push_back(p);
            std::vector<std::vector<int>> cells;
            for (const auto& item:atoms) cells.push_back(item.second);
            for (const auto& item:edge_cells) cells.push_back(item.second);
            global_group=nauty->generators(cells);
            for (auto& action:global_group) action.resize(n);
            group_prepared=true;
        }
        return global_group;
    }
};

// Source topology dependencies for completed-fragment cache validation.
// WBO values/elements never change between its cuts. Target canonicalization
// depends on P, not on reactant topology or reactant automorphism orbits.
struct SourceReads {
    std::vector<uint8_t> rows, pairs;
    explicit SourceReads(int n) : rows(n, 0), pairs((size_t)n*n, 0) {}
};

struct Source : Graph {
    mutable SourceReads* reads = nullptr;
    bool has_edge(int a, int b) const {
        if (reads) {
            uint8_t& first = reads->pairs[(size_t)std::min(a,b)*n+std::max(a,b)];
            first = 1;
        }
        return Graph::has_edge(a,b);
    }
    const std::vector<int>& neighbors(int a) const {
        if (reads) reads->rows[a] = 1;
        return adj[a];
    }
};

struct SourceReadScope {
    const Source& source;
    SourceReads* previous;
    SourceReadScope(const Source& s, SourceReads* reads) : source(s), previous(s.reads) {
        source.reads = reads;
    }
    ~SourceReadScope() { source.reads = previous; }
};

// ---------------------------------------------------------------------------
// candidate state (mappa.matcher.state._SymCand / _SymBlock)
// ---------------------------------------------------------------------------

struct Block {
    std::vector<int> r;  // sorted unique
    std::vector<int> p;  // sorted unique
    bool extendable = true;
    bool open() const { return r.size() < p.size(); }
    bool complete() const { return r.size() >= p.size(); }
};

inline void sort_unique(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

inline Block make_block(std::vector<int> r, std::vector<int> p, bool extendable) {
    sort_unique(r);
    sort_unique(p);
    return Block{std::move(r), std::move(p), extendable};
}

inline bool contains(const std::vector<int>& sorted, int x) {
    return std::binary_search(sorted.begin(), sorted.end(), x);
}


struct Cand {
    std::vector<int> img;     // NR entries, -1 = unmapped
    std::vector<int> mapped;  // ascending mapped R atoms
    std::vector<Block> blocks;
    std::vector<int> exact_fixed;  // sorted unique

    bool has(int r) const { return img[r] >= 0; }
    bool has_open_choice() const {
        for (const auto& b : blocks)
            if (b.open()) return true;
        return false;
    }
};

bool pool_witness(const std::vector<int>& raw, const std::vector<Block>& blocks,
                  int NP, std::vector<int>& img) {
    const int NR = (int)raw.size();
    std::vector<uint8_t> block_r(NR, 0);
    for (const auto& b : blocks)
        for (int r : b.r) block_r[r] = 1;
    img.assign(NR, -1);
    std::vector<uint8_t> used(NP, 0);
    for (int r = 0; r < NR; ++r) {
        if (raw[r] >= 0 && !block_r[r]) {
            img[r] = raw[r];
            used[raw[r]] = 1;
        }
    }
    for (const auto& b : blocks) {
        if (b.r.size() > b.p.size()) return false;
        for (int r : b.r) {
            if (raw[r] < 0) continue;
            int p = raw[r];
            if (!contains(b.p, p) || used[p]) return false;
            img[r] = p;
            used[p] = 1;
        }
        std::vector<int> available;
        for (int p : b.p)
            if (!used[p]) available.push_back(p);
        std::vector<int> missing;
        for (int r : b.r)
            if (img[r] < 0) missing.push_back(r);
        if (available.size() < missing.size()) return false;
        for (size_t i = 0; i < missing.size(); ++i) {
            img[missing[i]] = available[i];
            used[available[i]] = 1;
        }
    }
    return true;
}

// _SymCand.__init__: validate + back-fill the witness inside blocks; nullopt
// stands for the Python ValueError (callers return None).
bool make_cand(const std::vector<int>& raw, const std::vector<Block>& blocks,
               std::vector<int> exact_fixed, int NP, Cand& out) {
    const int NR = (int)raw.size();
    std::vector<int> img;
    if (!pool_witness(raw, blocks, NP, img)) return false;
    out.img = std::move(img);
    out.mapped.clear();
    for (int r = 0; r < NR; ++r)
        if (out.img[r] >= 0) out.mapped.push_back(r);
    out.blocks = blocks;
    sort_unique(exact_fixed);
    out.exact_fixed = std::move(exact_fixed);
    return true;
}

#include "fragment_query.h"

struct GrowthObservation {
    size_t count = 0;
    Cand witness, alternative;
    std::vector<uint8_t> fragment;
    std::vector<Pair> deferred;
    void record(const std::vector<Cand>& candidates, const std::vector<uint8_t>& atoms,
                const std::vector<Pair>& boundary) {
        if (candidates.empty() || candidates.size() <= count) return;
        count = candidates.size(); witness = candidates.front();
        if (count > 1) alternative = candidates[1];
        fragment = atoms; deferred = boundary;
    }
};

// ---------------------------------------------------------------------------
// growth loop (mappa.growth.island.grow_island)
// ---------------------------------------------------------------------------

struct HeapItem { double neg_w; int u, n; };

struct GrowthRandom {
    uint64_t state;
    size_t below(size_t n) {
        const uint64_t threshold = -uint64_t(n) % n;
        for (;;) {
            uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15));
            z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
            z ^= z >> 31;
            if (z >= threshold) return z % n;
        }
    }
};

struct GrowProfile {
    long feasibility_calls=0,feasibility_nodes=0,feasibility_searches=0,completion_queries=0,completion_pruned=0;
    long completion_generated=0,completion_peak_candidates=0;
    double completion_seconds=0;
    long seed_targets = 0, init_cands = 0, initial_heap = 0, heap_pops = 0,
         stale_pops = 0, fragment_skip_pops = 0, extend_calls = 0, max_cands_before = 0,
         max_cands_after = 0, max_fragment_size = 0, max_heap_len = 0, commits = 0,
         deferred = 0, merge_calls = 0, free_extend_calls = 0,
         final_cands = 0, final_fragment_size = 0, branches = 0;
    std::string result;
};

struct BranchCap {
    long count, limit;
};

struct IsoOut {
    std::vector<Pair> mapping;             // sorted (r, p)
    std::vector<Pair> deferred_edges;      // sorted
    std::vector<int> fragment;             // sorted
    std::vector<Pair> preserved_bonds;     // accepted growth edges
    Cand cand;
};

struct GrowResult {
    std::vector<IsoOut> isos;
    GrowProfile profile;
    bool capped = false;
    BranchCap cap{0, 0};
};

#include "fragment_completion.h"

void finish_fragment(const Source& R, const Target& P, const std::vector<int>& mapping,
                     std::vector<Cand> cands, const std::vector<uint8_t>& fragment,
                     const std::vector<Pair>& deferred, const std::vector<Pair>& preserved_bonds,
                     double iso_tol, int min_lock_size,
                     long max_branches,
                     GrowResult& result) {
    auto& prof=result.profile;
    int fragment_size=std::count(fragment.begin(),fragment.end(),uint8_t{1});
    prof.final_fragment_size=fragment_size;
    prof.final_cands=(long)cands.size();
    if (cands.empty() || fragment_size<min_lock_size) {
        prof.result=cands.empty() ? "no_cands" : "too_small";return;
    }
    bool set_unique=true;
    for (const auto& c:cands) if (c.has_open_choice()) {set_unique=false;break;}
    if (set_unique && cands.size()>1)
        for (size_t i=1;i<cands.size();++i)
            if (cands[i].img!=cands[0].img) {set_unique=false;break;}
    if (!set_unique && (long)cands.size()>max_branches) {
        prof.result="subtree_branch_cap";prof.branches=(long)cands.size();
        result.capped=true;result.cap={(long)cands.size(),max_branches};return;
    }
    prof.result=set_unique ? "success" : "branched";
    prof.branches=set_unique ? 1 : (long)cands.size();
    for (size_t i=0;i<(set_unique ? 1 : cands.size());++i) {
        IsoOut iso;
        for (int r:cands[i].mapped) iso.mapping.push_back({r,cands[i].img[r]});
        iso.deferred_edges=deferred;
        iso.preserved_bonds=preserved_bonds;
        std::sort(iso.preserved_bonds.begin(), iso.preserved_bonds.end());
        for (int r=0;r<R.n;++r) if (fragment[r]) iso.fragment.push_back(r);
        iso.cand=std::move(cands[i]);result.isos.push_back(std::move(iso));
    }
}

#include "growth_dependencies.h"

GrowResult grow_island(const Source& R, const Target& P, int seed, const std::vector<int>& mapping,
                       double graph_floor, double iso_tol, int min_lock_size, long max_branches,
                       const std::vector<Pair>* islands, const std::vector<Pair>& prior_deferred,
                       bool allow_mapped_seed, uint64_t random_seed = 42,
                       GrowthObservation* observation = nullptr,
                       GrowthDependencies* dependencies = nullptr,
                       SourceReads* source_reads = nullptr) {
    SourceReadScope dependency_scope(R, source_reads);
    GrowResult result;
    GrowProfile& prof = result.profile;
    const int NR = R.n, NP = P.n;
    std::vector<uint8_t> locked_p(NP, 0);
    for (int r = 0; r < NR; ++r)
        if (mapping[r] >= 0) locked_p[mapping[r]] = 1;
    std::vector<int> island_of(NR, -1);
    if (islands)
        for (const auto& item : *islands) island_of[item.first] = item.second;

    if (dependencies) dependencies->mapping_at(seed);
    if (mapping[seed]>=0 && !allow_mapped_seed) {prof.result="already_mapped";return result;}
    WitnessGrowth witness(R,P,mapping,locked_p,iso_tol,seed);
    prof.seed_targets=witness.seed_targets;
    const long initial_count=witness.seed_targets ? 1 : 0;
    if (!initial_count) {
        prof.result = "no_initial_cands";
        return result;
    }
    std::vector<uint8_t> fragment(NR, 0);
    fragment[seed] = 1;
    int fragment_size = 1;
    std::vector<uint8_t> used_edge(NR * (size_t)NR, 0);
    auto edge_used = [&](int a, int b) -> uint8_t& { return used_edge[(size_t)std::min(a, b) * NR + std::max(a, b)]; };
    std::vector<Pair> preserved_bonds;
    std::vector<Pair> deferred = prior_deferred;
    for (auto& e : deferred)
        if (e.first > e.second) std::swap(e.first, e.second);
    std::sort(deferred.begin(), deferred.end());
    deferred.erase(std::unique(deferred.begin(), deferred.end()), deferred.end());

    std::vector<HeapItem> heap;
    GrowthRandom rng{random_seed};
    auto push_edges_from = [&](int atom) {
        for (int nb : R.neighbors(atom)) {
            if (fragment[nb]) continue;
            if (edge_used(atom, nb)) continue;
            double w = R.wbo(atom, nb);
            if (w >= graph_floor) heap.push_back({-w, atom, nb});
        }
    };
    push_edges_from(seed);
    prof.init_cands = initial_count;
    prof.initial_heap = (long)heap.size();
    prof.max_cands_before = prof.max_cands_after = initial_count;
    prof.max_fragment_size = 1;
    prof.max_heap_len = (long)heap.size();
    // Embeddings during growth are not completed branches. Apply the branch
    // budget only after direct completion in finish_fragment.

    if (observation) observation->record({witness.observation_candidate()},fragment,deferred);

    while (!heap.empty()) {
        // Sample distinct boundary atoms uniformly, then one incident edge.
        // Neither bond order nor the number of incident edges biases the atom.
        heap.erase(std::remove_if(heap.begin(), heap.end(), [&](const auto& e) {
            return fragment[e.n] || edge_used(e.u,e.n);
        }), heap.end());
        if (heap.empty()) break;
        std::vector<int> atoms;
        for (const auto& e:heap) atoms.push_back(e.n);
        std::sort(atoms.begin(),atoms.end());
        atoms.erase(std::unique(atoms.begin(),atoms.end()),atoms.end());
        int selected = atoms[rng.below(atoms.size())];
        std::vector<size_t> edges;
        for (size_t i=0;i<heap.size();++i) if (heap[i].n==selected) edges.push_back(i);
        std::sort(edges.begin(),edges.end(),[&](size_t a,size_t b){return heap[a].u<heap[b].u;});
        size_t index = edges[rng.below(edges.size())];
        HeapItem item = heap[index];
        heap.erase(heap.begin()+index);
        prof.heap_pops++;
        int u = item.u, n = item.n;
        double wbo = -item.neg_w;
        if (edge_used(u, n)) { prof.stale_pops++; continue; }
        edge_used(u, n) = 1;
        if (fragment[n]) { prof.fragment_skip_pops++; continue; }
        bool n_in_mapping = mapping[n] >= 0;
        if (dependencies) {
            dependencies->mapping_at(n);
            if (n_in_mapping) dependencies->island_at(n, islands, island_of);
        }
        prof.extend_calls++;
        if (n_in_mapping) prof.merge_calls++; else prof.free_extend_calls++;
        long old_count = 1;
        prof.max_cands_before = std::max(prof.max_cands_before, old_count);

        bool matched=witness.extend(u,n,wbo,fragment,preserved_bonds,islands,island_of);
        if (matched) {
            preserved_bonds.push_back({std::min(u,n), std::max(u,n)});
            prof.commits++;
            prof.max_cands_after = std::max(prof.max_cands_after, 1L);
            std::vector<int> added{n};
            if (n_in_mapping && islands && island_of[n] >= 0) {
                int iid = island_of[n];
                for (const auto& it : *islands)
                    if (it.second == iid && !fragment[it.first] && it.first != n) added.push_back(it.first);
            }
            for (int a : added) {
                if (!fragment[a]) {
                    fragment[a] = 1;
                    fragment_size++;
                }
            }
            for (int a : added) push_edges_from(a);
            prof.max_fragment_size = std::max(prof.max_fragment_size, (long)fragment_size);
            prof.max_heap_len = std::max(prof.max_heap_len, (long)heap.size());
        } else {
            prof.deferred++;
            Pair e{std::min(u, n), std::max(u, n)};
            if (dependencies) dependencies->local_deferred.push_back(e);
            auto pos = std::lower_bound(deferred.begin(), deferred.end(), e);
            if (pos == deferred.end() || *pos != e) deferred.insert(pos, e);
        }
    }
    prof.feasibility_calls=witness.calls;prof.feasibility_nodes=witness.nodes;
    prof.feasibility_searches=witness.searches;
    auto cands=complete_fragment(R,P,mapping,fragment,preserved_bonds,iso_tol,prof);
    if (observation) observation->record(cands,fragment,deferred);
    if (dependencies) dependencies->boundary(fragment);
    finish_fragment(R,P,mapping,std::move(cands),fragment,deferred,preserved_bonds,iso_tol,min_lock_size,
                    max_branches,result);
    return result;
}

// ---------------------------------------------------------------------------
// Python conversions
// ---------------------------------------------------------------------------

std::unordered_map<std::string, int>& element_table() {
    static std::unordered_map<std::string, int> table;
    return table;
}

int element_code(const std::string& e) {
    auto& table = element_table();
    auto it = table.find(e);
    if (it == table.end()) it = table.emplace(e, (int)table.size()).first;
    return it->second;
}

void fill_graph(Graph& g, const std::vector<std::string>& elements,
                const std::vector<std::vector<double>>& wbo, double bond_cut,
                const std::vector<Pair>& edges) {
    g.n = (int)elements.size();
    g.elem = elements;
    g.ecode.resize(g.n);
    for (int i = 0; i < g.n; ++i) g.ecode[i] = element_code(elements[i]);
    g.adj.assign(g.n, {});
    g.edge.assign((size_t)g.n * g.n, 0);
    g.w.assign((size_t)g.n * g.n, 0.0);
    for (int i = 0; i < g.n; ++i)
        for (int j = 0; j < g.n; ++j) g.w[(size_t)i * g.n + j] = wbo[i][j];
    g.bond_cut = bond_cut;
    for (const auto& e : edges) {
        g.edge[(size_t)e.first * g.n + e.second] = 1;
        g.edge[(size_t)e.second * g.n + e.first] = 1;
        g.adj[e.first].push_back(e.second);
        g.adj[e.second].push_back(e.first);
    }
    for (auto& list : g.adj) std::sort(list.begin(), list.end());
}

struct PySource {
    Source g;
    PySource(const std::vector<std::string>& elements, const std::vector<std::vector<double>>& wbo,
             double bond_cut, const std::vector<Pair>& edges,
             const std::vector<std::tuple<int, int, double>>& tolerances, bool inclusive) {
        fill_graph(g, elements, wbo, bond_cut, edges);
        g.tolerance_inclusive = inclusive;
        if (!tolerances.empty()) {
            g.bond_tolerances.assign((size_t)g.n * g.n, -1.0);
            for (auto [a,b,value] : tolerances) {
                if (a < 0 || b < 0 || a >= g.n || b >= g.n || !std::isfinite(value) || value < 0)
                    throw std::invalid_argument("invalid source bond tolerance");
                g.bond_tolerances[(size_t)a*g.n+b] = g.bond_tolerances[(size_t)b*g.n+a] = value;
            }
            for (auto [a,b] : edges)
                if (g.tolerance(a,b,0) < 0) throw std::invalid_argument("missing source bond tolerance");
        }
    }
};

struct PyTarget {
    Target g;
    PyTarget(const std::vector<std::string>& elements, const std::vector<std::vector<double>>& wbo,
             double bond_cut, const std::vector<Pair>& edges, const std::vector<int>& orbits) {
        fill_graph(g, elements, wbo, bond_cut, edges);
        g.orbit = orbits;
        g.component.assign(g.n, -1);
        for (int root = 0; root < g.n; ++root) {
            if (g.component[root] >= 0) continue;
            std::vector<int> pending{root};
            g.component[root] = root;
            for (size_t i = 0; i < pending.size(); ++i)
                for (int neighbor : g.adj[pending[i]])
                    if (g.component[neighbor] < 0) {
                        g.component[neighbor] = root;
                        pending.push_back(neighbor);
                    }
        }
        // Preserve the full weighted matrix and active adjacency, including
        // weak pairs below the bond floor. No tolerance buckets prove symmetry.
        std::vector<Pair> nauty_edges;
        int next = g.n;
        for (int a = 0; a < g.n; ++a) for (int b = a + 1; b < g.n; ++b) {
            const double weight = g.wbo(a, b);
            if (weight == 0.0 && !g.has_edge(a, b)) continue;
            nauty_edges.push_back({a, next});
            nauty_edges.push_back({b, next});
            int64_t bits;
            std::memcpy(&bits, &weight, sizeof(bits));
            g.edge_cells[{bits, g.has_edge(a, b)}].push_back(next);
            ++next;
        }
        g.n_vertices = next;
        g.nauty = std::make_unique<DenseGraph>(next, nauty_edges);
        int max_code = 0;
        for (int c : g.ecode) max_code = std::max(max_code, c);
        g.same_element.assign((size_t)max_code + 1, {});
        for (int v = 0; v < g.n; ++v) g.same_element[g.ecode[v]].push_back(v);
    }
};

py::dict block_dict(const Block& b) {
    py::dict d;
    py::list r, p;
    for (int x : b.r) r.append(x);
    for (int x : b.p) p.append(x);
    d["r_atoms"] = r;
    d["p_atoms"] = p;
    d["extendable"] = b.extendable;
    d["open"] = b.open();
    size_t n = b.p.size(), k = b.r.size();
    std::string expr;
    if (k <= 0 || n <= 1) expr = "1";
    else if (k == 1) expr = std::to_string(n);
    else if (k == n) expr = std::to_string(n) + "!";
    else expr = "P(" + std::to_string(n) + "," + std::to_string(k) + ")";
    d["assignments"] = expr;
    return d;
}

py::dict symmetry_state(const Cand& c) {
    py::dict item;
    py::dict witness;
    for (int r : c.mapped) witness[py::int_(r)] = c.img[r];
    item["witness"] = witness;
    item["compression_version"] = 3;
    py::list blocks;
    item["blocks"] = blocks;
    py::list exact;
    for (int x : c.exact_fixed) exact.append(x);
    item["exact_fixed"] = exact;
    for (const auto& b : c.blocks) blocks.append(block_dict(b));
    return item;
}

py::object grow_result_dict(const GrowResult& res) {
    py::dict out;
    py::list isos;
    for (const auto& iso : res.isos) {
        py::dict d;
        py::dict m;
        for (const auto& pr : iso.mapping) m[py::int_(pr.first)] = pr.second;
        d["mapping"] = m;
        d["deferred_edges"] = iso.deferred_edges;
        d["preserved_bonds"] = iso.preserved_bonds;
        d["fragment"] = iso.fragment;
        d["symmetry"] = symmetry_state(iso.cand);
        isos.append(d);
    }
    out["isos"] = isos;
    out["capped"] = res.capped;
    out["cap_count"] = res.cap.count;
    out["cap_limit"] = res.cap.limit;
    py::dict prof;
    const GrowProfile& p = res.profile;
    prof["feasibility_calls"]=p.feasibility_calls;prof["feasibility_nodes"]=p.feasibility_nodes;
    prof["feasibility_searches"]=p.feasibility_searches;
    prof["completion_seconds"]=p.completion_seconds;
    prof["completion_queries"]=p.completion_queries;prof["completion_pruned"]=p.completion_pruned;
    prof["completion_generated"]=p.completion_generated;
    prof["completion_peak_candidates"]=p.completion_peak_candidates;
    prof["seed_targets"] = p.seed_targets;
    prof["init_cands"] = p.init_cands;
    prof["initial_heap"] = p.initial_heap;
    prof["heap_pops"] = p.heap_pops;
    prof["stale_pops"] = p.stale_pops;
    prof["fragment_skip_pops"] = p.fragment_skip_pops;
    prof["extend_calls"] = p.extend_calls;
    prof["max_cands_before"] = p.max_cands_before;
    prof["max_cands_after"] = p.max_cands_after;
    prof["max_fragment_size"] = p.max_fragment_size;
    prof["max_heap_len"] = p.max_heap_len;
    prof["commits"] = p.commits;
    prof["deferred"] = p.deferred;
    prof["merge_calls"] = p.merge_calls;
    prof["free_extend_calls"] = p.free_extend_calls;
    prof["final_cands"] = p.final_cands;
    prof["final_fragment_size"] = p.final_fragment_size;
    prof["branches"] = p.branches;
    prof["result"] = p.result;
    out["profile"] = prof;
    return out;
}

py::object py_grow_island(const PySource& source, const PyTarget& target, int seed,
                          const std::vector<int>& mapping, double graph_floor, double iso_tol,
                          int min_lock_size, long max_branches, py::object islands_obj,
                          const std::vector<Pair>& prior_deferred, bool allow_mapped_seed, uint64_t random_seed, bool observe) {
    std::vector<Pair> islands;
    const std::vector<Pair>* islands_ptr = nullptr;
    if (!islands_obj.is_none()) {
        islands = islands_obj.cast<std::vector<Pair>>();
        islands_ptr = &islands;
    }
    GrowResult res;
    GrowthObservation observation;
    {
        py::gil_scoped_release release;
        res = grow_island(source.g, target.g, seed, mapping, graph_floor, iso_tol, min_lock_size,
                          max_branches, islands_ptr, prior_deferred, allow_mapped_seed,
                          random_seed,observe ? &observation : nullptr);
    }
    auto out = grow_result_dict(res).cast<py::dict>();
    if (observe) {
        py::dict snapshot;
        snapshot["candidate_count"] = observation.count;
        snapshot["symmetry"] = symmetry_state(observation.witness);
        if (observation.count > 1)
            snapshot["alternative"] = symmetry_state(observation.alternative);
        std::vector<int> atoms;
        for (size_t i=0;i<observation.fragment.size();++i)
            if (observation.fragment[i]) atoms.push_back(i);
        snapshot["fragment"] = atoms;
        snapshot["deferred"] = observation.deferred;
        out["observation"] = snapshot;
    }
    return out;
}

py::dict py_find_fragment_witness(const PySource& source, const PyTarget& target,
                                 const std::vector<int>& atoms, const std::vector<Pair>& bonds,
                                 const std::vector<int>& locked, double tolerance) {
    if (locked.size() != (size_t)source.g.n)
        throw std::invalid_argument("locked mapping must cover the source indices");
    std::vector<uint8_t> fragment(source.g.n, 0);
    for (int atom : atoms) fragment.at(atom) = 1;
    for (auto [a, b] : bonds)
        if (!fragment.at(a) || !fragment.at(b))
            throw std::invalid_argument("accepted bond lies outside the fragment");
    for (int image : locked)
        if (image < -1 || image >= target.g.n)
            throw std::invalid_argument("locked image lies outside the target");
    std::vector<int> witness;
    long nodes;
    {
        py::gil_scoped_release release;
        FragmentQuery query(source.g, target.g, fragment, bonds, locked, tolerance);
        witness = query.find();
        nodes = query.nodes;
    }
    py::dict result;
    result["witness"] = py::none();
    if (!witness.empty()) {
        py::dict mapping;
        for (int atom : atoms) mapping[py::int_(atom)] = witness[atom];
        result["witness"] = mapping;
    }
    result["nodes"] = nodes;
    return result;
}

#include "fragment_repair.h"
#include "fragment_scheduler.h"
#include "fragment_scheduler_binding.h"

}  // namespace

// native/src/autgrp.cpp
void register_autgrp(pybind11::module_&);

PYBIND11_MODULE(_engine, mod) {
    mod.doc() = "Native AAM growth engine";
    register_autgrp(mod);
    py::class_<PySource>(mod, "SourceGraph")
        .def(py::init<const std::vector<std::string>&, const std::vector<std::vector<double>>&,
                      double, const std::vector<Pair>&,
                      const std::vector<std::tuple<int,int,double>>&, bool>(),
             py::arg("elements"), py::arg("wbo"), py::arg("bond_cut"), py::arg("edges"),
             py::arg("tolerances")=std::vector<std::tuple<int,int,double>>{},
             py::arg("inclusive")=true);
    py::class_<PyTarget>(mod, "TargetGraph")
        .def(py::init<const std::vector<std::string>&, const std::vector<std::vector<double>>&,
                      double, const std::vector<Pair>&, const std::vector<int>&>())
        .def("symmetry_stats", [](const PyTarget& t) {
            py::dict result;
            result["generator_calls"] = t.g.nauty->generator_calls;
            return result;
        });
    mod.def("grow_island", &py_grow_island, py::arg("source"), py::arg("target"), py::arg("seed"),
            py::arg("mapping"), py::arg("graph_floor"), py::arg("iso_tol"),
            py::arg("min_lock_size"), py::arg("max_branches"), py::arg("islands"),
            py::arg("prior_deferred_edges"), py::arg("allow_mapped_seed"), py::arg("random_seed")=42, py::arg("observe")=false);
    mod.def("find_fragment_witness", &py_find_fragment_witness,
            py::arg("source"), py::arg("target"), py::arg("atoms"), py::arg("bonds"),
            py::arg("locked"), py::arg("tolerance"));
    register_fragment_repair(mod);
    register_fragment_scheduler(mod);
}
