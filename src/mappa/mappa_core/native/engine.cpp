// Native compressed matching and fragment growth. Python matching provides an
// independent scientific oracle. Element matching uses this backend by default;
// custom node policies and diagnostic traces use the Python implementation.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
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

#include "autgrp.h"
#include "set_orbit.h"
#include <numeric>

namespace py = pybind11;

namespace {

using Pair = std::pair<int, int>;
using PlacementGenerators=permutation_native::Generators;
using PlacementGroup=std::shared_ptr<const PlacementGenerators>;

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

struct PreparedSymmetry {
    std::vector<int> atoms, anchors;
    std::vector<std::string> elements;
    double bond_threshold=0.2;
    std::vector<double> weights;
    std::vector<std::vector<int>> generators;
};

struct Target : Graph {
    std::vector<std::vector<int>> same_element;
    mutable std::shared_ptr<const PreparedSymmetry> global_self;
    mutable PlacementGroup placement_generators;
    mutable std::unique_ptr<permutation_native::SetOrbit> placement_orbits;
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

inline bool contains(const std::vector<int>& sorted, int x) {
    return std::binary_search(sorted.begin(), sorted.end(), x);
}

#include "fragment_query.h"

struct GrowthObservation {
    size_t count = 0;
    std::vector<int> witness, alternative;
    std::vector<uint8_t> fragment;
    std::vector<Pair> deferred;
    void record(const std::vector<std::vector<int>>& candidates, const std::vector<uint8_t>& atoms,
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

struct FragmentGrowth {
    std::vector<uint8_t> source_atoms;
    std::vector<Pair> accepted_edges, deferred_edges;
    std::vector<HeapItem> frontier;
    WitnessGrowth witness;
    GrowthRandom random;
    int atom_count = 0;

    FragmentGrowth(const Source& source, const Target& target,
                   const std::vector<int>& mapping, const std::vector<uint8_t>& occupied,
                   double tolerance, int seed, const std::vector<Pair>& deferred,
                   uint64_t random_seed)
        :witness(source,target,mapping,occupied,tolerance,seed),random{random_seed} {
        if (!witness.seed_targets) return;
        source_atoms.assign(source.n,0);
        source_atoms[seed]=1;atom_count=1;
        examined_edges.assign((size_t)source.n*source.n,0);
        deferred_edges=deferred;
        for (auto& [a,b]:deferred_edges) if (a>b) std::swap(a,b);
        std::sort(deferred_edges.begin(),deferred_edges.end());
        deferred_edges.erase(std::unique(deferred_edges.begin(),deferred_edges.end()),deferred_edges.end());
    }

    uint8_t& edge_used(int a, int b) {
        return examined_edges[(size_t)std::min(a,b)*source_atoms.size()+std::max(a,b)];
    }

    void queue_edges(const Source& source, int atom, double graph_floor) {
        for (int neighbor:source.neighbors(atom)) {
            if (source_atoms[neighbor] || edge_used(atom,neighbor)) continue;
            double weight=source.wbo(atom,neighbor);
            if (weight>=graph_floor) frontier.push_back({-weight,atom,neighbor});
        }
    }

private:
    std::vector<uint8_t> examined_edges;
};

struct GrowProfile {
    long feasibility_calls=0,feasibility_nodes=0,feasibility_searches=0,completion_queries=0,completion_pruned=0;
    long completion_generated=0,completion_peak_candidates=0;
    long completion_global_generators=0,completion_global_pruned=0;
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

struct Fragment {
    std::shared_ptr<const PreparedSymmetry> source, target;
    std::vector<std::vector<int>> domains;
    std::vector<int> witness;
    std::vector<Pair> anchors, accepted_edges, target_edges;
    std::vector<double> edge_tolerances;
    bool inclusive=true;
    double target_bond_floor=0.2;
};

struct GrowResult {
    std::vector<std::shared_ptr<const Fragment>> fragments;
    PlacementGroup placement_generators=std::make_shared<const PlacementGenerators>();
    std::vector<Pair> deferred_edges;
    GrowProfile profile;
    bool capped=false;
    BranchCap cap{0,0};
};

#include "fragment_completion.h"

#include "growth_dependencies.h"

GrowResult grow_island(const Source& R, const Target& P, int seed, const std::vector<int>& mapping,
                       double graph_floor, double iso_tol, int min_lock_size, long max_branches,
                       const std::vector<Pair>* islands, const std::vector<Pair>& prior_deferred,
                       bool allow_mapped_seed, uint64_t random_seed = 42,
                       GrowthObservation* observation = nullptr,
                       GrowthDependencies* dependencies = nullptr,
                       SourceReads* source_reads = nullptr,
                       const std::vector<int>& placement_fixed_targets = {}) {
    SourceReadScope dependency_scope(R, source_reads);
    GrowResult result;
    GrowProfile& prof = result.profile;
    const int NR = R.n, NP = P.n;
    for (int p:placement_fixed_targets) if (p<0 || p>=NP)
        throw std::invalid_argument("Placement fixed target is outside the target graph");
    std::vector<uint8_t> locked_p(NP, 0);
    for (int r = 0; r < NR; ++r)
        if (mapping[r] >= 0) locked_p[mapping[r]] = 1;
    std::vector<int> island_of(NR, -1);
    if (islands)
        for (const auto& item : *islands) island_of[item.first] = item.second;

    if (dependencies) dependencies->mapping_at(seed);
    if (mapping[seed]>=0 && !allow_mapped_seed) {prof.result="already_mapped";return result;}
    FragmentGrowth growth(R,P,mapping,locked_p,iso_tol,seed,prior_deferred,random_seed);
    prof.seed_targets=growth.witness.seed_targets;
    const long initial_count=growth.witness.seed_targets ? 1 : 0;
    if (!initial_count) {
        prof.result = "no_initial_cands";
        return result;
    }
    growth.queue_edges(R,seed,graph_floor);
    prof.init_cands = initial_count;
    prof.initial_heap = (long)growth.frontier.size();
    prof.max_cands_before = prof.max_cands_after = initial_count;
    prof.max_fragment_size = 1;
    prof.max_heap_len = (long)growth.frontier.size();
    // Embeddings during growth are not completed branches. Apply the branch
    // budget only after distinct-image completion.

    if (observation) observation->record({growth.witness.assignment()},growth.source_atoms,growth.deferred_edges);

    while (!growth.frontier.empty()) {
        // Sample distinct boundary atoms uniformly, then one incident edge.
        // Neither bond order nor the number of incident edges biases the atom.
        growth.frontier.erase(std::remove_if(growth.frontier.begin(), growth.frontier.end(), [&](const auto& e) {
            return growth.source_atoms[e.n] || growth.edge_used(e.u,e.n);
        }), growth.frontier.end());
        if (growth.frontier.empty()) break;
        std::vector<int> atoms;
        for (const auto& e:growth.frontier) atoms.push_back(e.n);
        std::sort(atoms.begin(),atoms.end());
        atoms.erase(std::unique(atoms.begin(),atoms.end()),atoms.end());
        int selected = atoms[growth.random.below(atoms.size())];
        std::vector<size_t> edges;
        for (size_t i=0;i<growth.frontier.size();++i) if (growth.frontier[i].n==selected) edges.push_back(i);
        std::sort(edges.begin(),edges.end(),[&](size_t a,size_t b){return growth.frontier[a].u<growth.frontier[b].u;});
        size_t index = edges[growth.random.below(edges.size())];
        HeapItem item = growth.frontier[index];
        growth.frontier.erase(growth.frontier.begin()+index);
        prof.heap_pops++;
        int u = item.u, n = item.n;
        double wbo = -item.neg_w;
        if (growth.edge_used(u, n)) { prof.stale_pops++; continue; }
        growth.edge_used(u, n) = 1;
        if (growth.source_atoms[n]) { prof.fragment_skip_pops++; continue; }
        bool n_in_mapping = mapping[n] >= 0;
        if (dependencies) {
            dependencies->mapping_at(n);
            if (n_in_mapping) dependencies->island_at(n, islands, island_of);
        }
        prof.extend_calls++;
        if (n_in_mapping) prof.merge_calls++; else prof.free_extend_calls++;
        long old_count = 1;
        prof.max_cands_before = std::max(prof.max_cands_before, old_count);

        bool matched=growth.witness.extend(u,n,wbo,growth.source_atoms,growth.accepted_edges,islands,island_of);
        if (matched) {
            growth.accepted_edges.push_back({std::min(u,n), std::max(u,n)});
            prof.commits++;
            prof.max_cands_after = std::max(prof.max_cands_after, 1L);
            std::vector<int> added{n};
            if (n_in_mapping && islands && island_of[n] >= 0) {
                int iid = island_of[n];
                for (const auto& it : *islands)
                    if (it.second == iid && !growth.source_atoms[it.first] && it.first != n) added.push_back(it.first);
            }
            for (int a : added) {
                if (!growth.source_atoms[a]) {
                    growth.source_atoms[a] = 1;
                    growth.atom_count++;
                }
            }
            for (int a : added) growth.queue_edges(R,a,graph_floor);
            prof.max_fragment_size = std::max(prof.max_fragment_size, (long)growth.atom_count);
            prof.max_heap_len = std::max(prof.max_heap_len, (long)growth.frontier.size());
        } else {
            prof.deferred++;
            Pair e{std::min(u, n), std::max(u, n)};
            if (dependencies) dependencies->local_deferred.push_back(e);
            auto pos = std::lower_bound(growth.deferred_edges.begin(), growth.deferred_edges.end(), e);
            if (pos == growth.deferred_edges.end() || *pos != e) growth.deferred_edges.insert(pos, e);
        }
    }
    prof.feasibility_calls=growth.witness.calls;prof.feasibility_nodes=growth.witness.nodes;
    prof.feasibility_searches=growth.witness.searches;
    if (dependencies) dependencies->boundary(growth.source_atoms);
    complete_fragment(R,P,mapping,growth,graph_floor,iso_tol,
                      min_lock_size,max_branches,result,observation,placement_fixed_targets);
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
             double bond_cut, const std::vector<Pair>& edges) {
        fill_graph(g, elements, wbo, bond_cut, edges);
        int max_code = 0;
        for (int c : g.ecode) max_code = std::max(max_code, c);
        g.same_element.assign((size_t)max_code + 1, {});
        for (int v = 0; v < g.n; ++v) g.same_element[g.ecode[v]].push_back(v);
    }
};

class FragmentExport {
    std::unordered_map<const PreparedSymmetry*,py::dict> prepared;
    std::unordered_map<const PlacementGenerators*,py::tuple> placements;
    std::map<std::vector<int>,py::tuple> generator_rows;
public:
    py::tuple placement(const PlacementGroup& value,const std::vector<int>& target_ids={}) {
        if (!value) return py::tuple();
        auto found=placements.find(value.get());
        if (found!=placements.end()) return found->second;
        py::tuple result(value->size());
        for (size_t i=0;i<value->size();++i) {
            auto translated=(*value)[i];
            if (!target_ids.empty()) {
                translated=permutation_native::identity(*std::max_element(target_ids.begin(),target_ids.end())+1);
                for (size_t p=0;p<target_ids.size();++p) translated[target_ids[p]]=target_ids[(*value)[i][p]];
            }
            auto row=generator_rows.find(translated);
            if (row==generator_rows.end()) row=generator_rows.emplace(translated,py::tuple(py::cast(translated))).first;
            result[i]=row->second;
        }
        placements.emplace(value.get(),result);
        return result;
    }
    py::dict self(const std::shared_ptr<const PreparedSymmetry>& value) {
        auto found=prepared.find(value.get());
        if (found!=prepared.end()) return found->second;
        py::dict out;
        out["self_id"]=prepared.size();
        out["atoms"]=value->atoms;out["anchors"]=value->anchors;
        out["elements"]=value->elements;
        out["bond_threshold"]=value->bond_threshold;
        out["generators"]=value->generators;
        const size_t n=value->atoms.size();
        py::array_t<double> weights({(py::ssize_t)n,(py::ssize_t)n});
        std::copy(value->weights.begin(),value->weights.end(),weights.mutable_data());
        weights.attr("setflags")(false);
        out["weights"]=weights;
        prepared.emplace(value.get(),out);
        return out;
    }
    py::dict relation(const Fragment& fragment) {
        py::dict out;
        out["source"]=self(fragment.source);out["target"]=self(fragment.target);
        out["domains"]=fragment.domains;out["witness"]=fragment.witness;
        out["anchors"]=fragment.anchors;out["accepted_edges"]=fragment.accepted_edges;
        out["edge_tolerances"]=fragment.edge_tolerances;out["inclusive"]=fragment.inclusive;
        out["target_edges"]=fragment.target_edges;out["target_bond_floor"]=fragment.target_bond_floor;
        return out;
    }
};

py::object grow_result_dict(const GrowResult& res) {
    py::dict out;
    FragmentExport exported;
    py::list fragments;
    for (const auto& fragment:res.fragments) fragments.append(exported.relation(*fragment));
    out["fragments"]=fragments;
    out["placement_generators"]=exported.placement(res.placement_generators);
    out["deferred_edges"] = res.deferred_edges;
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
    prof["completion_global_generators"]=p.completion_global_generators;
    prof["completion_global_pruned"]=p.completion_global_pruned;
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
                          const std::vector<Pair>& prior_deferred, bool allow_mapped_seed, uint64_t random_seed, bool observe,
                          const std::vector<int>& placement_fixed_targets) {
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
                          random_seed,observe ? &observation : nullptr,nullptr,nullptr,placement_fixed_targets);
    }
    auto out = grow_result_dict(res).cast<py::dict>();
    if (observe) {
        py::dict snapshot;
        snapshot["candidate_count"] = observation.count;
        auto assignments=[](const std::vector<int>& image) {
            py::dict out;
            for (size_t r=0;r<image.size();++r) if (image[r]>=0) out[py::int_(r)]=image[r];
            return out;
        };
        snapshot["witness"]=assignments(observation.witness);
        if (observation.count>1) snapshot["alternative"]=assignments(observation.alternative);
        std::vector<Pair> anchors;
        for (size_t r=0;r<observation.witness.size();++r)
            if (observation.witness[r]>=0 && mapping[r]>=0) anchors.push_back({(int)r,mapping[r]});
        snapshot["anchors"]=anchors;
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
                      double, const std::vector<Pair>&>());
    mod.def("grow_island", &py_grow_island, py::arg("source"), py::arg("target"), py::arg("seed"),
            py::arg("mapping"), py::arg("graph_floor"), py::arg("iso_tol"),
            py::arg("min_lock_size"), py::arg("max_branches"), py::arg("islands"),
            py::arg("prior_deferred_edges"), py::arg("allow_mapped_seed"), py::arg("random_seed")=42, py::arg("observe")=false,
            py::arg("placement_fixed_targets")=std::vector<int>{});
    mod.def("find_fragment_witness", &py_find_fragment_witness,
            py::arg("source"), py::arg("target"), py::arg("atoms"), py::arg("bonds"),
            py::arg("locked"), py::arg("tolerance"));
    register_fragment_repair(mod);
    register_fragment_scheduler(mod);
}
