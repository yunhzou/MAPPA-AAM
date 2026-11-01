// Reuse complete conditional fragment results across changed earlier histories.
// This is NOT witness-only repair: the whole compressed candidate result must
// have identical constraints. Different histories can pass that dependency test.
#include "fragment_cache_io.h"

class FragmentRepair {
    Source base, active_source;
    const Target& target;
    std::vector<Pair> active_cuts;
    size_t budget, resident=0, peak_resident=0;
    bool existence;
    std::unique_ptr<ExtensionCache> extensions;
    struct Entry {
        std::string key;
        GrowResult result;
        GrowthDependencies dependencies;
        std::vector<int> mapping;
        std::vector<Pair> islands, cuts, deferred, pair_reads;
        std::vector<uint8_t> row_reads;
        size_t bytes=0;
        bool in_snapshot=false;
        Entry(int n) : dependencies(n) {}
    };
    using Entries=std::list<Entry>;
    Entries lru;
    std::unordered_map<std::string,std::vector<Entries::iterator>> cache;
    long calls=0, hits=0, history_hits=0, topology_hits=0, boundary_rebases=0;
    long checks=0, evictions=0, oversized=0, mapping_misses=0, island_misses=0, topology_misses=0, boundary_misses=0;
    long logical_extensions=0, reused_extensions=0, logical_certificates=0, reused_certificates=0;
    uint64_t revision=0;

    static std::vector<Pair> normalized(std::vector<Pair> pairs) {
        for (auto& [a,b]:pairs) if (a>b) std::swap(a,b);
        std::sort(pairs.begin(),pairs.end());
        pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
        return pairs;
    }
    static std::vector<Pair> difference(const std::vector<Pair>& a,const std::vector<Pair>& b) {
        std::vector<Pair> out;
        std::set_symmetric_difference(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(out));
        return out;
    }
    template<class Value> static void append(std::string& key,const Value& value) {
        key.append(reinterpret_cast<const char*>(&value),sizeof(value));
    }
    static size_t entry_bytes(const Entry& e) {
        // Charge owned vector capacity, both key copies, and conservative
        // list/hash-index overhead. The budget does not include caller outputs.
        size_t out=sizeof(e)+2*(e.key.capacity()+1)+256+e.mapping.capacity()*sizeof(int)
            +e.row_reads.capacity()+e.dependencies.mapping_reads.capacity()+e.dependencies.boundary_atoms.capacity();
        out+=e.result.isos.capacity()*sizeof(IsoOut)
            +e.result.profile.result.capacity()+1+e.result.profile.growth_policy.capacity()+1
            +e.dependencies.island_reads.capacity()*sizeof(decltype(e.dependencies.island_reads)::value_type);
        for (const auto* v:{&e.islands,&e.cuts,&e.deferred,&e.pair_reads,&e.dependencies.local_deferred})
            out+=v->capacity()*sizeof(Pair);
        for (const auto& [r,members]:e.dependencies.island_reads)
            out+=members.capacity()*sizeof(int);
        for (const auto& iso:e.result.isos) {
            out+=(iso.mapping.capacity()+iso.deferred_edges.capacity()+iso.preserved_bonds.capacity())*sizeof(Pair)
                +iso.fragment.capacity()*sizeof(int);
            const auto& c=iso.cand;
            out+=(c.img.capacity()+c.mapped.capacity()+c.exact_fixed.capacity())*sizeof(int);
            for (const auto* blocks:{&c.blocks,&c.automorph}) {
                out+=blocks->capacity()*sizeof(Block);
                for (const auto& b:*blocks) out+=(b.r.capacity()+b.p.capacity())*sizeof(int);
            }
        }
        return out;
    }
    void insert(Entry entry) {
        while (!lru.empty() && resident+entry.bytes>budget) {
            auto old=lru.begin();
            auto bucket=cache.find(old->key);
            auto& entries=bucket->second;
            entries.erase(std::remove(entries.begin(),entries.end(),old),entries.end());
            if (entries.empty()) cache.erase(bucket);
            resident-=old->bytes;lru.erase(old);++evictions;
        }
        resident+=entry.bytes;
        peak_resident=std::max(peak_resident,resident);
        lru.push_back(std::move(entry));
        auto added=std::prev(lru.end());cache[added->key].push_back(added);
        ++revision;
    }
    bool valid(const Entry& entry,const std::vector<int>& mapping,
               const std::vector<Pair>& islands,const std::vector<Pair>& cuts,
               const std::vector<Pair>& deferred) {
        ++checks;
        for (auto [a,b]:difference(entry.cuts,cuts)) {
            if (entry.row_reads[a] || entry.row_reads[b] ||
                std::binary_search(entry.pair_reads.begin(),entry.pair_reads.end(),Pair{a,b})) {
                ++topology_misses;return false;
            }
        }
        for (int r=0;r<base.n;++r) if (entry.dependencies.mapping_reads[r] && entry.mapping[r]!=mapping[r]) {
            ++mapping_misses;return false;
        }
        for (const auto& [r,old_members]:entry.dependencies.island_reads) {
            int label=-1;
            for (auto [atom,iid]:islands) if (atom==r) {label=iid;break;}
            std::vector<int> members;
            if (label>=0) for (auto [atom,iid]:islands) if (iid==label) members.push_back(atom);
            if (members!=old_members) {++island_misses;return false;}
        }
        for (auto [a,b]:difference(entry.deferred,deferred))
            if (entry.dependencies.boundary_atoms[a] || entry.dependencies.boundary_atoms[b]) {
                ++boundary_misses;return false;
            }
        return true;
    }
public:
    py::bytes snapshot_context() const {
        fragment_cache_io::Writer out;
        out(std::string("mappa.fragment-cache/1"),uint64_t(sizeof(long)),existence);
        for (const Graph* g:{static_cast<const Graph*>(&base),static_cast<const Graph*>(&target)})
            // Element codes are process-local interned IDs; names are portable.
            out(g->n,g->elem,g->adj,g->edge,g->w,g->bond_cut,g->bond_tolerances,g->tolerance_inclusive);
        out(target.orbit,target.bucket,target.zero_bucket);
        return py::bytes(out.data);
    }
    py::bytes dump(const std::string& identity) const {
        fragment_cache_io::Writer out;
        out(identity,uint64_t(lru.size()));
        for (const auto& e:lru)
            out(e.key,e.result,e.dependencies,e.mapping,e.islands,e.cuts,e.deferred,e.pair_reads,e.row_reads);
        return py::bytes(out.data);
    }
    bool load(const std::string& data,const std::string& identity) {
        fragment_cache_io::Reader in(data);
        std::string observed; in(observed);
        if (observed!=identity) throw std::invalid_argument("Fragment cache context differs");
        auto count=in.size();
        Entries parsed;
        for (uint64_t i=0;i<count;++i) {
            Entry e(base.n);
            in(e.key,e.result,e.dependencies,e.mapping,e.islands,e.cuts,e.deferred,e.pair_reads,e.row_reads);
            if (e.mapping.size()!=size_t(base.n) || e.row_reads.size()!=size_t(base.n)
                || e.dependencies.mapping_reads.size()!=size_t(base.n) || e.dependencies.boundary_atoms.size()!=size_t(base.n))
                throw std::invalid_argument("Fragment cache dependency size differs");
            e.bytes=entry_bytes(e);
            parsed.push_back(std::move(e));
        }
        in.finish(); // Failed reads must leave the existing cache usable.
        // Import knowledge without discarding this worker's unpublished results.
        // Identical requests share an entry; compressed answers are never merged.
        for (auto& e:lru) e.in_snapshot=false;
        for (auto& e:parsed) {
            auto bucket=cache.find(e.key);
            bool present=false;
            if (bucket!=cache.end()) for (auto previous:bucket->second)
                if (previous->mapping==e.mapping && previous->islands==e.islands
                    && previous->cuts==e.cuts && previous->deferred==e.deferred) {
                    previous->in_snapshot=true;present=true;break;
                }
            if (!present && e.bytes<=budget) {e.in_snapshot=true;insert(std::move(e));}
        }
        return std::any_of(lru.begin(),lru.end(),[](const Entry& e){return !e.in_snapshot;});
    }
    FragmentRepair(const PySource& source,const PyTarget& product,size_t cache_bytes,size_t extension_bytes=0,
                   const std::string& growth_policy="incremental")
        :base(source.g),active_source(source.g),target(product.g),budget(cache_bytes),
         existence(existence_policy(growth_policy)) {
        base.reads=nullptr;active_source.reads=nullptr;
        if (existence && extension_bytes)
            throw std::invalid_argument("existence growth does not support a candidate extension cache");
        if (extension_bytes) extensions=std::make_unique<ExtensionCache>(extension_bytes);
    }
    bool uses_existence() const { return existence; }
    std::pair<GrowResult,long> compute(const std::vector<Pair>& cuts_in,int seed,const std::vector<int>& mapping,
                    double graph_floor,double iso_tol,int min_lock_size,long max_branches,
                    const std::vector<Pair>* islands_ptr,const std::vector<Pair>& prior_deferred,bool allow_mapped_seed, uint64_t random_seed=42) {
        const std::vector<Pair> empty_islands;
        const auto& islands=islands_ptr ? *islands_ptr : empty_islands;
        GrowResult result;
        long saved_extensions=0;
        {
            ++calls;
            auto cuts=normalized(cuts_in),deferred=normalized(prior_deferred);
            if (cuts!=active_cuts) {
                for (auto [a,b]:cuts)
                    if (a<0 || b<0 || a>=base.n || b>=base.n || !base.has_edge(a,b))
                        throw std::invalid_argument("cut is not an active edge of the repair source");
                active_source=base;active_cuts=cuts;
                for (auto [a,b]:cuts) {
                    active_source.edge[(size_t)a*base.n+b]=active_source.edge[(size_t)b*base.n+a]=0;
                    for (auto [x,y]:{Pair{a,b},Pair{b,a}}) {
                        auto& neighbors=active_source.adj[x];
                        neighbors.erase(std::remove(neighbors.begin(),neighbors.end(),y),neighbors.end());
                    }
                }
            }
            std::string key;
            append(key,random_seed);append(key,seed);append(key,graph_floor);append(key,iso_tol);
            append(key,min_lock_size);append(key,max_branches);append(key,allow_mapped_seed);
            std::vector<uint8_t> occupied(target.n,0);
            for (int p:mapping) if (p>=0) occupied[p]=1;
            key.append(reinterpret_cast<const char*>(occupied.data()),occupied.size());
            // Occupancy fixes the same pointwise target stabilizer. Relabeling
            // untouched locked R atoms changes color names, not the equivalence
            // relation between candidates. Touched mapping values are checked.
            auto bucket=cache.find(key);
            auto chosen=lru.end();
            if (bucket!=cache.end()) for (auto item=bucket->second.rbegin();item!=bucket->second.rend();++item)
                if (valid(**item,mapping,islands,cuts,deferred)) {
                    chosen=*item;
                    std::iter_swap(item,bucket->second.rbegin());
                    break;
                }
            if (chosen!=lru.end()) {
                ++hits;
                history_hits+=(chosen->mapping!=mapping || chosen->islands!=islands);
                topology_hits+=(chosen->cuts!=cuts);
                result=chosen->result;
                if (chosen->deferred!=deferred) {
                    ++boundary_rebases;
                    auto rebased=deferred;
                    rebased.insert(rebased.end(),chosen->dependencies.local_deferred.begin(),chosen->dependencies.local_deferred.end());
                    rebased=normalized(std::move(rebased));
                    for (auto& iso:result.isos) iso.deferred_edges=rebased;
                }
                saved_extensions=result.profile.extend_calls;
                reused_extensions+=saved_extensions;
                reused_certificates+=result.profile.certificate_calls;
                lru.splice(lru.end(),lru,chosen);
            } else {
                Entry entry(base.n);
                // Complete-result reuse needs dependency reads, not candidate
                // checkpoints. Both growth policies record these same inputs.
                SourceReads reads(base.n);
                entry.result=grow_island(active_source,target,seed,mapping,graph_floor,iso_tol,min_lock_size,
                    max_branches,islands_ptr,prior_deferred,allow_mapped_seed,nullptr,nullptr,-1,
                    &entry.dependencies,extensions.get(),random_seed,nullptr,existence,&reads);
                entry.key=std::move(key);entry.mapping=mapping;
                entry.islands=islands;entry.cuts=cuts;entry.deferred=deferred;
                entry.row_reads.resize(base.n);
                for (int r=0;r<base.n;++r) entry.row_reads[r]=(reads.rows[r]!=-2);
                for (int a=0;a<base.n;++a) for (int b=a+1;b<base.n;++b)
                    if (!entry.row_reads[a] && !entry.row_reads[b] && reads.pairs[(size_t)a*base.n+b]!=-2)
                        entry.pair_reads.push_back({a,b});
                entry.bytes=entry_bytes(entry);
                if (entry.bytes>budget) {
                    ++oversized;
                    result=std::move(entry.result);
                } else {
                    result=entry.result;
                    insert(std::move(entry));
                }
            }
            logical_extensions+=result.profile.extend_calls;
            logical_certificates+=result.profile.certificate_calls;
        }
        return {std::move(result),saved_extensions};
    }
    py::object grow(const std::vector<Pair>& cuts,int seed,const std::vector<int>& mapping,
                    double graph_floor,double iso_tol,int min_lock_size,long max_branches,
                    py::object islands_obj,const std::vector<Pair>& deferred,bool allow_mapped_seed, uint64_t random_seed=42) {
        std::vector<Pair> islands;
        const std::vector<Pair>* islands_ptr=nullptr;
        if (!islands_obj.is_none()) {islands=islands_obj.cast<std::vector<Pair>>();islands_ptr=&islands;}
        std::pair<GrowResult,long> computed;
        {
            py::gil_scoped_release release;
            computed=compute(cuts,seed,mapping,graph_floor,iso_tol,min_lock_size,max_branches,
                             islands_ptr,deferred,allow_mapped_seed,random_seed);
        }
        auto out=grow_result_dict(computed.first).cast<py::dict>();
        out["reused_extensions"]=computed.second;
        return out;
    }
    py::dict stats() const {
        py::dict out;
        out["calls"]=calls;out["full_hits"]=hits;out["prefix_hits"]=0;
        out["history_hits"]=history_hits;out["topology_hits"]=topology_hits;out["boundary_rebases"]=boundary_rebases;
        out["dependency_checks"]=checks;out["mapping_misses"]=mapping_misses;out["island_misses"]=island_misses;
        out["topology_misses"]=topology_misses;out["boundary_misses"]=boundary_misses;
        out["logical_extensions"]=logical_extensions;out["reused_extensions"]=reused_extensions;
        out["logical_certificates"]=logical_certificates;out["reused_certificates"]=reused_certificates;
        out["entries"]=lru.size();out["resident_bytes"]=resident;out["peak_resident_bytes"]=peak_resident;
        out["cache_budget_bytes"]=budget;out["evictions"]=evictions;out["oversized_results"]=oversized;
        out["revision"]=revision;
        if (extensions) out.attr("update")(extensions->stats());
        return out;
    }
};

void register_fragment_repair(py::module_& mod) {
    py::class_<FragmentRepair>(mod,"FragmentRepair")
        .def(py::init<const PySource&,const PyTarget&,size_t,size_t,const std::string&>(),
            py::arg("source"),py::arg("target"),py::arg("cache_bytes"),py::arg("extension_cache_bytes")=0,
            py::arg("growth_policy")="incremental",py::keep_alive<1,3>())
        .def_property_readonly("growth_policy",[](const FragmentRepair& cache) {
            return cache.uses_existence() ? "existence" : "incremental";
        })
        .def("snapshot_context",&FragmentRepair::snapshot_context)
        .def("dump",&FragmentRepair::dump)
        .def("load",&FragmentRepair::load)
        .def("grow",&FragmentRepair::grow,py::arg("cuts"),py::arg("seed"),py::arg("mapping"),
             py::arg("graph_floor"),py::arg("iso_tol"),py::arg("min_lock_size"),py::arg("max_branches"),
             py::arg("islands"),py::arg("prior_deferred_edges"),py::arg("allow_mapped_seed"),py::arg("random_seed")=42)
        .def("stats",&FragmentRepair::stats);
}
