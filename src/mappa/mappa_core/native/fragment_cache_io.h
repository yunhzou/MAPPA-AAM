// Trusted, versioned snapshots of completed fragments. No search or scheduling.
#include <type_traits>

namespace fragment_cache_io {
inline const std::string schema="mappa.fragment-cache/12";

struct Writer {
    std::string data;
    std::unordered_map<const PreparedSymmetry*,uint64_t> prepared;
    std::unordered_map<const Fragment*,uint64_t> fragments;
    std::unordered_map<const PlacementGenerators*,uint64_t> placements;
    void write(const std::shared_ptr<const PreparedSymmetry>& value);
    void write(const std::shared_ptr<const Fragment>& value);
    void write(const PlacementGroup& value);
    template<class... Values> void operator()(const Values&... values) { (write(values), ...); }
    template<class Value> void write(const Value& value) {
        if constexpr (std::is_arithmetic_v<Value>)
            data.append(reinterpret_cast<const char*>(&value), sizeof(value));
        else fields(*this, value);
    }
    void write(const std::string& value) {
        write(uint64_t(value.size())); data.append(value);
    }
    template<class Value> void write(const std::vector<Value>& values) {
        write(uint64_t(values.size()));
        if constexpr (std::is_arithmetic_v<Value>) {
            if (!values.empty()) data.append(reinterpret_cast<const char*>(values.data()), values.size()*sizeof(Value));
        } else for (const auto& value: values) write(value);
    }
    template<class A, class B> void write(const std::pair<A,B>& value) { (*this)(value.first, value.second); }
};

struct Reader {
    const std::string& data;
    size_t position=0;
    std::vector<std::shared_ptr<const PreparedSymmetry>> prepared;
    std::vector<std::shared_ptr<const Fragment>> fragments;
    std::vector<PlacementGroup> placements;
    void read(std::shared_ptr<const PreparedSymmetry>& value);
    void read(std::shared_ptr<const Fragment>& value);
    void read(PlacementGroup& value);
    explicit Reader(const std::string& bytes): data(bytes) {}
    void require(size_t size) const {
        if (size>data.size()-position) throw std::invalid_argument("Truncated fragment cache");
    }
    uint64_t size() {
        uint64_t count; read(count); require(count); return count;
    }
    template<class... Values> void operator()(Values&... values) { (read(values), ...); }
    template<class Value> void read(Value& value) {
        if constexpr (std::is_arithmetic_v<Value>) {
            require(sizeof(value)); std::memcpy(&value, data.data()+position, sizeof(value)); position+=sizeof(value);
        } else fields(*this, value);
    }
    void read(std::string& value) {
        auto count=size(); value.assign(data.data()+position, count); position+=count;
    }
    template<class Value> void read(std::vector<Value>& values) {
        auto count=size();
        if constexpr (std::is_arithmetic_v<Value>) {
            if (count>(data.size()-position)/sizeof(Value)) throw std::invalid_argument("Truncated fragment cache vector");
            values.resize(count);
            if (count) std::memcpy(values.data(), data.data()+position, count*sizeof(Value));
            position+=count*sizeof(Value);
        } else {
            values.resize(count);
            for (auto& value: values) read(value);
        }
    }
    template<class A, class B> void read(std::pair<A,B>& value) { (*this)(value.first, value.second); }
    void finish() const {
        if (position!=data.size()) throw std::invalid_argument("Trailing fragment cache data");
    }
};

template<class IO, class T> void self(IO& io,T& v) {
    io(v.atoms,v.anchors,v.elements,v.bond_threshold,v.weights,v.generators);
}
template<class IO> void fields(IO& io,PreparedSymmetry& v) { self(io,v); }
template<class IO> void fields(IO& io,const PreparedSymmetry& v) { self(io,v); }
inline void Writer::write(const std::shared_ptr<const PreparedSymmetry>& value) {
    auto [entry,added]=prepared.emplace(value.get(),prepared.size());
    write(entry->second);
    if (added) write(*value);
}
inline void Reader::read(std::shared_ptr<const PreparedSymmetry>& value) {
    uint64_t index;read(index);
    if (index>prepared.size()) throw std::invalid_argument("Invalid shared fragment self reference");
    if (index==prepared.size()) {
        auto item=std::make_shared<PreparedSymmetry>();read(*item);prepared.push_back(std::move(item));
    }
    value=prepared[index];
}
template<class IO, class T> void fragment(IO& io,T& v) {
    io(v.source,v.target,v.domains,v.witness,v.anchors,v.accepted_edges,v.target_edges,
       v.edge_tolerances,v.inclusive,v.target_bond_floor);
}
template<class IO> void fields(IO& io,Fragment& v) { fragment(io,v); }
template<class IO> void fields(IO& io,const Fragment& v) { fragment(io,v); }
inline void Writer::write(const std::shared_ptr<const Fragment>& value) {
    auto [entry,added]=fragments.emplace(value.get(),fragments.size());
    write(entry->second);
    if (added) write(*value);
}
inline void Reader::read(std::shared_ptr<const Fragment>& value) {
    uint64_t index;read(index);
    if (index>fragments.size()) throw std::invalid_argument("Invalid shared fragment reference");
    if (index==fragments.size()) {
        auto item=std::make_shared<Fragment>();read(*item);fragments.push_back(std::move(item));
    }
    value=fragments[index];
}
inline void Writer::write(const PlacementGroup& value) {
    auto [entry,added]=placements.emplace(value.get(),placements.size());
    write(entry->second);
    if (added) write(*value);
}
inline void Reader::read(PlacementGroup& value) {
    uint64_t index;read(index);
    if (index>placements.size()) throw std::invalid_argument("Invalid shared placement group reference");
    if (index==placements.size()) {
        auto item=std::make_shared<PlacementGenerators>();read(*item);placements.push_back(std::move(item));
    }
    value=placements[index];
}
template<class IO, class T> void profile(IO& io, T& v) {
    io(v.feasibility_calls,v.feasibility_nodes,v.feasibility_searches,
       v.completion_queries,v.completion_pruned,v.completion_seconds,
       v.completion_generated,v.completion_peak_candidates,
       v.completion_global_generators,v.completion_global_pruned,
       v.seed_targets,v.init_cands,v.initial_heap,v.heap_pops,v.stale_pops,
       v.fragment_skip_pops,v.extend_calls,v.max_cands_before,v.max_cands_after,
       v.max_fragment_size,v.max_heap_len,v.commits,v.deferred,v.merge_calls,
       v.free_extend_calls,v.final_cands,v.final_fragment_size,v.branches,v.result);
}
template<class IO> void fields(IO& io, GrowProfile& v) { profile(io,v); }
template<class IO> void fields(IO& io, const GrowProfile& v) { profile(io,v); }
template<class IO, class T> void result(IO& io, T& v) {
    io(v.fragments,v.placement_generators,v.deferred_edges,v.profile,v.capped,v.cap.count,v.cap.limit);
}
template<class IO> void fields(IO& io, GrowResult& v) { result(io,v); }
template<class IO> void fields(IO& io, const GrowResult& v) { result(io,v); }
template<class IO, class T> void dependencies(IO& io, T& v) { io(v.mapping_reads,v.boundary_atoms,v.island_reads,v.local_deferred); }
template<class IO> void fields(IO& io, GrowthDependencies& v) { dependencies(io,v); }
template<class IO> void fields(IO& io, const GrowthDependencies& v) { dependencies(io,v); }
}
