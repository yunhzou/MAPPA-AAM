// Trusted, versioned snapshots of completed fragments. No search or scheduling.
#include <type_traits>

namespace fragment_cache_io {
struct Writer {
    std::string data;
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

template<class IO, class T> void block(IO& io, T& v) { io(v.r,v.p,v.extendable); }
template<class IO> void fields(IO& io, Block& v) { block(io,v); }
template<class IO> void fields(IO& io, const Block& v) { block(io,v); }
template<class IO, class T> void candidate(IO& io, T& v) {
    // extension_key is a lazily reconstructed optimization, not matching data.
    io(v.img,v.mapped,v.blocks,v.exact_fixed,v.mult,v.automorph);
}
template<class IO> void fields(IO& io, Cand& v) { candidate(io,v); }
template<class IO> void fields(IO& io, const Cand& v) { candidate(io,v); }
template<class IO, class T> void iso(IO& io, T& v) { io(v.mapping,v.deferred_edges,v.fragment,v.preserved_bonds,v.cand); }
template<class IO> void fields(IO& io, IsoOut& v) { iso(io,v); }
template<class IO> void fields(IO& io, const IsoOut& v) { iso(io,v); }
template<class IO, class T> void profile(IO& io, T& v) {
    io(v.growth_policy,v.feasibility_calls,v.feasibility_nodes,v.feasibility_searches,
       v.completion_steps,v.completion_queries,v.completion_pruned,v.completion_seconds,
       v.seed_targets,v.seed_groups,v.init_cands,v.initial_heap,v.heap_pops,v.stale_pops,
       v.fragment_skip_pops,v.extend_calls,v.max_cands_before,v.max_cands_after,
       v.max_fragment_size,v.max_heap_len,v.commits,v.deferred,v.merge_calls,
       v.free_extend_calls,v.certificate_calls,v.final_cands,v.final_fragment_size,v.branches,v.result);
}
template<class IO> void fields(IO& io, GrowProfile& v) { profile(io,v); }
template<class IO> void fields(IO& io, const GrowProfile& v) { profile(io,v); }
template<class IO, class T> void result(IO& io, T& v) { io(v.isos,v.profile,v.capped,v.cap.count,v.cap.limit); }
template<class IO> void fields(IO& io, GrowResult& v) { result(io,v); }
template<class IO> void fields(IO& io, const GrowResult& v) { result(io,v); }
template<class IO, class T> void dependencies(IO& io, T& v) { io(v.mapping_reads,v.boundary_atoms,v.island_reads,v.local_deferred); }
template<class IO> void fields(IO& io, GrowthDependencies& v) { dependencies(io,v); }
template<class IO> void fields(IO& io, const GrowthDependencies& v) { dependencies(io,v); }
}
