// Python transport for the native fragment scheduler. No search policy lives here.
// Export shared immutable blocks directly, never a full temporary record per state.
class SharedStatePairs {
    std::unordered_map<uint64_t,py::tuple> rows;
    std::unordered_map<std::string,py::tuple> blocks;
    std::unordered_map<std::string,py::object> sequences;
    py::object sequence_type=py::module_::import("mappa.mappa_core.core.search.pairs").attr("SharedPairs");
public:
    py::tuple pairs(const std::vector<Pair>& values) {
        std::string key;
        key.reserve(values.size()*2*sizeof(int));
        for (auto [a,b]:values) {
            key.append(reinterpret_cast<const char*>(&a),sizeof(a));
            key.append(reinterpret_cast<const char*>(&b),sizeof(b));
        }
        auto found=blocks.find(key);
        if (found!=blocks.end()) return found->second;
        py::tuple result(values.size());
        for (size_t i=0;i<values.size();++i) {
            auto [a,b]=values[i];
            uint64_t row_key=(uint64_t(uint32_t(a))<<32)|uint32_t(b);
            auto row=rows.find(row_key);
            if (row==rows.end()) row=rows.emplace(row_key,py::make_tuple(a,b)).first;
            result[i]=row->second;
        }
        blocks.emplace(std::move(key),result);
        return result;
    }
    py::object state(const std::vector<int>& values,const std::vector<int>& source_ids,
                     const std::vector<int>* target_ids=nullptr) {
        constexpr size_t width=32;
        py::list parts;
        std::string key;
        size_t count=0;
        for (size_t start=0;start<values.size();start+=width) {
            std::vector<Pair> block;
            for (size_t r=start;r<std::min(start+width,values.size());++r) if (values[r]>=0)
                block.push_back({source_ids[r],target_ids ? (*target_ids)[values[r]] : values[r]});
            count+=block.size();
            if (!block.empty()) {
                auto shared=pairs(block);
                parts.append(shared);
                // Blocks already have unique identities in this export. Reuse
                // whole sequences without hashing or copying every atom again.
                auto identity=shared.ptr();
                key.append(reinterpret_cast<const char*>(&identity),sizeof(identity));
            }
        }
        auto found=sequences.find(key);
        if (found!=sequences.end()) return found->second;
        auto result=sequence_type(py::tuple(parts),count);
        sequences.emplace(std::move(key),result);
        return result;
    }
};

py::tuple scheduled_graph_result(const ScheduledGraph& graph,const std::vector<int>& r_ids,
                                 const std::vector<int>& p_ids,const py::object& context) {
    SharedStatePairs shared;
    auto types=py::module_::import("mappa.mappa_core.core.search.graph");
    auto state_type=types.attr("SearchState"), transition_type=types.attr("FragmentTransition");
    auto stop_type=types.attr("SearchStop");
    auto r_pairs=[&](const std::vector<Pair>& pairs) {
        std::vector<Pair> out;
        for (auto [a,b]:pairs) out.push_back({r_ids[a],r_ids[b]});
        return shared.pairs(out);
    };
    FragmentExport exported;
    py::dict self_pool;
    auto fragment_type=py::module_::import("mappa.mappa_core.core.matching.representation.records").attr("fragment_from_native");
    std::unordered_map<const Fragment*,py::object> fragments;
    py::tuple states(graph.states.size()),transitions(graph.transitions.size()),stops(graph.stops.size());
    for (size_t i=0;i<graph.states.size();++i) {
        const auto& state=graph.states[i];
        states[i]=state_type(i,0,shared.state(state.mapping,r_ids,&p_ids),
                            shared.state(state.labels,r_ids),r_pairs(state.deferred));
    }
    for (size_t i=0;i<graph.transitions.size();++i) {
        const auto& edge=graph.transitions[i];
        py::object fragment=py::none();
        if (edge.fragment) {
            auto found=fragments.find(edge.fragment.get());
            if (found==fragments.end()) found=fragments.emplace(edge.fragment.get(),
                fragment_type(exported.relation(*edge.fragment),r_ids,p_ids,self_pool)).first;
            fragment=found->second;
        }
        transitions[i]=transition_type(i,edge.source,edge.target,
            edge.seed<0 ? py::object(py::none()) : py::object(py::int_(r_ids[edge.seed])),
            py::make_tuple(edge.pass,edge.position),fragment,edge.island,
            exported.placement(edge.placement_generators,p_ids));
    }
    for (size_t i=0;i<graph.stops.size();++i) {
        const auto& stop=graph.stops[i];
        stops[i]=stop_type(stop.state<0 ? py::object(py::none()) : py::object(py::int_(stop.state)),stop.reason,
            stop.seed<0 ? py::object(py::none()) : py::object(py::int_(r_ids[stop.seed])),
            stop.stage,stop.count,stop.limit,py::make_tuple(stop.pass,stop.position),
            stop.state<0 ? py::object(py::int_(0)) : py::object(py::none()),stop.occurrences);
    }
    auto result=types.attr("AAMSearchGraph")(py::make_tuple(context),py::make_tuple(0),states,transitions,stops);
    return py::make_tuple(result,graph.growth_calls,graph.max_growth_candidates);
}

void register_fragment_scheduler(py::module_& mod) {
    mod.def("search_fragments",[](const PySource& r,const PyTarget& p,const std::vector<int>& order,
                                  double floor,double tol,long cap,const std::vector<Pair>& anchors,
                                  const std::vector<int>& core,bool stop_at_core,
                                  const std::vector<int>& r_ids,const std::vector<int>& p_ids,const py::object& context,
                                  FragmentRepair* repair,const std::vector<Pair>& cuts,uint64_t random_seed) {
        ScheduledGraph graph;
        {
            py::gil_scoped_release release;
            graph=schedule_fragments(r.g,p.g,order,floor,tol,cap,anchors,core,stop_at_core,repair,cuts,random_seed);
        }
        return scheduled_graph_result(graph,r_ids,p_ids,context);
    },py::arg("source"),py::arg("target"),py::arg("seed_order"),py::arg("graph_floor"),
       py::arg("iso_tolerance"),py::arg("branch_cap"),py::arg("anchors"),py::arg("core_atoms"),
       py::arg("stop_at_core"),py::arg("source_ids"),py::arg("target_ids"),py::arg("context"),
       py::arg("repair")=nullptr,py::arg("cuts")=std::vector<Pair>{},py::arg("random_seed")=42);
}
