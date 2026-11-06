// Python transport for the unique-fragment pool. Every call grows afresh.
void register_fragment_pool(py::module_& mod) {
    py::class_<FragmentPool>(mod,"FragmentPool")
        .def(py::init([](const PySource& source,const PyTarget& target,size_t bytes) {
            return std::make_unique<FragmentPool>(source.g,target.g,bytes);
        }),py::arg("source"),py::arg("target"),py::arg("cache_bytes"),py::keep_alive<1,3>())
        .def("snapshot_context",&FragmentPool::snapshot_context)
        .def("dump",&FragmentPool::dump)
        .def("load",&FragmentPool::load)
        .def("stats",&FragmentPool::stats)
        .def("grow",[](FragmentPool& pool,const std::vector<Pair>& cuts,int seed,
                const std::vector<int>& mapping,double floor,double tolerance,int minimum,long cap,
                py::object islands_obj,const std::vector<Pair>& deferred,bool anchored,uint64_t random_seed,
                const std::vector<int>& fixed) {
            std::vector<Pair> islands;
            const bool has_islands=!islands_obj.is_none();
            if (has_islands) islands=islands_obj.cast<std::vector<Pair>>();
            GrowResult result;
            {
                py::gil_scoped_release release;
                pool.cut(cuts);
                result=grow_island(pool.active_source,pool.target,seed,mapping,floor,tolerance,
                    minimum,cap,has_islands ? &islands : nullptr,deferred,anchored,
                    random_seed,nullptr,&pool,fixed);
            }
            return grow_result_dict(result);
        },py::arg("cuts"),py::arg("seed"),py::arg("mapping"),py::arg("graph_floor"),py::arg("iso_tol"),
          py::arg("min_lock_size"),py::arg("max_branches"),py::arg("islands"),
          py::arg("prior_deferred_edges"),py::arg("allow_mapped_seed"),py::arg("random_seed")=42,
          py::arg("placement_fixed_targets")=std::vector<int>{})
        .def("exists",[](FragmentPool& pool,const std::vector<int>& atoms,const std::vector<Pair>& edges,
                const std::vector<int>& locked,double tolerance) {
            if (locked.size()!=size_t(pool.base.n)) throw std::invalid_argument("Locked mapping size differs");
            for (int image:locked)
                if (image < -1 || image >= pool.target.n)
                    throw std::invalid_argument("Locked image lies outside the target");
            std::vector<uint8_t> selected(pool.base.n,0);
            for (int atom:atoms) selected.at(atom)=1;
            for (auto [a,b]:edges) if (!selected.at(a) || !selected.at(b))
                throw std::invalid_argument("Growth edge outside fragment");
            py::gil_scoped_release release;
            if (!pool.find(selected,edges,locked,tolerance).empty()) return true;
            FragmentQuery query(pool.base,pool.target,selected,edges,locked,tolerance);
            return !query.find().empty();
        },py::arg("atoms"),py::arg("edges"),py::arg("locked"),py::arg("tolerance"));
}
