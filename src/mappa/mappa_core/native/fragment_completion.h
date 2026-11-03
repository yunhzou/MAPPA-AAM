// Saturation discovers target atom sets and retains the full conditional relation.
std::shared_ptr<const PreparedSelf> prepare_self(const Graph& graph,
        std::vector<int> atoms,std::vector<int> anchors,double threshold,double tolerance) {
    auto result=std::make_shared<PreparedSelf>();
    result->atoms=std::move(atoms);result->anchors=std::move(anchors);
    result->bond_threshold=threshold;result->tolerance=tolerance;
    const size_t n=result->atoms.size();
    result->weights.reserve(n*n);
    std::vector<Pair> edges;
    std::map<std::pair<std::string,int>,size_t> colors;
    std::vector<std::vector<int>> cells;
    for (size_t i=0;i<n;++i) {
        int atom=result->atoms[i];
        result->elements.push_back(graph.elem[atom]);
        auto key=std::make_pair(graph.elem[atom],contains(result->anchors,atom) ? atom : -1);
        auto [entry,added]=colors.emplace(key,cells.size());
        if (added) cells.emplace_back();
        cells[entry->second].push_back((int)i);
        for (size_t j=0;j<n;++j) {
            double weight=graph.wbo(atom,result->atoms[j]);
            result->weights.push_back(weight);
            if (i<j && weight>threshold) edges.push_back({(int)i,(int)j});
        }
    }
    if (n) result->generators=autgrp_native::AutGraph((int)n,edges).generators(cells);
    return result;
}

void complete_fragment(const Source& source,const Target& target,const std::vector<int>& locked,
        FragmentGrowth& growth,double threshold,double tolerance,double self_tolerance,
        int min_lock_size,long max_branches,GrowResult& result,GrowthObservation* observation) {
    auto started=std::chrono::steady_clock::now();
    auto& profile=result.profile;
    FragmentQuery query(source,target,growth.source_atoms,growth.accepted_edges,locked,tolerance);
    const PreparedSelf* global=nullptr;
    if (query.needs_global_symmetry()) {
        if (!target.global_self || target.global_self->bond_threshold!=threshold ||
            target.global_self->tolerance!=self_tolerance) {
            std::vector<int> atoms(target.n);std::iota(atoms.begin(),atoms.end(),0);
            target.global_self=prepare_self(target,std::move(atoms),{},threshold,self_tolerance);
            target.global_group_allowed=group_within_tolerance(*target.global_self,target.global_self->generators);
        }
        global=target.global_self.get();
    }
    auto witnesses=query.images(max_branches+1,&growth.witness.assignment(),global,target.global_group_allowed);
    if (witnesses.empty()) throw std::runtime_error("completed fragment lost its proven matching witness");
    profile.completion_queries=query.nodes;profile.completion_pruned=query.rejected;
    profile.completion_global_generators=query.global_generators;
    profile.completion_global_images=query.global_images;
    profile.completion_global_pruned=query.global_pruned;
    profile.completion_generated=(long)witnesses.size();
    profile.completion_peak_candidates=(long)witnesses.size();
    profile.final_cands=(long)witnesses.size();profile.final_fragment_size=growth.atom_count;
    if (observation) observation->record(witnesses,growth.source_atoms,growth.deferred_edges);
    if (growth.atom_count<min_lock_size) profile.result="too_small";
    else if ((long)witnesses.size()>max_branches) {
        profile.result="subtree_branch_cap";profile.branches=(long)witnesses.size();
        result.capped=true;result.cap={(long)witnesses.size(),max_branches};
    } else {
        std::vector<int> anchors;
        for (int r:query.atoms) if (locked[r]>=0) anchors.push_back(r);
        auto prepared_source=prepare_self(source,query.atoms,anchors,threshold,self_tolerance);
        auto bonds=growth.accepted_edges;
        std::sort(bonds.begin(),bonds.end());
        std::vector<int> source_positions(source.n,-1);
        for (size_t i=0;i<query.atoms.size();++i) source_positions[query.atoms[i]]=(int)i;
        for (const auto& witness:witnesses) {
            std::vector<int> target_atoms,target_anchors;
            for (int r:query.atoms) target_atoms.push_back(witness[r]);
            std::sort(target_atoms.begin(),target_atoms.end());
            for (int r:anchors) target_anchors.push_back(locked[r]);
            std::sort(target_anchors.begin(),target_anchors.end());
            Fragment fragment;
            fragment.source=prepared_source;
            fragment.target=prepare_self(target,std::move(target_atoms),std::move(target_anchors),threshold,self_tolerance);
            fragment.inclusive=source.tolerance_inclusive;fragment.target_bond_floor=target.bond_cut;
            std::vector<int> target_positions(target.n,-1);
            for (size_t i=0;i<fragment.target->atoms.size();++i)
                target_positions[fragment.target->atoms[i]]=(int)i;
            auto domains=query.domains;
            for (int r:query.atoms) {
                auto& values=domains[r];
                values.erase(std::remove_if(values.begin(),values.end(),
                    [&](int p){return target_positions[p]<0;}),values.end());
            }
            if (!query.reduce(domains))
                throw std::runtime_error("fixed fragment image lost its proven matching witness");
            for (int r:query.atoms) {
                std::vector<int> domain;
                for (int p:domains[r]) domain.push_back(target_positions[p]);
                fragment.domains.push_back(std::move(domain));
                fragment.witness.push_back(target_positions[witness[r]]);
                if (locked[r]>=0) fragment.anchors.push_back({source_positions[r],target_positions[locked[r]]});
            }
            for (auto [a,b]:bonds) {
                fragment.accepted_edges.push_back({source_positions[a],source_positions[b]});
                fragment.edge_tolerances.push_back(source.tolerance(a,b,tolerance));
            }
            const auto& values=fragment.target->atoms;
            for (size_t i=0;i<values.size();++i) for (size_t j=i+1;j<values.size();++j)
                if (target.has_edge(values[i],values[j])) fragment.target_edges.push_back({(int)i,(int)j});
            result.fragments.push_back(std::move(fragment));
        }
        result.deferred_edges=std::move(growth.deferred_edges);
        profile.branches=(long)result.fragments.size();
        profile.result=result.fragments.size()==1 ? "success" : "branched";
    }
    profile.completion_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
