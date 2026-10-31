// Complete a proven fragment with the established compressed matching rules.
struct AcceptedGrowthStep {
    int anchor,atom;
    double weight;
    std::vector<Pair> deferred;
};

std::vector<Cand> complete_fragment(
        const Source& source,const Target& target,int seed,const std::vector<int>& locked,
        const std::vector<uint8_t>& occupied,const std::vector<Pair>* islands,
        const std::vector<int>& island_of,const std::vector<AcceptedGrowthStep>& steps,
        const std::vector<uint8_t>& fragment,const std::vector<Pair>& bonds,
        double graph_floor,double tolerance,const Canonicalizer& canonicalizer,
        long& certificate_calls,GrowProfile& profile) {
    auto started=std::chrono::steady_clock::now();
    std::unique_ptr<FragmentQuery> query;
    auto candidates=seed_candidates(source,target,seed,locked);
    std::vector<uint8_t> prefix(source.n,0);prefix[seed]=1;
    auto retain_completable=[&]() {
        // Equivalent children have the same feasibility. A sole candidate
        // contains the proven match; full candidates need no further check.
        if (candidates.size()<=1 || prefix==fragment) return;
        if (!query) query=std::make_unique<FragmentQuery>(source,target,fragment,bonds,locked,tolerance);
        candidates.erase(std::remove_if(candidates.begin(),candidates.end(),
            [&](const Cand& candidate){return query->find(&candidate).empty();}),candidates.end());
    };
    retain_completable();
    auto context=fragment_context(source,target,locked,occupied,islands,island_of,tolerance);
    for (const auto& step:steps) {
        prepare_growth_step(context,prefix,step.anchor,step.atom,step.weight,step.deferred,graph_floor);
        candidates=extend_sym_cands(candidates,context,canonicalizer,certificate_calls);
        prefix.swap(context.sig_fragment);
        retain_completable();
    }
    if (candidates.empty()) throw std::runtime_error("completed fragment lost its proven matching witness");
    profile.completion_steps=(long)steps.size();
    profile.completion_queries=query ? query->calls : 0;
    profile.completion_pruned=query ? query->rejected : 0;
    profile.completion_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    return candidates;
}
