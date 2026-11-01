// Reaction branch search and decision DAG construction, without worker scheduling
// or Python transport. The Python search is the full-graph equivalence oracle.
struct ScheduledBranch {
    std::vector<int> mapping, labels;
    std::vector<Pair> islands, deferred;
    int node=0, mapped=0;
    std::string key;
    size_t key_hash=0;
    void snapshot_key() {
        key.assign(reinterpret_cast<const char*>(mapping.data()),mapping.size()*sizeof(int));
        key.append(reinterpret_cast<const char*>(labels.data()),labels.size()*sizeof(int));
        for (auto [a,b]:deferred) {
            key.append(reinterpret_cast<const char*>(&a),sizeof(a));
            key.append(reinterpret_cast<const char*>(&b),sizeof(b));
        }
        key_hash=std::hash<std::string>()(key);
    }
};

struct ScheduledKeyHash {
    size_t operator()(const ScheduledBranch* b) const { return b->key_hash; }
};
struct ScheduledKeyEqual {
    bool operator()(const ScheduledBranch* a,const ScheduledBranch* b) const {
        return a==b || (a->key_hash==b->key_hash && a->key==b->key);
    }
};

struct ScheduledState {
    std::vector<int> mapping, labels;
    std::vector<Pair> deferred;
};
struct ScheduledTransition {
    int source,target,seed,pass,position,island=0;
    std::shared_ptr<IsoOut> match;
    std::vector<Pair> preserved;
};
struct ScheduledStop {
    int state,seed,pass,position;
    std::string reason,stage;
    long count=0,limit=0;
};
struct ScheduledGraph {
    std::vector<ScheduledState> states;
    std::vector<ScheduledTransition> transitions;
    std::vector<ScheduledStop> stops;
    int seed=-1,pass=0,position=0;
    long growth_calls=0;
    long max_growth_candidates=0;
    int state(const ScheduledBranch& b) {
        states.push_back({b.mapping,b.labels,b.deferred});
        return (int)states.size()-1;
    }
    void join(ScheduledBranch& kept,const ScheduledBranch& other) {
        if (kept.node==other.node) return;
        int node=state(kept);
        transitions.push_back({kept.node,node,seed,pass,position,0,nullptr,{}});
        transitions.push_back({other.node,node,seed,pass,position,0,nullptr,{}});
        kept.node=node;
    }
    void stop(const ScheduledBranch& b,std::string reason,std::string stage="",long count=0,long limit=0) {
        stops.push_back({b.node,seed,pass,position,std::move(reason),std::move(stage),count,limit});
    }
};

using BranchPtr=std::shared_ptr<ScheduledBranch>;
using ScheduledFrontier=std::unordered_map<const ScheduledBranch*,BranchPtr,ScheduledKeyHash,ScheduledKeyEqual>;

BranchPtr commit_scheduled_fragment(const BranchPtr& parent,IsoOut iso,
                                    const Source& source,ScheduledGraph& graph) {
    auto child=std::make_shared<ScheduledBranch>(*parent);
    for (auto [r,p]:iso.mapping) if (child->mapping[r]<0) {
        child->mapping[r]=p;++child->mapped;
    }
    child->deferred.insert(child->deferred.end(),iso.deferred_edges.begin(),iso.deferred_edges.end());
    std::sort(child->deferred.begin(),child->deferred.end());
    child->deferred.erase(std::unique(child->deferred.begin(),child->deferred.end()),child->deferred.end());
    std::vector<uint8_t> in_merge(source.n,0);
    std::unordered_set<int> touched;
    for (auto [r,p]:iso.mapping) {
        in_merge[r]=1;
        if (parent->labels[r]>=0) touched.insert(parent->labels[r]);
    }
    std::map<int,std::vector<int>> groups;
    for (int r=0;r<source.n;++r) if (parent->labels[r]>=0) {
        if (touched.count(parent->labels[r])) in_merge[r]=1;
        else groups[parent->labels[r]].push_back(r);
    }
    std::vector<std::vector<int>> components;
    for (auto& item:groups) components.push_back(std::move(item.second));
    std::vector<int> merged;
    for (int r=0;r<source.n;++r) if (in_merge[r]) merged.push_back(r);
    if (!merged.empty()) components.push_back(std::move(merged));
    std::sort(components.begin(),components.end());
    child->islands.clear();child->labels.assign(source.n,-1);
    for (size_t i=0;i<components.size();++i) for (int r:components[i]) {
        child->labels[r]=(int)i+1;child->islands.push_back({r,(int)i+1});
    }
    child->snapshot_key();
    if (child->key==parent->key) return parent;
    auto preserved = iso.preserved_bonds;
    int island=child->labels[iso.mapping.front().first];
    child->node=graph.state(*child);
    graph.transitions.push_back({parent->node,child->node,graph.seed,graph.pass,graph.position,island,
        std::make_shared<IsoOut>(std::move(iso)),std::move(preserved)});
    return child;
}

ScheduledGraph schedule_fragments(const Source& source,const Target& target,
                                  const std::vector<int>& order,double graph_floor,double iso_tol,
                                  long max_branches,const std::vector<Pair>& anchors,
                                  const std::vector<int>& core,bool stop_at_core,
                                  FragmentRepair* repair,const std::vector<Pair>& cuts,uint64_t random_seed,
                                  bool existence=false) {
    ScheduledGraph graph;
    auto root=std::make_shared<ScheduledBranch>();
    root->mapping.assign(source.n,-1);root->labels.assign(source.n,-1);
    std::vector<uint8_t> anchored(source.n,0);
    bool compatible=true;
    for (auto [r,p]:anchors) if (source.ecode[r]!=target.ecode[p]) compatible=false;
    if (compatible) {
        for (auto [r,p]:anchors) {
            root->mapping[r]=p;root->labels[r]=++root->mapped;
            root->islands.push_back({r,root->mapped});anchored[r]=1;
        }
    }
    root->snapshot_key();root->node=graph.state(*root);
    if (!compatible) {graph.stop(*root,"incompatible_anchors");return graph;}
    auto core_complete=[&](const ScheduledBranch& b) {
        if (core.empty()) return false;
        for (int r:core) if (b.mapping[r]<0) return false;
        return true;
    };
    auto mapped_anchor_can_seed=[&](const ScheduledBranch& b,int seed) {
        for (int nb:source.adj[seed])
            if (source.wbo(seed,nb)>=graph_floor && b.labels[nb]!=b.labels[seed]) return true;
        return false;
    };
    std::vector<BranchPtr> branches{root};
    bool progressed=true;
    int pass_no=0;
    while (progressed) {
        progressed=false;++pass_no;
        for (size_t position=0;position<order.size();++position) {
            int seed=order[position];graph.seed=seed;graph.position=(int)position;graph.pass=pass_no;
            bool any_work=std::any_of(branches.begin(),branches.end(),[&](const auto& b) {
                if (stop_at_core && core_complete(*b)) return false;
                return b->mapping[seed]<0 || (anchored[seed] && mapped_anchor_can_seed(*b,seed));
            });
            // The admitted frontier is already unique and sorted. Carrying it
            // unchanged cannot create joins or cap stops; only the trace step
            // and requested-core stop need advancing.
            if (!any_work && (long)branches.size()<=max_branches) {
                if (stop_at_core && !core.empty() && !branches.empty() &&
                    std::all_of(branches.begin(),branches.end(),[&](const auto& b){return core_complete(*b);})) {
                    for (const auto& b:branches) graph.stop(*b,"objective_met");
                    return graph;
                }
                continue;
            }
            std::vector<BranchPtr> next;
            ScheduledFrontier seen;
            seen.reserve(branches.size());next.reserve(branches.size());
            auto admit=[&](const std::vector<BranchPtr>& subtree,bool made_progress=false) {
                std::vector<BranchPtr> additions;
                ScheduledFrontier local;
                for (const auto& b:subtree) {
                    auto prior=seen.find(b.get());
                    if (prior!=seen.end()) {graph.join(*prior->second,*b);continue;}
                    auto sibling=local.find(b.get());
                    if (sibling!=local.end()) {graph.join(*sibling->second,*b);continue;}
                    local.emplace(b.get(),b);additions.push_back(b);
                }
                if ((long)(next.size()+additions.size())>max_branches) {
                    for (const auto& b:additions)
                        graph.stop(*b,"capped","combined_live_leaves",next.size()+additions.size(),max_branches);
                    return;
                }
                for (const auto& b:additions) {seen.emplace(b.get(),b);next.push_back(b);}
                if (!additions.empty() && made_progress) progressed=true;
            };
            auto carry=[&](const BranchPtr& b) {
                auto prior=seen.find(b.get());
                if (prior!=seen.end()) {graph.join(*prior->second,*b);return;}
                if ((long)next.size()+1>max_branches) {
                    graph.stop(*b,"capped","combined_live_leaves",next.size()+1,max_branches);return;
                }
                seen.emplace(b.get(),b);next.push_back(b);
            };
            for (const auto& branch:branches) {
                if (stop_at_core && core_complete(*branch)) {carry(branch);continue;}
                bool anchor=branch->mapping[seed]>=0 && anchored[seed];
                if (branch->mapping[seed]>=0 && !anchor) {carry(branch);continue;}
                if (anchor && !mapped_anchor_can_seed(*branch,seed)) {carry(branch);continue;}
                ++graph.growth_calls;
                auto result=repair
                    ? repair->compute(cuts,seed,branch->mapping,graph_floor,iso_tol,1,max_branches,
                                      &branch->islands,branch->deferred,anchor,random_seed).first
                    : grow_island(source,target,seed,branch->mapping,graph_floor,iso_tol,1,max_branches,
                                  &branch->islands,branch->deferred,anchor,nullptr,nullptr,-1,nullptr,nullptr,random_seed,nullptr,existence);
                graph.max_growth_candidates=std::max(graph.max_growth_candidates,result.profile.max_cands_before);
                if (result.capped) {
                    graph.stop(*branch,"capped","fragment_growth",result.cap.count,result.cap.limit);
                    continue;
                }
                if (result.isos.empty()) {carry(branch);continue;}
                std::vector<BranchPtr> subtree;
                bool changed=false;
                for (auto& iso:result.isos) {
                    auto child=commit_scheduled_fragment(branch,std::move(iso),source,graph);
                    changed|=(child->key!=branch->key);subtree.push_back(std::move(child));
                }
                admit(subtree,changed);
            }
            std::stable_sort(next.begin(),next.end(),[](const auto& a,const auto& b){return a->mapped>b->mapped;});
            branches=std::move(next);
            if (stop_at_core && !core.empty() && !branches.empty() &&
                std::all_of(branches.begin(),branches.end(),[&](const auto& b){return core_complete(*b);})) {
                for (const auto& b:branches) graph.stop(*b,"objective_met");
                return graph;
            }
        }
    }
    for (const auto& b:branches) graph.stop(*b,b->mapped==source.n ? "objective_met" : "stalled");
    return graph;
}
