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

using BranchPtr=std::shared_ptr<ScheduledBranch>;

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
    std::shared_ptr<const Fragment> fragment;
    PlacementGroup placement_generators;
};
struct ScheduledStop {
    int state,seed,pass,position;
    std::string reason,stage;
    long count=0,limit=0;
    long occurrences=1;
};
struct ScheduledGraph {
    std::vector<ScheduledState> states;
    std::vector<ScheduledTransition> transitions;
    std::vector<ScheduledStop> stops;
    int seed=-1,pass=0,position=0;
    long growth_calls=0;
    long max_growth_candidates=0;
    std::map<std::tuple<int,int,int,std::string,std::string,long,long>,size_t> cap_summaries;
    bool needs_pruning=false;
    int state(const ScheduledBranch& b) {
        states.push_back({b.mapping,b.labels,b.deferred});
        return (int)states.size()-1;
    }
    void join(ScheduledBranch& kept,const ScheduledBranch& other) {
        if (kept.node==other.node) return;
        int node=state(kept);
        transitions.push_back({kept.node,node,seed,pass,position,0,nullptr});
        transitions.push_back({other.node,node,seed,pass,position,0,nullptr});
        kept.node=node;
    }
    void stop(const ScheduledBranch& b,std::string reason,std::string stage="",long count=0,long limit=0) {
        int state=b.node;
        if (reason=="capped") {
            needs_pruning=true;
            auto key=std::make_tuple(seed,pass,position,reason,stage,count,limit);
            auto previous=cap_summaries.find(key);
            if (previous!=cap_summaries.end()) {
                ++stops[previous->second].occurrences;
                return;
            }
            cap_summaries.emplace(std::move(key),stops.size());
            state=-1;
        }
        stops.push_back({state,seed,pass,position,std::move(reason),std::move(stage),count,limit});
    }
    bool cap_frontier(std::vector<BranchPtr>& frontier,long limit) {
        if ((long)frontier.size()<=limit) return false;
        std::vector<std::vector<size_t>> outgoing(states.size());
        for (size_t i=0;i<transitions.size();++i) outgoing[transitions[i].source].push_back(i);
        std::vector<int> leaf(states.size(),-1);
        for (size_t i=0;i<frontier.size();++i) leaf[frontier[i]->node]=(int)i;
        using Leaves=std::shared_ptr<const std::vector<int>>;
        auto empty=std::make_shared<const std::vector<int>>();
        std::vector<Leaves> descendants(states.size(),empty);
        std::vector<uint8_t> removed(transitions.size(),0);
        bool changed=false;
        // States are topological. Joined histories share one continuation and
        // therefore one leaf identity, rather than multiplying ancestry paths.
        for (size_t node=states.size();node-->0;) {
            std::vector<size_t> children;
            for (size_t edge:outgoing[node])
                if (!descendants[transitions[edge].target]->empty()) children.push_back(edge);
            if (leaf[node]<0 && children.size()==1) {
                descendants[node]=descendants[transitions[children[0]].target];
                continue;
            }
            std::vector<int> retained;
            if (leaf[node]>=0) retained.push_back(leaf[node]);
            for (size_t edge:children) {
                const auto& values=*descendants[transitions[edge].target];
                retained.insert(retained.end(),values.begin(),values.end());
            }
            if (retained.empty()) continue;
            std::sort(retained.begin(),retained.end());
            retained.erase(std::unique(retained.begin(),retained.end()),retained.end());
            if ((long)retained.size()>limit) {
                std::unordered_map<int,size_t> references;
                references.reserve(retained.size());
                if (leaf[node]>=0) ++references[leaf[node]];
                for (size_t edge:children)
                    for (int value:*descendants[transitions[edge].target]) ++references[value];
                std::sort(children.begin(),children.end(),[&](size_t a,size_t b) {
                    size_t left=descendants[transitions[a].target]->size();
                    size_t right=descendants[transitions[b].target]->size();
                    return left!=right ? left>right : a>b;
                });
                long count=(long)retained.size();
                ScheduledBranch stopped;stopped.node=(int)node;
                for (size_t edge:children) {
                    if (count<=limit) break;
                    stop(stopped,"capped","combined_live_leaves",count,limit);
                    removed[edge]=1;changed=true;
                    for (int value:*descendants[transitions[edge].target])
                        if (--references[value]==0) --count;
                }
                retained.erase(std::remove_if(retained.begin(),retained.end(),
                    [&](int value){return references[value]==0;}),retained.end());
            }
            descendants[node]=std::make_shared<const std::vector<int>>(std::move(retained));
        }
        std::vector<uint8_t> reachable(states.size(),0);reachable[0]=1;
        for (size_t i=0;i<transitions.size();++i)
            if (!removed[i] && reachable[transitions[i].source]) reachable[transitions[i].target]=1;
        frontier.erase(std::remove_if(frontier.begin(),frontier.end(),
            [&](const auto& branch){return !reachable[branch->node];}),frontier.end());
        size_t kept=0;
        for (size_t i=0;i<transitions.size();++i)
            if (!removed[i] && reachable[transitions[i].source]) {
                if (kept!=i) transitions[kept]=std::move(transitions[i]);
                ++kept;
            }
        transitions.resize(kept);
        return changed;
    }
    void prune(const std::vector<BranchPtr>& frontier) {
        std::vector<uint8_t> reachable(states.size(),0);reachable[0]=1;
        for (const auto& edge:transitions)
            if (reachable[edge.source]) reachable[edge.target]=1;
        std::vector<uint8_t> retained(states.size(),0);
        retained[0]=1;
        for (const auto& branch:frontier) if (reachable[branch->node]) retained[branch->node]=1;
        for (const auto& stop:stops)
            if (stop.state>=0 && reachable[stop.state]) retained[stop.state]=1;
        // Every edge points to a newly appended state. Reverse traversal keeps
        // every incoming path through joins, not just one representative path.
        for (auto edge=transitions.rbegin();edge!=transitions.rend();++edge)
            if (reachable[edge->source] && retained[edge->target]) retained[edge->source]=1;
        std::vector<int> indices(states.size(),-1);
        size_t kept=0;
        for (size_t i=0;i<states.size();++i) if (retained[i]) {
            indices[i]=(int)kept;
            if (kept!=i) states[kept]=std::move(states[i]);
            ++kept;
        }
        states.resize(kept);
        kept=0;
        for (size_t i=0;i<transitions.size();++i)
            if (reachable[transitions[i].source] && retained[transitions[i].target]) {
                auto& edge=transitions[i];
                edge.source=indices[edge.source];edge.target=indices[edge.target];
                if (kept!=i) transitions[kept]=std::move(edge);
                ++kept;
            }
        transitions.resize(kept);
        kept=0;cap_summaries.clear();
        for (size_t i=0;i<stops.size();++i) {
            auto& stop=stops[i];
            if (stop.state>=0) {
                if (!reachable[stop.state]) continue;
                stop.state=indices[stop.state];
            } else if (stop.reason=="capped") {
                cap_summaries.emplace(std::make_tuple(stop.seed,stop.pass,stop.position,
                    stop.reason,stop.stage,stop.count,stop.limit),kept);
            }
            if (kept!=i) stops[kept]=std::move(stops[i]);
            ++kept;
        }
        stops.resize(kept);
        for (const auto& branch:frontier) branch->node=indices[branch->node];
        needs_pruning=false;
    }
};

using ScheduledFrontier=std::unordered_map<const ScheduledBranch*,BranchPtr,ScheduledKeyHash,ScheduledKeyEqual>;

BranchPtr commit_scheduled_fragment(const BranchPtr& parent,const std::shared_ptr<const Fragment>& saved,
                                    const std::vector<Pair>& deferred,const PlacementGroup& placement_generators,
                                    const Source& source,ScheduledGraph& graph) {
    const auto& fragment=*saved;
    auto child=std::make_shared<ScheduledBranch>(*parent);
    for (size_t i=0;i<fragment.witness.size();++i) {
        int r=fragment.source->atoms[i],p=fragment.target->atoms[fragment.witness[i]];
        if (child->mapping[r]<0) {child->mapping[r]=p;++child->mapped;}
    }
    child->deferred=deferred;
    std::vector<uint8_t> in_merge(source.n,0);
    std::unordered_set<int> touched;
    for (int r:fragment.source->atoms) {
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
    int island=child->labels[fragment.source->atoms.front()];
    child->node=graph.state(*child);
    graph.transitions.push_back({parent->node,child->node,graph.seed,graph.pass,graph.position,island,
        saved,placement_generators});
    return child;
}

ScheduledGraph schedule_fragments(const Source& source,const Target& target,
                                  const std::vector<int>& order,double graph_floor,double iso_tol,
                                  long max_branches,const std::vector<Pair>& anchors,
                                  const std::vector<int>& core,bool stop_at_core,
                                  FragmentRepair* repair,const std::vector<Pair>& cuts,uint64_t random_seed,double self_tolerance) {
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
            std::unordered_set<const ScheduledBranch*> advanced;
            seen.reserve(branches.size());next.reserve(branches.size());
            auto admit=[&](const std::vector<BranchPtr>& subtree,const ScheduledBranch& parent) {
                std::vector<BranchPtr> additions;
                ScheduledFrontier local;
                for (const auto& b:subtree) {
                    bool changed=b->key!=parent.key;
                    auto prior=seen.find(b.get());
                    if (prior!=seen.end()) {
                        graph.join(*prior->second,*b);
                        if (changed) advanced.insert(prior->second.get());
                        continue;
                    }
                    auto sibling=local.find(b.get());
                    if (sibling!=local.end()) {
                        graph.join(*sibling->second,*b);
                        if (changed) advanced.insert(sibling->second.get());
                        continue;
                    }
                    local.emplace(b.get(),b);additions.push_back(b);
                    if (changed) advanced.insert(b.get());
                }
                for (const auto& b:additions) {seen.emplace(b.get(),b);next.push_back(b);}
            };
            auto carry=[&](const BranchPtr& b) {
                auto prior=seen.find(b.get());
                if (prior!=seen.end()) {graph.join(*prior->second,*b);return;}
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
                                      &branch->islands,branch->deferred,anchor,random_seed,self_tolerance).first
                    : grow_island(source,target,seed,branch->mapping,graph_floor,iso_tol,1,max_branches,
                                  &branch->islands,branch->deferred,anchor,random_seed,nullptr,nullptr,nullptr,self_tolerance);
                graph.max_growth_candidates=std::max(graph.max_growth_candidates,result.profile.max_cands_before);
                if (result.capped) {
                    graph.stop(*branch,"capped","fragment_growth",result.cap.count,result.cap.limit);
                    continue;
                }
                if (result.fragments.empty()) {carry(branch);continue;}
                std::vector<BranchPtr> subtree;
                for (const auto& fragment:result.fragments) {
                    auto child=commit_scheduled_fragment(branch,fragment,result.deferred_edges,
                        result.placement_generators,source,graph);
                    subtree.push_back(std::move(child));
                }
                admit(subtree,*branch);
            }
            seen.clear();
            branches.clear();
            graph.cap_frontier(next,max_branches);
            if (std::any_of(next.begin(),next.end(),[&](const auto& b){return advanced.count(b.get());}))
                progressed=true;
            std::stable_sort(next.begin(),next.end(),[](const auto& a,const auto& b){return a->mapped>b->mapped;});
            branches=std::move(next);
            if (graph.needs_pruning) graph.prune(branches);
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
