// Construct disjoint compressed placements from the completed constraints.
class FragmentCompletion {
    const Source& source;
    const Target& target;
    const std::vector<int>& locked;
    FragmentQuery query;
    using Domains = FragmentQuery::Domains;
    using Generators = std::vector<std::vector<int>>;
    std::vector<std::vector<int>> classes, eligible;
    std::vector<int> class_of, locked_source;

    void prepare_classes(const Domains& domains) {
        // Each code describes which accepted-edge predicates allow this pair.
        // Intern only the codes on active target edges, not candidate graphs.
        std::map<std::pair<double,double>,int> predicates;
        for (const auto& edge:query.edges)
            predicates.emplace(std::make_pair(edge.weight,edge.tolerance),0);
        int index=0;
        for (auto& item:predicates) item.second=index++;
        const size_t words=(predicates.size()+63)/64;
        std::vector<uint8_t> allowed((size_t)source.n*target.n,0);
        eligible.resize(target.n);
        for (int r:query.atoms) for (int p:domains[r]) {
            allowed[(size_t)r*target.n+p]=1;
            eligible[p].push_back(r);
        }
        std::vector<int> predicate_bits;
        for (const auto& edge:query.edges) predicate_bits.push_back(predicates.at({edge.weight,edge.tolerance}));
        std::map<std::vector<uint64_t>,int> codes;
        codes[std::vector<uint64_t>(words,0)]=0;
        std::vector<int> pairs((size_t)target.n*target.n,0);
        for (int p=0;p<target.n;++p) if (!eligible[p].empty() && !query.edges.empty())
            for (int q:target.adj[p]) if (p<q && !eligible[q].empty()) {
            std::vector<uint64_t> code(words,0);
            for (size_t e=0;e<query.edges.size();++e) {
                const auto& edge=query.edges[e];
                bool eligible_pair=(allowed[(size_t)edge.a*target.n+p] && allowed[(size_t)edge.b*target.n+q]) ||
                                   (allowed[(size_t)edge.a*target.n+q] && allowed[(size_t)edge.b*target.n+p]);
                int bit=predicate_bits[e];
                if (eligible_pair && growth_edge_supported(edge.weight,target.wbo(p,q),edge.tolerance,
                                                          target.bond_cut,source.tolerance_inclusive))
                    code[bit/64]|=uint64_t{1}<<(bit%64);
            }
            auto [entry,inserted]=codes.emplace(code,(int)codes.size());
            pairs[(size_t)p*target.n+q]=pairs[(size_t)q*target.n+p]=entry->second;
        }
        class_of.assign(target.n,-1);
        using Color=std::tuple<std::vector<int>,int,int,std::string>;
        std::map<Color,std::vector<int>> by_color;
        for (int p=0;p<target.n;++p) {
            if (eligible[p].empty()) continue;
            auto& possibilities=by_color[{eligible[p],target.orbit[p],locked_source[p],target.elem[p]}];
            for (int group:possibilities) {
                int representative=classes[group][0];
                bool equal=true;
                for (int q=0;q<target.n;++q) if (q!=p && q!=representative &&
                    pairs[(size_t)p*target.n+q]!=pairs[(size_t)representative*target.n+q]) {
                    equal=false;break;
                }
                if (equal) {class_of[p]=group;classes[group].push_back(p);break;}
            }
            if (class_of[p]<0) {
                class_of[p]=(int)classes.size();possibilities.push_back(class_of[p]);
                classes.push_back({p});
            }
        }
    }

    Generators generators(const std::vector<int>& assigned) const {
        // If the full target group fixes every active class, all conditioned
        // subgroups do too. No further symmetry solve can change these choices.
        bool fixed=true;
        for (const auto& action:target.global_generators())
            for (size_t group=0;group<classes.size();++group)
                for (int p:classes[group]) if (class_of[action[p]]!=(int)group) fixed=false;
        if (fixed) return {};
        std::vector<std::vector<int>> roles(classes.size());
        for (int r:query.atoms) if (assigned[r]>=0) roles[assigned[r]].push_back(r);
        using Color=std::tuple<std::string,int,int,std::vector<int>,std::vector<int>>;
        std::map<Color,std::vector<int>> colored;
        const std::vector<int> no_roles;
        for (int p=0;p<target.n;++p)
            colored[{target.elem[p],target.component[p],locked_source[p],
                     eligible[p],class_of[p]<0 ? no_roles : roles[class_of[p]]}].push_back(p);
        std::vector<std::vector<int>> cells;
        for (auto& item:colored) cells.push_back(std::move(item.second));
        for (const auto& item:target.edge_cells) cells.push_back(item.second);
        Generators result;
        for (const auto& action:target.nauty->generators(cells)) {
            std::vector<int> projected;
            bool changed=false;
            for (size_t group=0;group<classes.size();++group) {
                int image=class_of[action[classes[group][0]]];
                for (int p:classes[group]) if (class_of[action[p]]!=image)
                    throw std::invalid_argument("supplied target orbit partition is incompatible with exact target symmetry");
                projected.push_back(image);
                changed|=image!=(int)group;
            }
            if (changed) result.push_back(std::move(projected));
        }
        return result;
    }

    bool moved(int group,const Generators& actions) const {
        return std::any_of(actions.begin(),actions.end(),[&](const auto& action) {
            return action[group]!=group;
        });
    }

    std::vector<int> representatives(const std::vector<int>& values,const Generators& actions) const {
        std::vector<uint8_t> remaining(classes.size(),0);
        for (int group:values) remaining[group]=1;
        std::vector<int> result;
        for (int first:values) if (remaining[first]) {
            result.push_back(first);remaining[first]=0;
            std::vector<int> orbit{first};
            for (size_t i=0;i<orbit.size();++i) for (const auto& action:actions) {
                int image=action[orbit[i]];
                if (remaining[image]) {remaining[image]=0;orbit.push_back(image);}
            }
        }
        return result;
    }

    void emit(const std::vector<int>& assigned) {
        std::vector<std::vector<int>> sources(classes.size());
        for (int r:query.atoms) sources[assigned[r]].push_back(r);
        std::vector<int> image(source.n,-1),fixed;
        std::vector<Block> pools;
        for (size_t group=0;group<classes.size();++group) if (!sources[group].empty()) {
            const auto& atoms=sources[group];const auto& values=classes[group];
            for (size_t i=0;i<atoms.size();++i) image[atoms[i]]=values[i];
            if (values.size()>1) pools.push_back(make_block(atoms,values,false));
        }
        for (int r:query.atoms) if (locked[r]>=0) fixed.push_back(r);
        Cand candidate;
        if (!make_cand(image,pools,fixed,target.n,candidate))
            throw std::runtime_error("invalid completed fragment pool");
        results.push_back(std::move(candidate));
    }

    void search(Domains domains,std::vector<int> assigned,const Generators* actions) {
        ++nodes;
        if (!query.reduce(domains)) {++rejected;return;}
        std::vector<std::vector<int>> choices(source.n);
        int variable=-1;
        for (int r:query.atoms) if (assigned[r]<0) {
            for (int p:domains[r]) if (!contains(choices[r],class_of[p])) choices[r].push_back(class_of[p]);
            if (choices[r].size()==1) assigned[r]=choices[r][0];
            else if (variable<0 || choices[r].size()<choices[variable].size() ||
                (choices[r].size()==choices[variable].size() && query.incident[r].size()>query.incident[variable].size()))
                variable=r;
        }
        if (variable<0) {emit(assigned);return;}
        Generators restricted;
        if (!actions || std::any_of(query.atoms.begin(),query.atoms.end(),[&](int r) {
                return assigned[r]>=0 && moved(assigned[r],*actions);
            })) {
            restricted=generators(assigned);actions=&restricted;
        }
        for (int group:representatives(choices[variable],*actions)) {
            assigned[variable]=group;
            domains[variable]=classes[group];
            search(domains,assigned,actions);
        }
    }

public:
    long nodes=0,rejected=0;
    std::vector<Cand> results;
    FragmentCompletion(const Source& source,const Target& target,const std::vector<int>& locked,
                       const std::vector<uint8_t>& fragment,const std::vector<Pair>& bonds,double tolerance)
        :source(source),target(target),locked(locked),query(source,target,fragment,bonds,locked,tolerance),
         locked_source(target.n,-1) {
        for (int r=0;r<source.n;++r) if (locked[r]>=0) locked_source[locked[r]]=r;
    }

    std::vector<Cand> run() {
        auto domains=query.domains;
        if (!query.reduce(domains)) {nodes=rejected=1;return {};}
        prepare_classes(domains);
        std::vector<int> assigned(source.n,-1);
        search(std::move(domains),assigned,nullptr);
        return std::move(results);
    }
};

std::vector<Cand> complete_fragment(const Source& source,const Target& target,
        const std::vector<int>& locked,const std::vector<uint8_t>& fragment,
        const std::vector<Pair>& bonds,double tolerance,GrowProfile& profile) {
    auto started=std::chrono::steady_clock::now();
    FragmentCompletion completion(source,target,locked,fragment,bonds,tolerance);
    auto result=completion.run();
    if (result.empty()) throw std::runtime_error("completed fragment lost its proven matching witness");
    profile.completion_queries=completion.nodes;profile.completion_pruned=completion.rejected;
    profile.completion_generated=(long)result.size();profile.completion_peak_candidates=(long)result.size();
    profile.completion_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    return result;
}
