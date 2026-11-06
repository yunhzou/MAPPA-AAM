// Exact first-solution matching of accepted growth constraints.
class FragmentQuery {
public:
    const Source& source;
    const Target& target;
    struct Edge { int a,b; double weight,tolerance; };
    using Domains = std::vector<std::vector<int>>;
    std::vector<int> atoms;
    std::vector<Edge> edges;
    Domains domains,incident;
    struct ImageComponents {
        bool connected=false;
        size_t eligible=0;
        std::vector<int> component,sizes;
    };
    std::unique_ptr<ImageComponents> prepared_components;

    const ImageComponents& image_components() {
        if (prepared_components) return *prepared_components;
        prepared_components=std::make_unique<ImageComponents>();
        auto& result=*prepared_components;
        result.component.assign(target.n,-1);
        for (int r:atoms) for (int p:domains[r]) result.component[p]=-2;
        result.eligible=std::count(result.component.begin(),result.component.end(),-2);
        if (!atoms.empty()) {
            std::vector<uint8_t> reached(source.n,0);
            std::vector<int> pending{atoms.front()};reached[atoms.front()]=1;
            for (size_t i=0;i<pending.size();++i) for (int index:incident[pending[i]]) {
                const auto& edge=edges[index];
                int other=edge.a==pending[i] ? edge.b : edge.a;
                if (!reached[other]) {reached[other]=1;pending.push_back(other);}
            }
            result.connected=pending.size()==atoms.size();
        }
        if (result.connected) for (int p=0;p<target.n;++p) if (result.component[p]==-2) {
            int index=(int)result.sizes.size();
            std::vector<int> pending{p};result.component[p]=index;
            for (size_t i=0;i<pending.size();++i) for (int q:target.adj[pending[i]])
                if (result.component[q]==-2) {result.component[q]=index;pending.push_back(q);}
            result.sizes.push_back((int)pending.size());
        }
        return result;
    }

    bool needs_global_symmetry() {
        return image_components().eligible>atoms.size();
    }

    bool preserves_domains(const std::vector<int>& generator,const Domains& allowed) const {
        for (int r:atoms) for (int p:allowed[r])
            if (!contains(allowed[r],generator[p])) return false;
        return true;
    }

    bool supports(const Edge& edge,int p,int q) const {
        return p!=q && target.has_edge(p,q) &&
            growth_edge_supported(edge.weight,target.wbo(p,q),edge.tolerance,
                                  target.bond_cut,source.tolerance_inclusive);
    }

    bool distinct_destinations(const Domains& allowed) const {
        std::vector<int> owner(target.n, -1), seen(target.n, -1), order;
        for (int r : atoms) {
            if (allowed[r].size() == 1) owner[allowed[r][0]] = r;
            else order.push_back(r);
        }
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return allowed[a].size() < allowed[b].size();
        });
        int visit = 0;
        std::function<bool(int)> augment = [&](int r) {
            for (int p : allowed[r]) {
                if (seen[p] == visit) continue;
                seen[p] = visit;
                if (owner[p] < 0 || augment(owner[p])) { owner[p] = r; return true; }
            }
            return false;
        };
        for (int r : order) {
            ++visit;
            if (!augment(r)) return false;
        }
        return true;
    }

    bool reduce(Domains& allowed) const {
        std::vector<size_t> domain_marks(target.n,0);
        size_t mark=0;
        bool changed=true;
        while (changed) {
            changed=false;
            std::vector<uint8_t> occupied(target.n,0);
            for (int r:atoms) if (allowed[r].size()==1) {
                int p=allowed[r][0];
                if (occupied[p]) return false;
                occupied[p]=1;
            }
            for (int r:atoms) {
                auto& values=allowed[r];
                size_t before=values.size();
                if (before>1) values.erase(std::remove_if(values.begin(),values.end(),
                    [&](int p){return occupied[p];}),values.end());
                if (values.empty()) return false;
                changed|=before!=values.size();
            }
            for (const auto& edge:edges) for (const auto& direction:{Pair{edge.a,edge.b},Pair{edge.b,edge.a}}) {
                int a=direction.first,b=direction.second;
                auto& values=allowed[a];
                size_t before=values.size();
                int fixed=allowed[b].size()==1 ? allowed[b][0] : -1;
                // Each stamp marks one domain without clearing the scratch array.
                ++mark;
                if (fixed<0) for (int q:allowed[b]) domain_marks[q]=mark;
                values.erase(std::remove_if(values.begin(),values.end(),[&](int p) {
                    if (fixed>=0) return !supports(edge,p,fixed);
                    const auto& neighbors=target.adj[p];
                    return std::none_of(neighbors.begin(),neighbors.end(),[&](int q) {
                        return domain_marks[q]==mark && supports(edge,p,q);
                    });
                }),values.end());
                if (values.empty()) return false;
                changed|=before!=values.size();
            }
        }
        ++mark;
        for (int r:atoms) for (int p:allowed[r]) domain_marks[p]=mark;
        return (size_t)std::count(domain_marks.begin(),domain_marks.end(),mark)>=atoms.size()
            && distinct_destinations(allowed);
    }

    int variable(const Domains& allowed) const {
        int selected=-1;
        for (int r:atoms) if (allowed[r].size()>1 && (selected<0 ||
            allowed[r].size()<allowed[selected].size() ||
            (allowed[r].size()==allowed[selected].size() && incident[r].size()>incident[selected].size())))
            selected=r;
        return selected;
    }

    std::vector<int> preferred_extension(const std::vector<int>& preferred,const Domains& allowed) const {
        std::vector<int> witness(source.n,-1);
        std::vector<uint8_t> used(target.n,0);
        for (int r:atoms) if (preferred[r]>=0) {
            int p=preferred[r];
            if (!contains(allowed[r],p) || used[p]) return {};
            witness[r]=p;used[p]=1;
        }
        for (const auto& edge:edges)
            if (witness[edge.a]>=0 && witness[edge.b]>=0 && !supports(edge,witness[edge.a],witness[edge.b]))
                return {};
        for (int r:atoms) if (witness[r]<0) {
            for (int p:allowed[r]) if (!used[p]) {
                bool valid=true;
                for (int index:incident[r]) {
                    const auto& edge=edges[index];
                    int other=edge.a==r ? edge.b : edge.a;
                    if (witness[other]>=0 && !supports(edge,p,witness[other])) {valid=false;break;}
                }
                if (valid) {witness[r]=p;used[p]=1;break;}
            }
            if (witness[r]<0) return {};
        }
        return witness;
    }

public:
    long calls=0,nodes=0,rejected=0;
    long global_generators=0,global_pruned=0;
    FragmentQuery(const Source& R,const Target& P,const std::vector<uint8_t>& fragment,
                  const std::vector<Pair>& bonds,const std::vector<int>& mapping,double tolerance)
        :source(R),target(P),domains(R.n),incident(R.n) {
        std::vector<uint8_t> occupied(P.n,0);
        for (int p:mapping) if (p>=0) occupied[p]=1;
        for (int r=0;r<R.n;++r) if (fragment[r]) {
            atoms.push_back(r);
            if (mapping[r]>=0) domains[r].push_back(mapping[r]);
            else if (R.ecode[r]>=0 && R.ecode[r]<(int)P.same_element.size())
                for (int p:P.same_element[R.ecode[r]]) if (!occupied[p]) domains[r].push_back(p);
        }
        for (auto [a,b]:bonds) {
            int index=(int)edges.size();
            edges.push_back({a,b,R.wbo(a,b),R.tolerance(a,b,tolerance)});
            incident[a].push_back(index);incident[b].push_back(index);
        }
        for (int r : atoms) {
            auto& values = domains[r];
            values.erase(std::remove_if(values.begin(), values.end(), [&](int p) {
                return target.adj[p].size() < incident[r].size();
            }), values.end());
        }
    }

    std::vector<int> find(const std::vector<int>* preferred=nullptr) {
        ++calls;
        auto allowed=domains;
        if (preferred) {
            auto witness=preferred_extension(*preferred,allowed);
            if (!witness.empty()) {++nodes;return witness;}
        }
        std::vector<int> witness;
        std::function<bool(Domains)> search=[&](Domains allowed) {
            ++nodes;
            if (!reduce(allowed)) return false;
            int r=variable(allowed);
            if (r<0) {
                witness.assign(source.n,-1);
                for (int a:atoms) witness[a]=allowed[a][0];
                return true;
            }
            const auto choices=allowed[r];
            for (int p:choices) {
                allowed[r]={p};
                if (search(allowed)) return true;
            }
            return false;
        };
        if (search(std::move(allowed))) return witness;
        ++rejected;return {};
    }

    // One representative per target-set orbit. Exclusions act on sets before
    // witness construction; internal assignment permutations are never emitted.
    std::vector<std::vector<int>> images(long limit,const std::vector<int>* preferred=nullptr,
            const PlacementGenerators& generators={}) {
        if (limit<1) throw std::invalid_argument("image witness limit must be positive");
        ++calls;
        std::vector<std::vector<int>> results;
        std::vector<std::vector<uint64_t>> excluded;
        std::vector<std::vector<int>> representative_sets;
        permutation_native::SetOrbit placement_orbits(target.n,generators);
        const auto& prepared=image_components();
        const auto& component=prepared.component;
        const auto& component_sizes=prepared.sizes;
        bool connected=prepared.connected;
        std::vector<uint8_t> forbidden_components;
        size_t forbidden_count=0;
        for (int size:component_sizes) {
            bool too_small=size<(int)atoms.size();
            forbidden_components.push_back(too_small);forbidden_count+=too_small;
        }
        global_generators=(long)generators.size();
        auto emit=[&](std::vector<int> witness) {
            std::vector<int> selected;
            for (int r:atoms) selected.push_back(witness[r]);
            std::sort(selected.begin(),selected.end());
            representative_sets.push_back(std::move(selected));
            // A connected placement filling an eligible target component has
            // no different atom set inside that component. Exclude it at once.
            int index=connected ? component[witness[atoms.front()]] : -1;
            if (index>=0 && component_sizes[index]==(int)atoms.size()) {
                forbidden_components[index]=1;++forbidden_count;
            } else {
                std::vector<uint64_t> image((target.n+63)/64,0);
                for (int r:atoms) {
                    int p=witness[r];
                    image[p/64]|=uint64_t{1}<<(p%64);
                }
                excluded.push_back(std::move(image));
            }
            results.push_back(std::move(witness));
            return results.size()<static_cast<size_t>(limit);
        };
        auto constrain=[&](Domains& allowed) {
            if (excluded.empty() && !forbidden_count) return reduce(allowed);
            for (;;) {
                size_t before=0;
                for (int r:atoms) before+=allowed[r].size();
                if (forbidden_count) for (int r:atoms) {
                    auto& values=allowed[r];
                    values.erase(std::remove_if(values.begin(),values.end(),
                        [&](int p){return forbidden_components[component[p]];}),values.end());
                    if (values.empty()) return false;
                }
                for (const auto& image:excluded) {
                    auto inside=[&](int p) {return (image[p/64]>>(p%64))&uint64_t{1};};
                    int outside=-1;
                    bool several=false;
                    for (int r:atoms) {
                        if (std::none_of(allowed[r].begin(),allowed[r].end(),
                                         [&](int p){return !inside(p);})) continue;
                        if (outside>=0) {several=true;break;}
                        outside=r;
                    }
                    // Equal-size injective images differ only if at least one
                    // source atom goes outside the whole previously emitted set.
                    if (outside<0) return false;
                    if (!several) {
                        auto& values=allowed[outside];
                        values.erase(std::remove_if(values.begin(),values.end(),inside),values.end());
                    }
                }
                if (!reduce(allowed)) return false;
                size_t after=0;
                for (int r:atoms) after+=allowed[r].size();
                if (before==after) break;
            }
            if (!generators.empty() && !representative_sets.empty()) {
                std::vector<uint8_t> used(target.n,0);
                for (int r:atoms) for (int p:allowed[r]) used[p]=1;
                if ((size_t)std::count(used.begin(),used.end(),1)==atoms.size()) {
                    std::vector<int> forced;
                    for (int p=0;p<target.n;++p) if (used[p]) forced.push_back(p);
                    for (const auto& image:representative_sets)
                        if (placement_orbits.contains(image,forced)) {++global_pruned;return false;}
                }
            }
            return true;
        };
        if (preferred) {
            auto witness=preferred_extension(*preferred,domains);
            if (!witness.empty()) {
                ++nodes;
                if (!emit(std::move(witness))) return results;
            }
        }
        std::function<bool(Domains)> search=[&](Domains allowed) {
            ++nodes;
            if (!constrain(allowed)) return true;
            int r=variable(allowed);
            if (r<0) {
                std::vector<int> witness(source.n,-1);
                for (int a:atoms) witness[a]=allowed[a][0];
                return emit(std::move(witness));
            }
            auto choices=allowed[r];
            std::vector<const std::vector<int>*> stabilizers;
            for (const auto& generator:generators)
                if (preserves_domains(generator,allowed)) stabilizers.push_back(&generator);
            if (!stabilizers.empty()) {
                std::vector<uint8_t> visited(target.n,0);
                std::vector<int> representatives;
                for (int p:choices) {
                    if (visited[p]) {++global_pruned;continue;}
                    representatives.push_back(p);
                    std::vector<int> orbit{p};visited[p]=1;
                    for (size_t i=0;i<orbit.size();++i) for (const auto* generator:stabilizers) {
                        int q=(*generator)[orbit[i]];
                        if (!visited[q]) {visited[q]=1;orbit.push_back(q);}
                    }
                }
                choices=std::move(representatives);
            }
            for (int p:choices) {
                if (!contains(allowed[r],p)) continue;
                auto child=allowed;
                child[r]={p};
                size_t before=results.size();
                if (!search(std::move(child))) return false;
                // A descendant's new image exclusion also constrains every
                // unexplored sibling. Continue this DFS instead of restarting.
                if (results.size()!=before && !constrain(allowed)) return true;
            }
            return true;
        };
        search(domains);
        if (results.empty()) ++rejected;
        return results;
    }
};

std::vector<int> FragmentPool::find(const FragmentQuery& query) {
    ++lookups;
    if (query.atoms.empty()) return {};
    int rare=query.atoms.front();
    for (int atom:query.atoms)
        if (by_atom[atom].size()<by_atom[rare].size()) rare=atom;
    for (auto* entry:by_atom[rare]) {
        ++checks;
        const auto& image=entry->image;
        bool valid=true;
        for (int atom:query.atoms)
            if (image[atom]<0 || !contains(query.domains[atom],image[atom])) {valid=false;break;}
        if (!valid) continue;
        for (const auto& edge:query.edges)
            if (!query.supports(edge,image[edge.a],image[edge.b])) {valid=false;break;}
        if (!valid) continue;
        ++hits;
        entries.splice(entries.end(),entries,index.at(entry->key));
        std::vector<int> witness(base.n,-1);
        for (int atom:query.atoms) witness[atom]=image[atom];
        return witness;
    }
    return {};
}

// One concrete assignment during growth. It is not a compressed candidate:
// accepted extensions update only new atoms; old edges remain proven.
class WitnessGrowth {
    const Source& source;
    const Target& target;
    const std::vector<int>& locked;
    const std::vector<uint8_t>& occupied;
    double tolerance;
    std::vector<uint8_t> used;
    std::vector<int> image;
    FragmentPool* pool;

    void assign(int atom,int value) { image[atom]=value;used[value]=1; }
public:
    long calls=0,nodes=0,searches=0,seed_targets=0;

    const std::vector<int>& assignment() const { return image; }

    WitnessGrowth(const Source& source,const Target& target,const std::vector<int>& locked,
                  const std::vector<uint8_t>& occupied,double tolerance,int seed,FragmentPool* pool)
        :source(source),target(target),locked(locked),occupied(occupied),tolerance(tolerance),
         used(target.n,0),image(source.n,-1),pool(pool) {
        if (locked[seed]>=0) {assign(seed,locked[seed]);seed_targets=1;}
        else for (int p=0;p<target.n;++p) if (!occupied[p] && source.ecode[seed]==target.ecode[p]) {
            if (!seed_targets) assign(seed,p);
            ++seed_targets;
        }
    }

    bool extend(int anchor,int atom,double weight,const std::vector<uint8_t>& fragment,
                const std::vector<Pair>& bonds,const std::vector<Pair>* islands,
                const std::vector<int>& island_of) {
        ++calls;
        const int previous=image[anchor];
        const double limit=source.tolerance(anchor,atom,tolerance);
        auto supports=[&](int value) {
            return target.has_edge(previous,value) && growth_edge_supported(
                weight,target.wbo(previous,value),limit,target.bond_cut,source.tolerance_inclusive);
        };
        std::vector<int> added;
        if (locked[atom]>=0) {
            added.push_back(atom);
            if (islands && island_of[atom]>=0)
                for (auto [r,label]:*islands)
                    if (r!=atom && !fragment[r] && label==island_of[atom]) added.push_back(r);
            bool valid=supports(locked[atom]);
            std::vector<int> values;
            for (int r:added) {
                int p=locked[r];
                if (used[p] || std::find(values.begin(),values.end(),p)!=values.end()) valid=false;
                values.push_back(p);
            }
            if (valid) {
                for (int r:added) assign(r,locked[r]);
                ++nodes;return true;
            }
        } else {
            for (int p:target.adj[previous])
                if (!used[p] && !occupied[p] && target.ecode[p]==source.ecode[atom] && supports(p)) {
                    assign(atom,p);++nodes;return true;
                }
        }
        auto enlarged=fragment;enlarged[atom]=1;
        for (int r:added) enlarged[r]=1;
        auto edges=bonds;edges.push_back({anchor,atom});
        FragmentQuery query(source,target,enlarged,edges,locked,tolerance);
        auto next=pool ? pool->find(query) : std::vector<int>{};
        if (next.empty()) {++searches;next=query.find();nodes+=query.nodes;}
        if (next.empty()) return false;
        image=std::move(next);
        std::fill(used.begin(),used.end(),0);
        for (int p:image) if (p>=0) used[p]=1;
        return true;
    }
};
