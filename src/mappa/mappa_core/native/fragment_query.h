// Certify the whole structural subgroup against the original weights. Testing
// generators separately would be invalid for a non-transitive tolerance.
bool group_within_tolerance(const PreparedSelf& self,
                           const std::vector<std::vector<int>>& generators) {
    if (generators.empty() || self.weights.empty()) return true;
    auto [low,high]=std::minmax_element(self.weights.begin(),self.weights.end());
    if (*high-*low<=self.tolerance) return true;
    const size_t n=self.atoms.size();
    std::vector<uint8_t> seen(n*n,0);
    for (size_t a=0;a<n;++a) for (size_t b=a+1;b<n;++b) if (!seen[a*n+b]) {
        std::vector<Pair> pending{{(int)a,(int)b}};
        seen[a*n+b]=1;
        double minimum=self.weights[a*n+b],maximum=minimum;
        for (size_t i=0;i<pending.size();++i) for (const auto& generator:generators) {
            int p=generator[pending[i].first],q=generator[pending[i].second];
            if (p>q) std::swap(p,q);
            if (seen[p*n+q]) continue;
            seen[p*n+q]=1;pending.push_back({p,q});
            minimum=std::min(minimum,self.weights[p*n+q]);
            maximum=std::max(maximum,self.weights[p*n+q]);
            if (maximum-minimum>self.tolerance) return false;
        }
    }
    return true;
}

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
        const auto& components=image_components();
        if (components.eligible<=atoms.size()) return false;
        return !components.connected || std::any_of(components.sizes.begin(),components.sizes.end(),
            [&](int size){return size>(int)atoms.size();});
    }

    bool preserves_domains(const std::vector<int>& generator,const Domains& allowed) const {
        for (int r:atoms) for (int p:allowed[r])
            if (!contains(allowed[r],generator[p])) return false;
        return true;
    }

    std::vector<std::vector<int>> image_generators(const PreparedSelf* global,bool whole_group_allowed) const {
        std::vector<std::vector<int>> generators;
        if (!global) return generators;
        for (const auto& generator:global->generators) {
            if (!preserves_domains(generator,domains)) continue;
            bool valid=true;
            for (const auto& edge:edges) {
                for (int p:domains[edge.a]) for (int q:target.adj[p])
                    if (contains(domains[edge.b],q) && supports(edge,p,q) &&
                        !supports(edge,generator[p],generator[q])) {valid=false;break;}
                if (!valid) break;
            }
            if (valid) generators.push_back(generator);
        }
        if (!whole_group_allowed && !group_within_tolerance(*global,generators)) generators.clear();
        return generators;
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
    long global_generators=0,global_images=0,global_pruned=0;
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

    // One witness per target atom set. Callers request cap+1 witnesses to
    // distinguish a complete result from a capped subtree.
    std::vector<std::vector<int>> images(long limit,const std::vector<int>* preferred=nullptr,
            const PreparedSelf* global=nullptr,bool whole_group_allowed=false) {
        if (limit<1) throw std::invalid_argument("image witness limit must be positive");
        ++calls;
        std::vector<std::vector<int>> results;
        std::vector<std::vector<uint64_t>> excluded;
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
        const auto generators=image_generators(global,whole_group_allowed);
        global_generators=(long)generators.size();
        std::set<std::vector<int>> emitted_images;
        auto emit_one=[&](std::vector<int> witness) {
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
        auto emit=[&](std::vector<int> witness) {
            if (generators.empty()) return emit_one(std::move(witness));
            std::vector<int> image;
            for (int r:atoms) image.push_back(witness[r]);
            std::sort(image.begin(),image.end());
            if (!emitted_images.insert(image).second) return true;
            std::vector<size_t> pending{results.size()};
            if (!emit_one(std::move(witness))) return false;
            for (size_t i=0;i<pending.size();++i) for (const auto& generator:generators) {
                image.clear();
                for (int r:atoms) image.push_back(generator[results[pending[i]][r]]);
                std::sort(image.begin(),image.end());
                if (!emitted_images.insert(image).second) continue;
                // The set is new before a witness is transformed. No orbit of
                // internal bijections is ever constructed.
                auto moved=results[pending[i]];
                for (int r:atoms) moved[r]=generator[moved[r]];
                pending.push_back(results.size());++global_images;
                if (!emit_one(std::move(moved))) return false;
            }
            return true;
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
                if (before==after) return true;
            }
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
                if (!search(std::move(child))) return false;
                // A descendant's new image exclusion also constrains every
                // unexplored sibling. Continue this DFS instead of restarting.
                if (!constrain(allowed)) return true;
            }
            return true;
        };
        search(domains);
        if (results.empty()) ++rejected;
        return results;
    }
};

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

    void assign(int atom,int value) { image[atom]=value;used[value]=1; }
public:
    long calls=0,nodes=0,searches=0,seed_targets=0;

    const std::vector<int>& assignment() const { return image; }

    WitnessGrowth(const Source& source,const Target& target,const std::vector<int>& locked,
                  const std::vector<uint8_t>& occupied,double tolerance,int seed)
        :source(source),target(target),locked(locked),occupied(occupied),tolerance(tolerance),
         used(target.n,0),image(source.n,-1) {
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
        ++searches;
        auto next=query.find();nodes+=query.nodes;
        if (next.empty()) return false;
        image=std::move(next);
        std::fill(used.begin(),used.end(),0);
        for (int p:image) if (p>=0) used[p]=1;
        return true;
    }
};
