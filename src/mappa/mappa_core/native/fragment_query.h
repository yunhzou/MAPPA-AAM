// Exact first-solution search over accepted growth edges. No symmetry expansion.
class FragmentQuery {
    const Source& source;
    const Target& target;
    struct Edge { int a,b; double weight,tolerance; };
    std::vector<int> atoms;
    std::vector<Edge> edges;
    std::vector<std::vector<int>> domains, incident;

    bool supports(const Edge& e,int p,int q) const {
        return p!=q && growth_edge_supported(e.weight,target.wbo(p,q),e.tolerance,
                                             target.bond_cut,source.tolerance_inclusive);
    }
public:
    long calls=0,nodes=0,rejected=0;
    FragmentQuery(const Source& R,const Target& P,const std::vector<uint8_t>& fragment,
                  const std::vector<Pair>& bonds,const std::vector<int>& locked,double tolerance)
        :source(R),target(P),domains(R.n),incident(R.n) {
        std::vector<uint8_t> occupied(P.n,0);
        for (int p:locked) if (p>=0) occupied[p]=1;
        for (int r=0;r<R.n;++r) if (fragment[r]) {
            atoms.push_back(r);
            if (locked[r]>=0) domains[r].push_back(locked[r]);
            else if (R.ecode[r]>=0 && R.ecode[r]<(int)P.same_element.size())
                for (int p:P.same_element[R.ecode[r]]) if (!occupied[p]) domains[r].push_back(p);
        }
        for (auto [a,b]:bonds) {
            int index=(int)edges.size();
            edges.push_back({a,b,R.wbo(a,b),R.tolerance(a,b,tolerance)});
            incident[a].push_back(index);incident[b].push_back(index);
        }
    }

    std::vector<int> find(const Cand* candidate=nullptr,const std::vector<int>* preferred=nullptr) {
        ++calls;
        auto allowed=domains;
        if (candidate) {
            const auto index=block_indexes(*candidate,source.n,target.n);
            for (int r:candidate->mapped) {
                int block=index.r_to_block[r];
                auto& values=allowed[r];
                values.erase(std::remove_if(values.begin(),values.end(),[&](int p) {
                    return block<0 ? p!=candidate->img[r] : !contains(candidate->blocks[block].p,p);
                }),values.end());
            }
            if (!preferred) preferred=&candidate->img;
        }
        // Full-fragment arc consistency prunes future constraints before
        // compressed final matching constructs the corresponding prefix.
        bool changed=true;
        while (changed) {
            changed=false;
            for (const auto& e:edges) for (const auto& direction:{Pair{e.a,e.b},Pair{e.b,e.a}}) {
                int a=direction.first,b=direction.second;
                auto& values=allowed[a];
                size_t before=values.size();
                values.erase(std::remove_if(values.begin(),values.end(),[&](int p) {
                    return std::none_of(allowed[b].begin(),allowed[b].end(),
                                        [&](int q){return supports(e,p,q);});
                }),values.end());
                if (values.empty()) {++rejected;return {};}
                changed|=values.size()!=before;
            }
        }
        std::vector<uint8_t> reachable(target.n,0);
        for (int r:atoms) {
            if (allowed[r].empty()) {++rejected;return {};}
            for (int p:allowed[r]) reachable[p]=1;
        }
        if ((size_t)std::count(reachable.begin(),reachable.end(),1)<atoms.size()) {
            ++rejected;return {};
        }
        std::vector<int> witness(source.n,-1);
        std::vector<uint8_t> used(target.n,0);
        std::function<bool(size_t)> search=[&](size_t assigned) {
            ++nodes;
            if (assigned==atoms.size()) return true;
            int selected=-1;
            std::vector<int> choices;
            for (int r:atoms) if (witness[r]<0) {
                std::vector<int> values;
                for (int p:allowed[r]) if (!used[p]) {
                    bool valid=true;
                    for (int index:incident[r]) {
                        const auto& e=edges[index];
                        int other=e.a==r ? e.b : e.a;
                        if (witness[other]>=0 && !supports(e,p,witness[other])) {valid=false;break;}
                    }
                    if (valid) values.push_back(p);
                }
                if (values.empty()) return false;
                if (selected<0 || values.size()<choices.size() ||
                    (values.size()==choices.size() && incident[r].size()>incident[selected].size())) {
                    selected=r;choices=std::move(values);
                }
            }
            if (preferred) {
                auto hint=std::find(choices.begin(),choices.end(),(*preferred)[selected]);
                if (hint!=choices.end()) std::rotate(choices.begin(),hint,hint+1);
            }
            for (int p:choices) {
                witness[selected]=p;used[p]=1;
                if (search(assigned+1)) return true;
                used[p]=0;witness[selected]=-1;
            }
            return false;
        };
        if (search(0)) return witness;
        ++rejected;return {};
    }
};
