// Exact set transport under compact permutation generators. No orbit of sets
// or complete permutation group is materialized.
#pragma once
#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace permutation_native {
using Permutation=std::vector<int>;
using Generators=std::vector<Permutation>;

inline Permutation identity(int degree) {
    Permutation result(degree);std::iota(result.begin(),result.end(),0);return result;
}
inline Permutation compose(const Permutation& left,const Permutation& right) {
    Permutation result(right.size());
    for (size_t i=0;i<right.size();++i) result[i]=left[right[i]];
    return result;
}
inline Permutation inverse(const Permutation& action) {
    Permutation result(action.size());
    for (size_t i=0;i<action.size();++i) result[action[i]]=(int)i;
    return result;
}

struct Transversal {
    std::vector<int> points;
    Generators actions;
    Transversal(int degree,const Generators& generators,int fixed):actions(degree) {
        points.push_back(fixed);actions[fixed]=identity(degree);
        for (size_t i=0;i<points.size();++i) for (const auto& action:generators) {
            int image=action[points[i]];
            if (actions[image].empty()) {
                actions[image]=compose(action,actions[points[i]]);points.push_back(image);
            }
        }
    }
    Generators stabilizer(const Generators& generators) const {
        std::set<Permutation> unique;
        const auto unit=identity((int)actions.size());
        Generators inverses(actions.size());
        for (int p:points) inverses[p]=inverse(actions[p]);
        for (int p:points) for (const auto& action:generators) {
            auto residual=compose(inverses[action[p]],compose(action,actions[p]));
            if (residual!=unit) unique.insert(std::move(residual));
        }
        return Generators(unique.begin(),unique.end());
    }
};

struct StabilizerChain {
    int fixed=-1;
    std::unique_ptr<Transversal> transversal;
    std::unique_ptr<StabilizerChain> next;
    StabilizerChain(int degree,const Generators& generators) {
        for (int p=0;p<degree && fixed<0;++p)
            for (const auto& action:generators) if (action[p]!=p) {fixed=p;break;}
        if (fixed<0) return;
        transversal=std::make_unique<Transversal>(degree,generators,fixed);
        next=std::make_unique<StabilizerChain>(degree,transversal->stabilizer(generators));
    }
    bool contains(const std::vector<uint8_t>& first,const std::vector<uint8_t>& second) const {
        if (first==second) return true;
        if (fixed<0) return false;
        for (int p:transversal->points) if (first[fixed]==second[p]) {
            std::vector<uint8_t> moved(second.size());
            const auto& action=transversal->actions[p];
            for (size_t i=0;i<second.size();++i) moved[i]=second[action[i]];
            if (next->contains(first,moved)) return true;
        }
        return false;
    }
};

class SetOrbit {
    struct Component {
        Permutation atoms;
        Generators generators;
        mutable std::unique_ptr<StabilizerChain> chain;
    };
    int degree;
    std::vector<Component> components;
    std::vector<uint8_t> moving;
    std::vector<uint8_t> mask(const Permutation& atoms) const {
        std::vector<uint8_t> result(degree,0);
        for (int atom:atoms) {
            if (atom<0 || atom>=degree || result[atom])
                throw std::invalid_argument("set must contain distinct atoms in its permutation frame");
            result[atom]=1;
        }
        return result;
    }
public:
    SetOrbit(int n,const Generators& generators):degree(n),moving(n>=0?n:0,0) {
        if (n<0) throw std::invalid_argument("negative permutation degree");
        auto parent=identity(n);
        auto root=[&](int p) {while (parent[p]!=p) p=parent[p];return p;};
        for (const auto& action:generators) {
            if ((int)action.size()!=n) throw std::invalid_argument("generator degree differs from its frame");
            auto sorted=action;std::sort(sorted.begin(),sorted.end());
            if (sorted!=identity(n)) throw std::invalid_argument("invalid permutation generator");
            int first=-1;
            for (int p=0;p<n;++p) if (action[p]!=p) {
                moving[p]=1;
                if (first<0) first=p;
                else parent[root(p)]=root(first);
            }
        }
        // A generator's entire support stays together: separate atom orbits
        // can still have correlated motion and must not be split independently.
        std::map<int,size_t> indices;
        for (int p=0;p<n;++p) if (moving[p]) {
            auto [entry,added]=indices.emplace(root(p),components.size());
            if (added) components.emplace_back();
            components[entry->second].atoms.push_back(p);
        }
        for (auto& component:components) {
            Permutation position(n,-1);
            for (size_t i=0;i<component.atoms.size();++i) position[component.atoms[i]]=(int)i;
            std::set<Permutation> unique;
            for (const auto& action:generators) {
                Permutation local;bool changed=false;
                for (int p:component.atoms) {local.push_back(position[action[p]]);changed|=action[p]!=p;}
                if (changed) unique.insert(std::move(local));
            }
            component.generators.assign(unique.begin(),unique.end());
        }
    }
    bool contains(const Permutation& first,const Permutation& second) const {
        if (first.size()!=second.size()) return false;
        auto left=mask(first),right=mask(second);
        for (int p=0;p<degree;++p) if (!moving[p] && left[p]!=right[p]) return false;
        for (const auto& component:components) {
            std::vector<uint8_t> a,b;
            for (int p:component.atoms) {a.push_back(left[p]);b.push_back(right[p]);}
            if (std::count(a.begin(),a.end(),1)!=std::count(b.begin(),b.end(),1)) return false;
            if (a==b) continue;
            if (!component.chain) component.chain=std::make_unique<StabilizerChain>(
                (int)component.atoms.size(),component.generators);
            if (!component.chain->contains(a,b)) return false;
        }
        return true;
    }
    Generators stabilizer(const Permutation& fixed) const {
        auto locked=mask(fixed);
        Generators result;
        for (const auto& component:components) {
            Permutation local_fixed;
            for (size_t i=0;i<component.atoms.size();++i)
                if (locked[component.atoms[i]]) local_fixed.push_back((int)i);
            if (local_fixed.size()==component.atoms.size()) continue;
            auto generators=component.generators;
            if (!local_fixed.empty()) {
                generators=Transversal((int)component.atoms.size(),generators,local_fixed.front()).stabilizer(generators);
                // Fixing a point can remove the action coupling several
                // supports. Repartition the residual generators before the
                // remaining locks, retaining every generator's joint support.
                if (local_fixed.size()>1 && !generators.empty()) {
                    local_fixed.erase(local_fixed.begin());
                    generators=SetOrbit((int)component.atoms.size(),generators).stabilizer(local_fixed);
                    std::sort(generators.begin(),generators.end());
                }
            }
            for (const auto& local:generators) {
                auto action=identity(degree);
                for (size_t i=0;i<component.atoms.size();++i) action[component.atoms[i]]=component.atoms[local[i]];
                result.push_back(std::move(action));
            }
        }
        return result;
    }
};
} // namespace permutation_native
