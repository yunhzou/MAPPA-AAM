// Compact permutation operations shared by event preparation and queries.
// The stabilizer chain stores point transversals, never all group elements.
#pragma once
#include "set_orbit.h"
#include <mutex>
#include <optional>

namespace permutation_native {

inline bool is_identity(const Permutation& action) {
    for (size_t i = 0; i < action.size(); ++i) if (action[i] != static_cast<int>(i)) return false;
    return true;
}

// Incremental Schreier closure. Sifting each residual before insertion keeps
// the generators compact even when a raw Schreier list contains many repeats.
class GroupChain {
    int degree, point = -1;
    Generators generators;
    std::unique_ptr<Transversal> transversal;
    std::unique_ptr<GroupChain> next;
public:
    explicit GroupChain(int n): degree(n) {}

    bool contains(const Permutation& action) const {
        if (is_identity(action)) return true;
        if (point < 0 || transversal->actions[action[point]].empty()) return false;
        auto remainder = compose(inverse(transversal->actions[action[point]]), action);
        return next ? next->contains(remainder) : is_identity(remainder);
    }

    bool insert(const Permutation& action) {
        if (contains(action)) return false;
        if (point < 0) {
            for (int i = 0; i < degree; ++i) if (action[i] != i) {point = i; break;}
        }
        generators.push_back(action);
        transversal = std::make_unique<Transversal>(degree, generators, point);
        Generators inverses(degree);
        for (int p : transversal->points) inverses[p] = inverse(transversal->actions[p]);
        for (int p : transversal->points) for (const auto& generator : generators) {
            auto remainder = compose(inverses[generator[p]], compose(generator, transversal->actions[p]));
            if (!is_identity(remainder)) {
                if (!next) next = std::make_unique<GroupChain>(degree);
                next->insert(remainder);
            }
        }
        return true;
    }

    void factors(std::vector<int>& sizes) const {
        if (!transversal) return;
        sizes.push_back(static_cast<int>(transversal->points.size()));
        if (next) next->factors(sizes);
    }

    void transversals(std::vector<std::vector<std::pair<int, Permutation>>>& result) const {
        if (!transversal) return;
        std::vector<std::pair<int, Permutation>> level;
        for (int p : transversal->points) level.emplace_back(p, transversal->actions[p]);
        result.push_back(std::move(level));
        if (next) next->transversals(result);
    }
};

inline Generators reduce_generators(int degree, const Generators& generators) {
    GroupChain chain(degree);
    Generators result;
    for (const auto& action : generators) if (chain.insert(action)) result.push_back(action);
    return result;
}

class Group {
    mutable std::once_flag chain_ready;
    mutable std::unique_ptr<GroupChain> chain;
    void prepare_chain() const {
        std::call_once(chain_ready, [&] {
            auto prepared = std::make_unique<GroupChain>(degree);
            for (const auto& action : generators) prepared->insert(action);
            chain = std::move(prepared);
        });
    }
public:
    const int degree;
    Generators generators;
    explicit Group(int n, Generators values): degree(n) {
        if (degree < 0) throw std::invalid_argument("negative permutation degree");
        const auto unit = identity(degree);
        std::set<Permutation> seen;
        for (auto& action : values) {
            auto sorted = action;
            std::sort(sorted.begin(), sorted.end());
            if (sorted != unit) throw std::invalid_argument("invalid permutation generator");
            if (action != unit && seen.insert(action).second) generators.push_back(std::move(action));
        }
    }

    Permutation support() const {
        std::vector<bool> moved(degree, false);
        for (const auto& action : generators) for (int p = 0; p < degree; ++p) moved[p] = moved[p] || action[p] != p;
        Permutation result;
        for (int p = 0; p < degree; ++p) if (moved[p]) result.push_back(p);
        return result;
    }

    std::vector<Permutation> orbits() const {
        auto parent = identity(degree);
        auto root = [&](int p) {
            while (parent[p] != p) {parent[p] = parent[parent[p]]; p = parent[p];}
            return p;
        };
        for (const auto& action : generators) for (int p = 0; p < degree; ++p) parent[root(p)] = root(action[p]);
        std::map<int, size_t> positions;
        std::vector<Permutation> result;
        for (int p = 0; p < degree; ++p) {
            auto [entry, added] = positions.emplace(root(p), result.size());
            if (added) result.emplace_back();
            result[entry->second].push_back(p);
        }
        return result;
    }

    Permutation orbit(int point) const {
        if (point < 0 || point >= degree) throw std::invalid_argument("point outside group frame");
        std::vector<bool> seen(degree, false);
        Permutation result{point}; seen[point] = true;
        for (size_t i = 0; i < result.size(); ++i) for (const auto& action : generators) {
            int p = action[result[i]];
            if (!seen[p]) {seen[p] = true; result.push_back(p);}
        }
        std::sort(result.begin(), result.end());
        return result;
    }

    std::vector<std::vector<std::pair<int, int>>> pair_orbits(
            std::optional<std::vector<std::pair<int, int>>> selected = std::nullopt) const {
        std::vector<std::pair<int, int>> starts;
        if (selected) starts = *selected;
        else for (int i = 0; i < degree; ++i) for (int j = i + 1; j < degree; ++j) starts.emplace_back(i, j);
        std::sort(starts.begin(), starts.end());
        std::vector<bool> seen(static_cast<size_t>(degree) * degree, false);
        std::vector<std::vector<std::pair<int, int>>> result;
        for (const auto& [first, second] : starts) {
            if (first < 0 || second <= first || second >= degree) throw std::invalid_argument("pair outside group frame");
            if (seen[static_cast<size_t>(first) * degree + second]) continue;
            std::vector<std::pair<int, int>> orbit{{first, second}};
            seen[static_cast<size_t>(first) * degree + second] = true;
            for (size_t i = 0; i < orbit.size(); ++i) for (const auto& action : generators) {
                int a = action[orbit[i].first], b = action[orbit[i].second];
                if (a > b) std::swap(a, b);
                const size_t key = static_cast<size_t>(a) * degree + b;
                if (!seen[key]) {seen[key] = true; orbit.emplace_back(a, b);}
            }
            std::sort(orbit.begin(), orbit.end());
            result.push_back(std::move(orbit));
        }
        return result;
    }

    std::vector<std::pair<int, Permutation>> transversal(int point) const {
        if (point < 0 || point >= degree) throw std::invalid_argument("point outside group frame");
        Transversal table(degree, generators, point);
        std::vector<std::pair<int, Permutation>> result;
        for (int p : table.points) result.emplace_back(p, table.actions[p]);
        return result;
    }

    std::shared_ptr<Group> stabilizer(const Permutation& fixed) const {
        auto remaining = generators;
        for (int point : fixed) {
            if (point < 0 || point >= degree) throw std::invalid_argument("point outside group frame");
            if (std::all_of(remaining.begin(), remaining.end(), [&](const auto& g) {return g[point] == point;})) continue;
            Transversal table(degree, remaining, point);
            remaining = reduce_generators(degree, table.stabilizer(remaining));
        }
        return std::make_shared<Group>(degree, std::move(remaining));
    }

    bool contains(const Permutation& action) const {
        auto sorted = action;
        std::sort(sorted.begin(), sorted.end());
        if (sorted != identity(degree)) return false;
        prepare_chain();
        return chain->contains(action);
    }

    std::vector<int> order_factors() const {
        prepare_chain();
        std::vector<int> sizes;
        chain->factors(sizes);
        return sizes;
    }

    std::vector<std::vector<std::pair<int, Permutation>>> transversals() const {
        prepare_chain();
        std::vector<std::vector<std::pair<int, Permutation>>> result;
        chain->transversals(result);
        return result;
    }
};

inline std::vector<Permutation> merged_blocks(const Permutation& atoms, const std::vector<Permutation>& joins) {
    const int degree = atoms.empty() ? 0 : *std::max_element(atoms.begin(), atoms.end()) + 1;
    auto parent = identity(degree);
    auto root = [&](int p) {while (parent[p] != p) {parent[p] = parent[parent[p]]; p = parent[p];} return p;};
    for (const auto& joined : joins) if (!joined.empty())
        for (int p : joined) parent[root(p)] = root(joined.front());
    std::map<int, size_t> positions;
    std::vector<Permutation> blocks;
    for (int p : atoms) {
        auto [entry, added] = positions.emplace(root(p), blocks.size());
        if (added) blocks.emplace_back();
        blocks[entry->second].push_back(p);
    }
    return blocks;
}

inline std::map<int, int> order_signature(const Group& group) {
    std::map<int, int> powers;
    for (int value : group.order_factors()) {
        for (int factor = 2; factor * factor <= value; ++factor)
            while (value % factor == 0) {++powers[factor]; value /= factor;}
        if (value > 1) ++powers[value];
    }
    return powers;
}

inline std::vector<std::shared_ptr<Group>> independent_factors(const std::shared_ptr<Group>& group,
        const std::vector<std::vector<double>>* adjacency) {
    auto moved = group->support();
    if (moved.empty()) return {};
    std::vector<Permutation> supports;
    for (const auto& action : group->generators) {
        Permutation support;
        for (int p : moved) if (action[p] != p) support.push_back(p);
        supports.push_back(std::move(support));
    }
    auto blocks = merged_blocks(moved, supports);
    if (blocks.size() > 1) {
        std::vector<std::shared_ptr<Group>> result;
        for (const auto& block : blocks) {
            Generators generators;
            for (size_t i = 0; i < supports.size(); ++i)
                if (std::binary_search(block.begin(), block.end(), supports[i].front())) generators.push_back(group->generators[i]);
            auto factors = independent_factors(std::make_shared<Group>(group->degree, std::move(generators)), adjacency);
            result.insert(result.end(), factors.begin(), factors.end());
        }
        return result;
    }
    auto orbits = group->orbits();
    orbits.erase(std::remove_if(orbits.begin(), orbits.end(), [](const auto& orbit) {return orbit.size() == 1;}), orbits.end());
    std::stable_sort(orbits.begin(), orbits.end(), [](const auto& a, const auto& b) {return a.size() < b.size();});
    if (orbits.size() > 1) for (const auto& orbit : orbits) {
        Permutation positions(group->degree, -1), outside;
        for (size_t i = 0; i < orbit.size(); ++i) positions[orbit[i]] = static_cast<int>(i);
        bool correlated = false;
        if (adjacency) for (int p = 0; p < group->degree && !correlated; ++p) if (positions[p] < 0)
            for (int q : orbit) if ((*adjacency)[p][q] != (*adjacency)[p][orbit.front()]) {correlated = true; break;}
        if (correlated) continue;
        for (int p : moved) if (positions[p] < 0) outside.push_back(p);
        auto isolated = group->stabilizer(outside);
        if (isolated->generators.empty()) continue;
        Generators projected;
        for (const auto& action : group->generators) {
            Permutation local;
            for (int p : orbit) local.push_back(positions[action[p]]);
            projected.push_back(std::move(local));
        }
        Group image(static_cast<int>(orbit.size()), std::move(projected));
        if (order_signature(*isolated) == order_signature(image)) {
            auto result = independent_factors(isolated, adjacency);
            auto remaining = independent_factors(group->stabilizer(orbit), adjacency);
            result.insert(result.end(), remaining.begin(), remaining.end());
            return result;
        }
    }
    return {group};
}

inline std::vector<std::shared_ptr<Group>> event_components(const std::vector<std::shared_ptr<Group>>& factors,
        const std::vector<Permutation>& scopes) {
    if (factors.size() <= 1) return factors;
    const int degree = factors.front()->degree;
    Permutation owner(degree, -1);
    for (size_t i = 0; i < factors.size(); ++i) for (int p : factors[i]->support()) owner[p] = static_cast<int>(i);
    std::vector<Permutation> joins;
    for (const auto& scope : scopes) {
        Permutation joined;
        for (int p : scope) {
            if (p < 0 || p >= degree) throw std::invalid_argument("event scope outside group frame");
            if (owner[p] >= 0) joined.push_back(owner[p]);
        }
        joins.push_back(std::move(joined));
    }
    std::vector<std::shared_ptr<Group>> result;
    for (const auto& block : merged_blocks(identity(static_cast<int>(factors.size())), joins)) {
        if (block.size() == 1) {result.push_back(factors[block.front()]); continue;}
        Generators generators;
        for (int i : block) generators.insert(generators.end(), factors[i]->generators.begin(), factors[i]->generators.end());
        result.push_back(std::make_shared<Group>(degree, std::move(generators)));
    }
    return result;
}
} // namespace permutation_native
