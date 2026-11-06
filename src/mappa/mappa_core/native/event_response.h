// Sparse original-weight responses over compact binary permutation groups.
// Only atom-pair orbits are visited; complete bijections are never expanded.
#pragma once
#include "assignment_search.h"
#include "permutation_group.h"

namespace event_native {
using Pair = std::pair<int, int>;
using Pairs = std::vector<Pair>;

class Deadline {
    std::chrono::steady_clock::time_point end;
public:
    explicit Deadline(std::optional<double> seconds) {
        end = seconds ? std::chrono::steady_clock::now() +
            std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(*seconds))
            : std::chrono::steady_clock::time_point::max();
    }
    void check() const {
        if (std::chrono::steady_clock::now() >= end) throw SearchExpired();
    }
};

inline std::vector<double> pair_values(const Matrix& weights, const Pairs& pairs) {
    std::vector<double> values;
    values.reserve(pairs.size());
    for (const auto& [a, b] : pairs) values.push_back(weights[a][b]);
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

inline std::optional<bool> constant_response(const std::vector<double>& left,
        const std::vector<double>& right, double threshold) {
    if (std::max(std::abs(left.front() - right.back()), std::abs(left.back() - right.front())) <= threshold) return false;
    if (left.front() - right.back() > threshold || right.front() - left.back() > threshold) return true;
    for (double value : left) {
        auto next = std::lower_bound(right.begin(), right.end(), value);
        if ((next != right.end() && std::abs(value - *next) <= threshold)
            || (next != right.begin() && std::abs(value - *std::prev(next)) <= threshold)) return std::nullopt;
    }
    return true;
}

struct EventTerms {
    int fixed = 0;
    size_t orbit_states = 0;
    std::vector<std::vector<int>> pairs;
};

inline EventTerms event_terms(const Matrix& source, const Matrix& target, const Matrix& limits,
        const Matrix& target_binary, const permutation_native::Group& source_group,
        const permutation_native::Group& target_group, bool two_sided, std::optional<double> seconds) {
    Deadline deadline(seconds);
    const int size = static_cast<int>(source.size());
    EventTerms result;
    deadline.check();
    if (two_sided) {
        std::vector<std::vector<int>> labels(2, std::vector<int>(size * size));
        std::vector<std::vector<std::vector<double>>> values(2);
        for (int side = 0; side < 2; ++side) {
            const auto& group = side ? target_group : source_group;
            const auto& weights = side ? target : source;
            for (const auto& orbit : group.pair_orbits()) {
                deadline.check();
                result.orbit_states += orbit.size();
                for (const auto& [a, b] : orbit) labels[side][a * size + b] = static_cast<int>(values[side].size());
                values[side].push_back(pair_values(weights, orbit));
            }
        }
        std::map<std::tuple<int, int, double>, std::optional<bool>> responses;
        for (int a = 0; a < size; ++a) for (int b = a + 1; b < size; ++b) {
            int left = labels[0][a * size + b], right = labels[1][a * size + b];
            auto key = std::make_tuple(left, right, limits[a][b]);
            auto [entry, added] = responses.try_emplace(key);
            if (added) entry->second = constant_response(values[0][left], values[1][right], limits[a][b]);
            if (entry->second) result.fixed += *entry->second;
            else result.pairs.push_back({a, b, size + a, size + b});
        }
    } else {
        std::vector<bool> remaining(size * size, true);
        auto response = [&](int a, int b, const std::vector<double>& values) {
            std::vector<uint8_t> flags;
            flags.reserve(values.size());
            for (double weight : values) flags.push_back(std::abs(source[a][b] - weight) > limits[a][b]);
            return flags;
        };
        auto constant = [](const auto& flags) {return std::all_of(flags.begin(), flags.end(), [&](auto value) {return value == flags.front();});};
        for (bool present : {false, true}) {
            Pairs pairs;
            for (int a = 0; a < size; ++a) for (int b = a + 1; b < size; ++b)
                if (bool(target_binary[a][b]) == present) pairs.emplace_back(a, b);
            if (pairs.empty()) continue;
            auto values = pair_values(target, pairs);
            for (const auto& [a, b] : pairs) {
                auto flags = response(a, b, values);
                if (constant(flags)) {result.fixed += flags.front(); remaining[a * size + b] = false;}
            }
        }
        Pairs active;
        for (int a = 0; a < size; ++a) for (int b = a + 1; b < size; ++b)
            if (remaining[a * size + b]) active.emplace_back(a, b);
        for (const auto& orbit : target_group.pair_orbits(active)) {
            deadline.check();
            result.orbit_states += orbit.size();
            auto values = pair_values(target, orbit);
            auto signature = response(orbit.front().first, orbit.front().second, values);
            bool invariant = true;
            for (const auto& [a, b] : orbit) if (response(a, b, values) != signature) {invariant = false; break;}
            for (const auto& [a, b] : orbit) if (remaining[a * size + b]) {
                if (invariant) result.fixed += std::abs(source[a][b] - target[a][b]) > limits[a][b];
                else {
                    auto flags = response(a, b, values);
                    if (constant(flags)) result.fixed += flags.front();
                    else result.pairs.push_back({a, b});
                }
            }
        }
        std::sort(result.pairs.begin(), result.pairs.end());
    }
    deadline.check();
    return result;
}

struct InvariantGraph {
    Pairs edges;
    std::vector<int> colors;
    std::vector<std::vector<uint8_t>> profiles;
};

inline InvariantGraph invariant_graph(const Matrix& weights, const Matrix& other,
        const Matrix& binary, const Matrix& limits, std::optional<double> seconds) {
    Deadline deadline(seconds);
    const int size = static_cast<int>(weights.size());
    Pairs pairs;
    for (int a = 0; a < size; ++a) for (int b = a + 1; b < size; ++b) pairs.emplace_back(a, b);
    auto values = pair_values(other, pairs);
    std::map<std::vector<uint8_t>, int> profiles;
    std::map<std::tuple<bool, double, double>, int> comparisons;
    std::set<int> nonedges;
    std::vector<int> ids;
    InvariantGraph result;
    for (const auto& [a, b] : pairs) {
        deadline.check();
        auto key = std::make_tuple(bool(binary[a][b]), weights[a][b], limits[a][b]);
        auto [comparison, fresh] = comparisons.try_emplace(key);
        if (fresh) {
            std::vector<uint8_t> profile{uint8_t(bool(binary[a][b]))};
            for (double value : values) profile.push_back(std::abs(value - weights[a][b]) > limits[a][b]);
            auto [entry, added] = profiles.emplace(profile, static_cast<int>(profiles.size()));
            if (added) result.profiles.push_back(std::move(profile));
            comparison->second = entry->second;
        }
        ids.push_back(comparison->second);
        if (!binary[a][b]) nonedges.insert(comparison->second);
    }
    for (size_t i = 0; i < pairs.size(); ++i) {
        auto [a, b] = pairs[i];
        if (!binary[a][b] && nonedges.size() <= 1) continue;
        int vertex = size + static_cast<int>(result.colors.size());
        result.colors.push_back(ids[i]);
        result.edges.emplace_back(a, vertex);
        result.edges.emplace_back(b, vertex);
    }
    return result;
}
} // namespace event_native
