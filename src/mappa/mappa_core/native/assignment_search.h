// Bounded starting-correspondence search. Domains and original-weight event
// bounds stay native between improvements; no complete assignment table is kept.
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace event_native {
using Domains = std::vector<std::vector<int>>;
using Matrix = std::vector<std::vector<double>>;

struct SearchExpired : std::exception {};

// Shortest augmenting paths for integer assignment costs; forbidden entries
// use the caller's impossible cost. Both event relaxations share this solver.
template<class Check>
int assignment_cost(const std::vector<std::vector<int>>& costs, int impossible, Check check) {
    const int n = static_cast<int>(costs.size());
    std::vector<int> row(n + 1), column(n + 1), owner(n + 1), previous(n + 1);
    for (int i = 1; i <= n; ++i) {
        check();
        owner[0] = i;
        int current = 0;
        std::vector<int> distance(n + 1, impossible);
        std::vector<bool> used(n + 1, false);
        do {
            used[current] = true;
            int next = 0, delta = impossible;
            for (int j = 1; j <= n; ++j) if (!used[j]) {
                int cost = costs[owner[current] - 1][j - 1] - row[owner[current]] - column[j];
                if (cost < distance[j]) {distance[j] = cost; previous[j] = current;}
                if (distance[j] < delta) {delta = distance[j]; next = j;}
            }
            if (!next) return impossible;
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {row[owner[j]] += delta; column[j] -= delta;}
                else distance[j] -= delta;
            }
            current = next;
        } while (owner[current]);
        do {
            int next = previous[current];
            owner[current] = owner[next];
            current = next;
        } while (current);
    }
    int total = 0;
    for (int j = 1; j <= n; ++j) {
        int value = costs[owner[j] - 1][j - 1];
        if (value >= impossible) return impossible;
        total += value;
    }
    return total;
}

// Relax neighbor occupancy separately for each possible atom correspondence,
// then assign atoms injectively. Each real event touches two atoms, so half
// this relaxed incident cost bounds the number of bond events from below.
inline int incidence_bound(const Matrix& source, const Matrix& target, const Matrix& thresholds,
        const std::vector<int>& source_types, const std::vector<int>& target_types,
        const std::vector<std::pair<int, int>>& anchors, double seconds) {
    const int n = static_cast<int>(source.size());
    if (source_types.size() != source.size() || target_types.size() != source.size())
        throw std::invalid_argument("event atom types differ from matrix size");
    for (const auto* matrix : {&source, &target, &thresholds}) {
        if (matrix->size() != source.size()) throw std::invalid_argument("event matrix sizes differ");
        for (const auto& row : *matrix)
            if (row.size() != source.size()) throw std::invalid_argument("event matrices must be square");
    }
    if (!std::isfinite(seconds) || seconds < 0)
        throw std::invalid_argument("bound allowance must be finite and nonnegative");
    const auto deadline = std::chrono::steady_clock::now() +
        std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(seconds));
    auto check = [&] {if (std::chrono::steady_clock::now() >= deadline) throw SearchExpired();};
    check();
    auto colors = source_types;
    std::sort(colors.begin(), colors.end());
    auto right_colors = target_types;
    std::sort(right_colors.begin(), right_colors.end());
    if (colors != right_colors) throw std::invalid_argument("fragment element multiplicities differ");
    colors.erase(std::unique(colors.begin(), colors.end()), colors.end());
    using Profiles = std::vector<std::vector<std::vector<double>>>;
    auto profiles = [&](const Matrix& weights, const std::vector<int>& types) {
        Profiles result(n, std::vector<std::vector<double>>(colors.size()));
        for (int i = 0; i < n; ++i) {
            check();
            for (int j = 0; j < n; ++j) if (i != j) {
                auto kind = std::lower_bound(colors.begin(), colors.end(), types[j]) - colors.begin();
                result[i][kind].push_back(weights[i][j]);
            }
            for (auto& values : result[i]) std::sort(values.begin(), values.end());
        }
        return result;
    };
    const auto left = profiles(source, source_types), right = profiles(target, target_types);
    Matrix limits(n, std::vector<double>(colors.size(), -1.));
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (i != j) {
        auto kind = std::lower_bound(colors.begin(), colors.end(), source_types[j]) - colors.begin();
        auto& value = limits[i][kind];
        if (value >= 0. && value != thresholds[i][j])
            throw std::invalid_argument("event threshold must be constant within each element pair");
        value = thresholds[i][j];
    }
    std::vector<int> fixed(n, -1), occupied(n, -1);
    for (const auto& [a, p] : anchors) {
        if (a < 0 || a >= n || p < 0 || p >= n || fixed[a] >= 0 || occupied[p] >= 0
                || source_types[a] != target_types[p])
            throw std::invalid_argument("invalid event anchor");
        fixed[a] = p; occupied[p] = a;
    }
    const int impossible = n * n + 1;
    std::vector<std::vector<int>> costs(n, std::vector<int>(n, impossible));
    for (int i = 0; i < n; ++i) {
        check();
        for (int p = 0; p < n; ++p) {
            if (source_types[i] != target_types[p] || (fixed[i] >= 0 && fixed[i] != p)
                    || (occupied[p] >= 0 && occupied[p] != i)) continue;
            int events = 0;
            for (size_t kind = 0; kind < colors.size(); ++kind) {
                const auto& a = left[i][kind];
                const auto& b = right[p][kind];
                size_t x = 0, y = 0, matched = 0;
                while (x < a.size() && y < b.size()) {
                    if (std::abs(a[x] - b[y]) <= limits[i][kind]) {++x; ++y; ++matched;}
                    else if (a[x] < b[y]) ++x;
                    else ++y;
                }
                events += static_cast<int>(a.size() - matched);
            }
            costs[i][p] = events;
        }
    }
    return (assignment_cost(costs, impossible, check) + 1) / 2;
}

class AssignmentSearch {
    struct Edge {
        int first, second;
        std::vector<std::pair<int, int>> pairs;
    };
    struct Frame {
        Domains domains;
        int atom = -1;
        std::vector<int> choices;
        size_t next = 0;
        bool entered = false;
    };
    int size;
    Matrix source, target, thresholds, ordering;
    std::vector<int> degrees;
    std::vector<Edge> edges;
    std::vector<Frame> stack;
    std::chrono::steady_clock::time_point deadline;

    void check() const {
        if (std::chrono::steady_clock::now() >= deadline) throw SearchExpired();
    }

    bool propagate(Domains& domains) const {
        bool changed = true;
        while (changed) {
            check();
            changed = false;
            // Accepted edge tables and injective occupancy retain their joint
            // constraints. Domain propagation alone never certifies a mapping.
            for (const auto& edge : edges) {
                std::vector<bool> left(size, false), right(size, false);
                for (const auto& [a, b] : edge.pairs) {
                    if (std::binary_search(domains[edge.first].begin(), domains[edge.first].end(), a)
                        && std::binary_search(domains[edge.second].begin(), domains[edge.second].end(), b)) {
                        left[a] = right[b] = true;
                    }
                }
                for (const auto& item : {std::make_pair(edge.first, &left), std::make_pair(edge.second, &right)}) {
                    auto& values = domains[item.first];
                    auto end = std::remove_if(values.begin(), values.end(), [&](int p) {return !(*item.second)[p];});
                    changed |= end != values.end();
                    values.erase(end, values.end());
                    if (values.empty()) return false;
                }
            }
            std::vector<bool> fixed(size, false);
            for (const auto& values : domains) {
                if (values.empty()) return false;
                if (values.size() == 1) {
                    if (fixed[values[0]]) return false;
                    fixed[values[0]] = true;
                }
            }
            for (auto& values : domains) if (values.size() > 1) {
                auto end = std::remove_if(values.begin(), values.end(), [&](int p) {return fixed[p];});
                changed |= end != values.end();
                values.erase(end, values.end());
                if (values.empty()) return false;
            }
        }
        return true;
    }

    void validate_domains(const Domains& domains) const {
        if (static_cast<int>(domains.size()) != size) throw std::invalid_argument("assignment domain count differs from fragment size");
        for (const auto& values : domains)
            if (!std::is_sorted(values.begin(), values.end())
                || std::adjacent_find(values.begin(), values.end()) != values.end()
                || (!values.empty() && (values.front() < 0 || values.back() >= size)))
                throw std::invalid_argument("assignment domains must be sorted distinct local atoms");
    }

public:
    size_t nodes = 0, improvements = 0;
    bool expired = false;

    AssignmentSearch(Matrix left, Matrix right, Matrix limits, Matrix priorities,
            Domains domains, std::vector<int> degree,
            const std::vector<std::tuple<int, int, std::vector<std::pair<int, int>>>>& constraints)
        : size(static_cast<int>(left.size())), source(std::move(left)), target(std::move(right)),
          thresholds(std::move(limits)), ordering(std::move(priorities)), degrees(std::move(degree)) {
        for (const auto* matrix : {&source, &target, &thresholds, &ordering}) {
            if (static_cast<int>(matrix->size()) != size) throw std::invalid_argument("event matrix sizes differ");
            for (const auto& values : *matrix)
                if (static_cast<int>(values.size()) != size) throw std::invalid_argument("event matrices must be square");
        }
        if (static_cast<int>(degrees.size()) != size) throw std::invalid_argument("degree count differs from fragment size");
        validate_domains(domains);
        for (const auto& [a, b, pairs] : constraints) {
            if (a < 0 || b < 0 || a >= size || b >= size || a == b)
                throw std::invalid_argument("invalid discovery edge");
            for (const auto& [p, q] : pairs)
                if (p < 0 || q < 0 || p >= size || q >= size || p == q)
                    throw std::invalid_argument("invalid discovery pair");
            edges.push_back({a, b, pairs});
        }
        stack.push_back({std::move(domains)});
        deadline = std::chrono::steady_clock::time_point::max();
    }

    int bound(const Domains& domains, int ceiling) const {
        validate_domains(domains);
        std::vector<int> fixed, free, targets;
        std::vector<bool> occupied(size, false);
        for (int a = 0; a < size; ++a) {
            if (domains[a].empty()) return ceiling;
            if (domains[a].size() == 1) {fixed.push_back(a); occupied[domains[a][0]] = true;}
            else free.push_back(a);
        }
        if (fixed.empty()) return 0;
        int cost = 0;
        for (size_t i = 0; i < fixed.size(); ++i) for (size_t j = i + 1; j < fixed.size(); ++j) {
            int a = fixed[i], b = fixed[j];
            cost += std::abs(source[a][b] - target[domains[a][0]][domains[b][0]]) > thresholds[a][b];
        }
        if (free.empty() || cost >= ceiling) return cost;
        for (int p = 0; p < size; ++p) if (!occupied[p]) targets.push_back(p);
        if (targets.size() != free.size()) return ceiling;
        const int impossible = size * size + 1;
        std::vector<std::vector<int>> costs(free.size(), std::vector<int>(targets.size(), impossible));
        for (size_t i = 0; i < free.size(); ++i) for (size_t j = 0; j < targets.size(); ++j) {
            const int a = free[i], p = targets[j];
            if (!std::binary_search(domains[a].begin(), domains[a].end(), p)) continue;
            int value = 0;
            for (int b : fixed) value += std::abs(source[a][b] - target[p][domains[b][0]]) > thresholds[a][b];
            costs[i][j] = value;
        }
        int remaining = assignment_cost(costs, impossible, [&] {check();});
        return remaining >= impossible ? ceiling : cost + remaining;
    }

    std::optional<std::vector<int>> next(int ceiling, double seconds) {
        if (!std::isfinite(seconds) || seconds < 0) throw std::invalid_argument("search allowance must be finite and nonnegative");
        deadline = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(seconds));
        try {
            while (!stack.empty()) {
                check();
                auto& frame = stack.back();
                if (!frame.entered) {
                    frame.entered = true;
                    ++nodes;
                    if (!propagate(frame.domains) || bound(frame.domains, ceiling) >= ceiling) {stack.pop_back(); continue;}
                    for (int a = 0; a < size; ++a) if (frame.domains[a].size() > 1
                        && (frame.atom < 0 || std::make_tuple(frame.domains[a].size(), -degrees[a], a)
                            < std::make_tuple(frame.domains[frame.atom].size(), -degrees[frame.atom], frame.atom))) frame.atom = a;
                    if (frame.atom < 0) {
                        std::vector<int> images;
                        for (const auto& values : frame.domains) images.push_back(values[0]);
                        stack.pop_back();
                        ++improvements;
                        return images;
                    }
                    frame.choices = frame.domains[frame.atom];
                    std::sort(frame.choices.begin(), frame.choices.end(), [&](int a, int b) {
                        return std::make_pair(ordering[frame.atom][a], a) < std::make_pair(ordering[frame.atom][b], b);
                    });
                }
                if (frame.next == frame.choices.size()) {stack.pop_back(); continue;}
                Domains child = frame.domains;
                child[frame.atom] = {frame.choices[frame.next++]};
                stack.push_back({std::move(child)});
            }
        } catch (const SearchExpired&) {expired = true;}
        return std::nullopt;
    }
};
} // namespace event_native
