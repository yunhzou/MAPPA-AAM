// Conservative domains for ordered actions. These bounds accelerate queries;
// they never replace the actions' correlated scientific relationship.
#pragma once
#include "set_orbit.h"
#include <array>
#include <list>
#include <mutex>

namespace permutation_native {
using Domain = std::shared_ptr<const Permutation>;

class DomainAction {
public:
    std::map<int, Domain> rows;
    explicit DomainAction(const std::map<int, Permutation>& values) {
        // Share equal rows without closing overlapping conservative domains.
        std::map<Permutation, Domain> shared;
        for (const auto& [atom, domain] : values) {
            if (atom < 0 || domain.empty() || domain.front() < 0
                || !std::is_sorted(domain.begin(), domain.end())
                || std::adjacent_find(domain.begin(), domain.end()) != domain.end()
                || !std::binary_search(domain.begin(), domain.end(), atom))
                throw std::invalid_argument("action bounds must contain their original atom and distinct sorted images");
            auto [entry, added] = shared.try_emplace(domain);
            if (added) entry->second = std::make_shared<const Permutation>(domain);
            rows.emplace(atom, entry->second);
        }
    }
};

struct DomainBinding {
    Permutation atoms;
    std::vector<Domain> domains;
    std::vector<Permutation> dependencies;
};

class DomainPrograms {
    using Key = std::pair<std::shared_ptr<DomainAction>, Permutation>;
    struct Entry {Domain domain; size_t bytes; std::list<Key>::iterator position;};
    int degree;
    std::vector<Domain> fixed;
    std::map<Key, Entry> cache;
    std::list<Key> recent;
    size_t retained_bytes = 0, hits = 0, misses = 0;
    std::mutex guard;

    Domain advance(const std::shared_ptr<DomainAction>& action, const Domain& domain) {
        if (domain->size() == 1) {
            auto found = action->rows.find(domain->front());
            return found == action->rows.end() ? domain : found->second;
        }
        Key key{action, *domain};
        auto found = cache.find(key);
        if (found != cache.end()) {
            ++hits;
            recent.splice(recent.end(), recent, found->second.position);
            return found->second.domain;
        }
        ++misses;
        auto values = *domain;
        for (int p : *domain) {
            auto row = action->rows.find(p);
            if (row != action->rows.end()) values.insert(values.end(), row->second->begin(), row->second->end());
        }
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
        auto result = std::make_shared<const Permutation>(std::move(values));
        // Account both map/list keys, result storage and tree/list links.
        const size_t bytes = sizeof(Entry) + 2 * sizeof(Key) + 6 * sizeof(void*)
            + (2 * domain->size() + result->capacity()) * sizeof(int);
        recent.push_back(key);
        cache.emplace(std::move(key), Entry{result, bytes, std::prev(recent.end())});
        retained_bytes += bytes;
        while (cache.size() > 4096 || retained_bytes > 16 * 1024 * 1024) {
            auto oldest = cache.find(recent.front());
            retained_bytes -= oldest->second.bytes;
            cache.erase(oldest);
            recent.pop_front();
        }
        return result;
    }
public:
    explicit DomainPrograms(int n): degree(n) {
        if (n < 0) throw std::invalid_argument("negative action frame");
        for (int p = 0; p < n; ++p) fixed.push_back(std::make_shared<const Permutation>(Permutation{p}));
    }

    DomainBinding bind(const std::vector<std::shared_ptr<DomainAction>>& actions,
                       const std::vector<std::pair<int, int>>& mapping) {
        std::lock_guard<std::mutex> lock(guard);
        DomainBinding result;
        result.dependencies.resize(mapping.size());
        std::vector<std::set<size_t>> reaching(degree);
        for (const auto& [atom, image] : mapping) {
            if (atom < 0 || image < 0 || image >= degree) throw std::invalid_argument("mapping outside action frame");
            reaching[image].insert(result.atoms.size());
            result.atoms.push_back(atom);
            result.domains.push_back(fixed[image]);
        }
        for (size_t position = 0; position < actions.size(); ++position) {
            const auto& action = actions[position];
            std::set<size_t> affected;
            for (const auto& [image, domain] : action->rows) {
                if (image >= degree || domain->back() >= degree) throw std::invalid_argument("action outside its frame");
                affected.insert(reaching[image].begin(), reaching[image].end());
            }
            // Freeze incidence before advancing: a non-transitive relationship
            // is applied once at this position, not repeatedly to its new images.
            for (size_t row : affected) {
                auto before = result.domains[row];
                auto after = advance(action, before);
                result.domains[row] = after;
                result.dependencies[row].push_back(static_cast<int>(position));
                if (*before != *after) for (int p : *after)
                    if (!std::binary_search(before->begin(), before->end(), p)) reaching[p].insert(row);
            }
        }
        return result;
    }

    std::array<size_t, 4> statistics() {
        std::lock_guard<std::mutex> lock(guard);
        return {cache.size(), retained_bytes, hits, misses};
    }
};
} // namespace permutation_native
