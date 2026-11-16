// Shared search ancestry and ordered fragment/placement programs. A join retains its
// parent links; this projection never expands all paths through the graph.
#pragma once
#include <algorithm>
#include <map>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace reduction_native {
struct HistoryEdge {int id, source, target, decision, key;};
using Parents = std::vector<std::pair<int, int>>;
using Choices = std::vector<std::pair<int, int>>;
struct HistoryProjection {
    std::vector<Parents> parents{{}};
    std::vector<std::pair<int, int>> programs{{-1, -1}};
    std::vector<std::tuple<int, int, int>> endings;
};

inline HistoryProjection collect_histories(int states, const std::vector<int>& roots,
        const std::vector<int>& terminals, const std::vector<HistoryEdge>& edges,
        const std::optional<std::vector<std::pair<int, std::vector<int>>>>& selected) {
    HistoryProjection result;
    std::map<std::pair<int, int>, int> program_ids;
    auto append = [&](int program, const HistoryEdge& edge) {
        if (edge.decision < 0) return program;
        auto key = std::make_pair(edge.key, program);
        auto [entry, added] = program_ids.emplace(key, static_cast<int>(result.programs.size()));
        // Keep one original reachable decision, without making its witness
        // another program on the same atom sets and placement context.
        if (added) result.programs.emplace_back(edge.decision, program);
        return entry->second;
    };
    auto add = [&](Parents parents) {
        int node = static_cast<int>(result.parents.size());
        result.parents.push_back(std::move(parents));
        return node;
    };
    auto valid_state = [&](int state) {
        if (state < 0 || state >= states) throw std::invalid_argument("history state outside graph");
    };
    if (states < 0) throw std::invalid_argument("negative history state count");
    for (int root : roots) valid_state(root);
    for (int terminal : terminals) valid_state(terminal);
    for (size_t i = 0; i < edges.size(); ++i) {
        valid_state(edges[i].source); valid_state(edges[i].target);
        if (edges[i].id != static_cast<int>(i)) throw std::invalid_argument("history edge IDs must index transitions");
    }
    if (selected) {
        std::map<std::pair<int, int>, int> links;
        std::vector<int> terminal_order;
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> endings;
        for (const auto& [terminal, path] : *selected) {
            valid_state(terminal);
            int program = 0, node = 0;
            for (int eid : path) {
                if (eid < 0 || eid >= static_cast<int>(edges.size())) throw std::invalid_argument("history path edge outside graph");
                program = append(program, edges[eid]);
                auto key = std::make_pair(eid, node);
                auto [entry, added] = links.emplace(key, static_cast<int>(result.parents.size()));
                if (added) add({key});
                node = entry->second;
            }
            auto [entry, added] = endings.try_emplace(terminal);
            if (added) terminal_order.push_back(terminal);
            auto& choices = entry->second;
            auto choice = std::find_if(choices.begin(), choices.end(), [&](const auto& item) {return item.first == program;});
            if (choice == choices.end()) {choices.emplace_back(program, std::vector<int>{}); choice = std::prev(choices.end());}
            choice->second.push_back(node);
        }
        for (int terminal : terminal_order) for (const auto& [program, nodes] : endings[terminal]) {
            Parents parents;
            for (int node : nodes) parents.emplace_back(-1, node);
            result.endings.emplace_back(terminal, program, add(std::move(parents)));
        }
    } else {
        std::vector<std::vector<int>> incoming(states);
        for (const auto& edge : edges) incoming[edge.target].push_back(edge.id);
        std::vector<bool> root(states, false), ready(states, false);
        for (int state : roots) root[state] = true;
        std::vector<int> order;
        for (int end : terminals) {
            std::vector<std::pair<int, bool>> stack{{end, false}};
            while (!stack.empty()) {
                auto [state, complete] = stack.back(); stack.pop_back();
                if (ready[state]) continue;
                if (complete || root[state]) {ready[state] = true; order.push_back(state);}
                else {
                    stack.emplace_back(state, true);
                    for (auto edge = incoming[state].rbegin(); edge != incoming[state].rend(); ++edge)
                        stack.emplace_back(edges[*edge].source, false);
                }
            }
        }
        std::vector<Choices> choices(states);
        for (int state : order) {
            if (root[state]) {choices[state].emplace_back(0, 0); continue;}
            std::map<int, size_t> positions;
            std::vector<std::pair<int, Parents>> parents;
            for (int eid : incoming[state]) {
                const auto& edge = edges[eid];
                for (const auto& [program, node] : choices[edge.source]) {
                    int next = append(program, edge);
                    auto [entry, added] = positions.emplace(next, parents.size());
                    if (added) parents.emplace_back(next, Parents{});
                    parents[entry->second].second.emplace_back(edge.id, node);
                }
            }
            for (auto& [program, links] : parents) choices[state].emplace_back(program, add(std::move(links)));
        }
        for (int terminal : terminals) for (const auto& [program, node] : choices[terminal])
            result.endings.emplace_back(terminal, program, node);
    }
    return result;
}
} // namespace reduction_native
