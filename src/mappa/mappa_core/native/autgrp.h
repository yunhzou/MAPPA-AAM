// Shared native graph symmetry solver. No Python transport or chemistry policy.
#pragma once
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
extern "C" {
#include "nauty.h"
}
namespace autgrp_native {
// nauty's userautomproc carries no user-data pointer; nautywrap.c uses a
// static global (GRAPH_PTR) for the same reason.
inline thread_local std::vector<std::vector<int>>* g_sink = nullptr;

inline void store_generator(int /*count*/, int* perm, int* /*orbits*/,
                            int /*numorbits*/, int /*stabvertex*/, int n) {
    if (g_sink != nullptr) g_sink->emplace_back(perm, perm + n);
}

struct AutGraph {
    int n = 0;
    int m = 0;
    std::vector<setword> adjacency;
    std::vector<int> lab, ptn, orbits;
    std::vector<setword> workspace;
    optionblk options;
    statsblk stats;

    AutGraph(int vertices, const std::vector<std::pair<int, int>>& edges) : n(vertices) {
        if (n < 1) throw std::invalid_argument("AutGraph needs at least one vertex");
        m = (n + WORDSIZE - 1) / WORDSIZE;
        nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
        adjacency.assign((size_t)n * m, 0);
        for (const auto& e : edges) {
            check_vertex(e.first);
            check_vertex(e.second);
            ADDELEMENT(GRAPHROW(adjacency.data(), e.first, m), e.second);
            ADDELEMENT(GRAPHROW(adjacency.data(), e.second, m), e.first);
        }
        lab.assign(n, 0);
        ptn.assign(n, 0);
        orbits.assign(n, 0);
        workspace.assign((size_t)66 * m, 0);
        static DEFAULTOPTIONS_GRAPH(defaults);
        options = defaults;
        options.digraph = FALSE;
        options.getcanon = FALSE;
        options.defaultptn = TRUE;
        options.writeautoms = FALSE;
        options.cartesian = TRUE;
        options.linelength = 0;
        options.userautomproc = store_generator;
    }

    void check_vertex(int v) const {
        if (v < 0 || v >= n)
            throw std::invalid_argument("vertex " + std::to_string(v) +
                                        " conflicts with n_vertices=" + std::to_string(n));
    }

    // ``cells``: the final pynauty ``vertex_coloring`` (after its completion
    // and single-cell dropping rules), each set already converted with
    // ``list(cell)``.  Empty -> nauty's default single-cell partition.
    std::vector<std::vector<int>> generators(const std::vector<std::vector<int>>& cells) {
        if (cells.empty()) {
            options.defaultptn = TRUE;
        } else {
            options.defaultptn = FALSE;
            std::vector<char> seen((size_t)n, 0);
            int k = 0;
            for (const auto& cell : cells) {
                for (int v : cell) {
                    check_vertex(v);
                    if (seen[(size_t)v])
                        throw std::invalid_argument("vertex " + std::to_string(v) +
                                                    " appears in two cells");
                    seen[(size_t)v] = 1;
                    lab[(size_t)k] = v;
                    ptn[(size_t)k] = 1;
                    ++k;
                }
                if (k > 0) ptn[(size_t)k - 1] = 0;
            }
            if (k != n) throw std::invalid_argument("colouring does not cover every vertex");
        }
        std::vector<std::vector<int>> out;
        g_sink = &out;
        nauty(adjacency.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options,
              &stats, workspace.data(), (int)workspace.size(), m, n, nullptr);
        g_sink = nullptr;
        return out;
    }
};


}  // namespace autgrp_native
