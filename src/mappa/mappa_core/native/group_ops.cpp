// Exact permutation bookkeeping. Search choices and chemical semantics live in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace py = pybind11;
using Images = std::vector<int>;
using Fragment = std::pair<int, Images>;
using Bond = std::pair<int, int>;
using Key = std::tuple<Images, Images, std::vector<Fragment>, std::vector<Bond>>;
using State = std::pair<Images, Images>;

// Only a closed, permanently unobservable set can be collapsed. Include every
// point that any recorded generator can transport to an observed point. This
// over-approximates chronological reachability and is safe for every stage.
static Images observation_labels(int degree,
        const std::vector<std::vector<Images>>& stages,
        const std::optional<Images>& observed) {
    Images labels(degree);
    std::iota(labels.begin(), labels.end(), 0);
    if (!observed) return labels;
    Images parent = labels;
    auto root = [&](int a) {
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };
    for (const auto& stage : stages) for (const auto& g : stage)
        for (int a = 0; a < degree; ++a) parent[root(a)] = root(g[a]);
    std::vector<bool> visible(degree, false);
    for (int a : *observed) {
        if (a < 0 || a >= degree) throw std::invalid_argument("observed atom outside frame");
        visible[root(a)] = true;
    }
    for (int a = 0; a < degree; ++a)
        if (!visible[root(a)]) labels[a] = -1;
    return labels;
}

static void check_permutation(const Images& p) {
    Images ordered = p;
    std::sort(ordered.begin(), ordered.end());
    for (size_t i = 0; i < ordered.size(); ++i)
        if (ordered[i] != static_cast<int>(i))
            throw std::invalid_argument("expected a permutation of its complete frame");
}

static Images select_sorted(const Images& images, const Images& positions) {
    Images selected;
    selected.reserve(positions.size());
    for (int position : positions) selected.push_back(images[position]);
    std::sort(selected.begin(), selected.end());
    return selected;
}

static Key occupation_key(const Images& images, const Images& attachments,
                          const std::vector<Fragment>& fragments,
                          const std::vector<Bond>& bonds) {
    Images covered = images;
    std::sort(covered.begin(), covered.end());
    std::vector<Fragment> parts;
    parts.reserve(fragments.size());
    for (const auto& [label, positions] : fragments)
        parts.emplace_back(label, select_sorted(images, positions));
    std::sort(parts.begin(), parts.end());
    std::vector<Bond> edges;
    edges.reserve(bonds.size());
    for (const auto& [a, b] : bonds)
        edges.push_back(std::minmax(images[a], images[b]));
    std::sort(edges.begin(), edges.end());
    return {std::move(covered), select_sorted(images, attachments), std::move(parts), std::move(edges)};
}

// Same chronological stage closure and first witness as the Python algorithm.
// This enumerates distinct fragment relations, never the complete bijection group.
static std::vector<State> occupation_orbit(
        const Images& witness, int degree, const std::vector<std::vector<Images>>& stages,
        const Images& attachments, const std::vector<Fragment>& fragments,
        const std::vector<Bond>& bonds, int limit,
        const std::optional<Images>& observed = std::nullopt) {
    if (degree < 0) throw std::invalid_argument("negative permutation degree");
    for (int image : witness)
        if (image < 0 || image >= degree) throw std::invalid_argument("image outside frame");
    auto check_position = [&](int i) {
        if (i < 0 || static_cast<size_t>(i) >= witness.size())
            throw std::invalid_argument("source position outside witness");
    };
    for (int i : attachments) check_position(i);
    for (const auto& fragment : fragments) for (int i : fragment.second) check_position(i);
    for (const auto& bond : bonds) { check_position(bond.first); check_position(bond.second); }
    for (const auto& stage : stages) for (const auto& g : stage) {
        if (g.size() != static_cast<size_t>(degree)) throw std::invalid_argument("generator frame mismatch");
        check_permutation(g);
    }
    Images identity(degree);
    std::iota(identity.begin(), identity.end(), 0);
    const Images labels = observation_labels(degree, stages, observed);
    auto projected = [&](const Images& images) {
        Images value;
        value.reserve(images.size());
        for (int image : images) value.push_back(labels[image]);
        return value;
    };
    std::vector<State> states{{witness, identity}};
    std::map<Key, size_t> seen;
    seen.emplace(occupation_key(projected(witness), attachments, fragments, bonds), 0);
    std::map<Images, size_t> evaluated;
    for (const auto& stage : stages) {
        std::vector<size_t> already;
        for (const auto& g : stage) already.push_back(evaluated[g]);
        for (size_t current = 0; current < states.size(); ++current) {
            // Appending states may reallocate; do not hold references into states.
            const State state = states[current];
            const Images state_projection = projected(state.first);
            for (size_t gi = 0; gi < stage.size(); ++gi) {
                if (current < already[gi]) continue;
                const auto& g = stage[gi];
                Images moved;
                moved.reserve(witness.size());
                for (int image : state.first) moved.push_back(g[image]);
                Images moved_projection = projected(moved);
                if (moved_projection == state_projection) continue;
                Key key = occupation_key(moved_projection, attachments, fragments, bonds);
                if (seen.find(key) != seen.end()) continue;
                if (limit >= 0 && states.size() >= static_cast<size_t>(limit))
                    throw std::length_error("fragment occupation limit exceeded");
                Images action;
                action.reserve(degree);
                for (int image : state.second) action.push_back(g[image]);
                seen.emplace(std::move(key), states.size());
                states.emplace_back(std::move(moved), std::move(action));
            }
        }
        for (const auto& g : stage) evaluated[g] = states.size();
    }
    // Preserve discovery order; the Python caller merges families and sorts final keys.
    return states;
}

static py::tuple project_generators(const std::vector<Images>& raw,
                                    const Images& atom_labels, const Images& vertex_indices) {
    if (atom_labels.size() != vertex_indices.size()) throw std::invalid_argument("atom index size mismatch");
    const int degree = atom_labels.empty() ? 0 : *std::max_element(atom_labels.begin(), atom_labels.end()) + 1;
    std::map<int, int> atoms;
    for (size_t i = 0; i < atom_labels.size(); ++i) atoms.emplace(vertex_indices[i], atom_labels[i]);
    Images identity(degree);
    std::iota(identity.begin(), identity.end(), 0);
    std::set<Images> seen;
    std::vector<Images> output;
    {
        py::gil_scoped_release release;
        for (const auto& g : raw) {
            Images p = identity;
            for (size_t i = 0; i < atom_labels.size(); ++i) {
                const int vertex = vertex_indices[i];
                if (vertex < 0 || static_cast<size_t>(vertex) >= g.size())
                    throw std::invalid_argument("atom vertex outside generator frame");
                auto image = atoms.find(g[vertex]);
                if (image == atoms.end()) throw std::invalid_argument("candidate automorphism mixed atom/edge vertices");
                p.at(atom_labels[i]) = image->second;
            }
            if (p != identity && seen.insert(p).second) output.push_back(std::move(p));
        }
    }
    py::tuple result(output.size());
    for (size_t i = 0; i < output.size(); ++i) result[i] = py::tuple(py::cast(output[i]));
    return result;
}

static py::tuple conjugate_generators(const std::vector<Images>& generators, const Images& action) {
    check_permutation(action);
    std::vector<Images> output;
    {
        py::gil_scoped_release release;
        for (const auto& g : generators) {
            check_permutation(g);
            size_t degree = std::max(g.size(), action.size());
            Images moved(degree);
            auto image = [&](size_t i) { return i < action.size() ? action[i] : static_cast<int>(i); };
            for (size_t i = 0; i < degree; ++i)
                moved[image(i)] = image(i < g.size() ? g[i] : i);
            output.push_back(std::move(moved));
        }
    }
    py::tuple result(output.size());
    for (size_t i = 0; i < output.size(); ++i) result[i] = py::tuple(py::cast(output[i]));
    return result;
}

// A compiled local action is shared by all branch programs that use it.
// Projection follows action order, retaining the two endpoints jointly. It
// visits at most (degree+1)^2 pair states, never complete permutation elements.
class PairAction {
public:
    int degree;
    Images pool;
    std::vector<unsigned char> member;
    std::vector<Images> generators;

    static std::shared_ptr<PairAction> from_pool(int degree, const Images& atoms) {
        auto action = std::make_shared<PairAction>(degree);
        for (int atom : atoms) {
            if (atom < 0 || atom >= degree || action->member[atom])
                throw std::invalid_argument("pool must contain distinct atoms in its frame");
            action->member[atom] = true;
        }
        action->pool = atoms;
        return action;
    }

    static std::shared_ptr<PairAction> from_group(int degree, const std::vector<Images>& raw) {
        auto action = std::make_shared<PairAction>(degree);
        std::set<Images> seen;
        for (const auto& permutation : raw) {
            if (permutation.size() != static_cast<size_t>(degree))
                throw std::invalid_argument("generator degree differs from its frame");
            check_permutation(permutation);
            bool moved = false;
            for (int i = 0; i < degree; ++i) moved |= permutation[i] != i;
            if (!moved || !seen.insert(permutation).second) continue;
            action->generators.push_back(permutation);
            action->generators.back().push_back(degree); // Unmatched endpoint is fixed.
        }
        return action;
    }

    explicit PairAction(int n) : degree(n), member(n >= 0 ? n+1 : 0, false) {
        if (n < 0) throw std::invalid_argument("negative permutation degree");
    }
};

class EndpointPairProgram {
    int degree;
    std::vector<std::shared_ptr<PairAction>> actions;
public:
    EndpointPairProgram(int n, std::vector<std::shared_ptr<PairAction>> program)
        : degree(n), actions(std::move(program)) {
        if (degree < 0) throw std::invalid_argument("negative permutation degree");
        for (const auto& action : actions)
            if (action->degree != degree) throw std::invalid_argument("action frame differs");
    }

    py::tuple project(Bond pair) const {
        for (int atom : {pair.first, pair.second})
            if (atom < -1 || atom >= degree) throw std::invalid_argument("endpoint outside frame");
        std::vector<Bond> pairs;
        {
            py::gil_scoped_release release;
            const size_t width = static_cast<size_t>(degree)+1;
            std::vector<unsigned char> seen(width*width, false);
            std::vector<size_t> states{static_cast<size_t>(pair.first < 0 ? degree : pair.first)*width
                                      + (pair.second < 0 ? degree : pair.second)};
            for (const auto& action : actions) {
                if (action->pool.empty() && action->generators.empty()) continue;
                std::fill(seen.begin(), seen.end(), false);
                std::vector<size_t> next;
                auto add = [&](int a, int b) {
                    const size_t key = static_cast<size_t>(a)*width+b;
                    if (!seen[key]) { seen[key] = true; next.push_back(key); }
                };
                if (!action->pool.empty()) {
                    bool both_expanded[2] = {false, false};
                    std::vector<unsigned char> rows(width, false), columns(width, false);
                    for (size_t key : states) {
                        const int a = key/width, b = key%width;
                        const bool left = action->member[a], right = action->member[b];
                        if (left && right) {
                            const bool diagonal = a == b;
                            if (both_expanded[diagonal]) continue;
                            both_expanded[diagonal] = true;
                            for (int x : action->pool) for (int y : action->pool)
                                if ((x == y) == diagonal) add(x, y);
                        } else if (left) {
                            if (columns[b]) continue;
                            columns[b] = true;
                            for (int x : action->pool) add(x, b);
                        } else if (right) {
                            if (rows[a]) continue;
                            rows[a] = true;
                            for (int y : action->pool) add(a, y);
                        } else {
                            add(a, b);
                        }
                    }
                } else {
                    for (size_t key : states) add(key/width, key%width);
                    for (size_t i = 0; i < next.size(); ++i) {
                        const int a = next[i]/width, b = next[i]%width;
                        for (const auto& g : action->generators) add(g[a], g[b]);
                    }
                }
                states.swap(next);
            }
            pairs.reserve(states.size());
            for (size_t key : states) {
                const int a = key/width, b = key%width;
                pairs.emplace_back(a == degree ? -1 : a, b == degree ? -1 : b);
            }
            std::sort(pairs.begin(), pairs.end());
        }
        py::tuple result(pairs.size());
        for (size_t i = 0; i < pairs.size(); ++i)
            result[i] = py::make_tuple(pairs[i].first, pairs[i].second);
        return result;
    }
};

PYBIND11_MODULE(_group_ops, m) {
    py::class_<PairAction, std::shared_ptr<PairAction>>(m, "PairAction")
        .def_static("pool", &PairAction::from_pool)
        .def_static("group", &PairAction::from_group);
    py::class_<EndpointPairProgram>(m, "EndpointPairProgram")
        .def(py::init<int, std::vector<std::shared_ptr<PairAction>>>())
        .def("project", &EndpointPairProgram::project);

    m.def("occupation_orbit", &occupation_orbit,
          py::arg("witness"), py::arg("degree"), py::arg("stages"),
          py::arg("attachments"), py::arg("fragments"), py::arg("bonds"),
          py::arg("limit"), py::arg("observed_atoms") = py::none(),
          py::call_guard<py::gil_scoped_release>());
    m.def("project_generators", &project_generators);
    m.def("conjugate_generators", &conjugate_generators);
    py::register_exception<std::length_error>(m, "OccupationLimitExceeded");
}
