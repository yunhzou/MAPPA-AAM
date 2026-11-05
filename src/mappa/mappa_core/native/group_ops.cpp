// Exact permutation bookkeeping. Search choices and chemical semantics live in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "set_orbit.h"

namespace py = pybind11;
using Images = std::vector<int>;
using Fragment = std::pair<int, Images>;
using Bond = std::pair<int, int>;
using Key = std::tuple<Images, Images, std::vector<Fragment>, std::vector<Bond>>;
using State = std::pair<Images, Images>;

static py::tuple atom_tuple(const Images& atoms) {
    py::tuple result(atoms.size());
    for (size_t i = 0; i < atoms.size(); ++i) result[i] = atoms[i];
    return result;
}

// Project completed states to interned fragment IDs in one call. This is only
// partition bookkeeping: it performs no matching, symmetry or event search.
static py::tuple fragment_compositions(py::sequence states, const Images& terminals,
                                       const Images& policies) {
    using AtomSets = std::pair<Images, Images>;
    std::map<std::pair<int, AtomSets>, int> fragments;
    std::map<std::pair<std::vector<Bond>, Images>, std::vector<Images>> partitions;
    py::list records, compositions;
    for (int terminal : terminals) {
        py::object state = states[terminal];
        const int policy = policies.at(py::cast<int>(state.attr("context")));
        std::vector<Bond> assignments, owners;
        Images source_atoms;
        for (py::handle item : py::iterable(state.attr("mapping"))) {
            const auto pair = py::cast<Bond>(item);
            if (pair.first < 0 || pair.second < 0)
                throw std::invalid_argument("completed atom indices must be nonnegative");
            assignments.push_back(pair);
            source_atoms.push_back(pair.first);
        }
        std::sort(source_atoms.begin(), source_atoms.end());
        Images mapping(source_atoms.empty() ? 0 : source_atoms.back()+1, -1);
        for (const auto& [atom, target] : assignments) mapping[atom] = target;
        for (py::handle item : py::iterable(state.attr("islands"))) {
            owners.push_back(py::cast<Bond>(item));
        }
        std::sort(owners.begin(), owners.end());
        const auto partition_key = std::make_pair(owners, source_atoms);
        auto found = partitions.find(partition_key);
        if (found == partitions.end()) {
            std::map<std::pair<bool, int>, Images> groups;
            auto owner = owners.begin();
            for (int atom : source_atoms) {
                while (owner != owners.end() && owner->first < atom) ++owner;
                const bool singleton = owner == owners.end() || owner->first != atom;
                groups[{singleton, singleton ? atom : owner->second}].push_back(atom);
            }
            std::vector<Images> partition;
            for (auto& [owner, atoms] : groups) partition.push_back(std::move(atoms));
            std::sort(partition.begin(), partition.end());
            found = partitions.emplace(partition_key, std::move(partition)).first;
        }
        Images ids;
        for (const Images& source : found->second) {
            Images target;
            for (int atom : source) target.push_back(mapping.at(atom));
            std::sort(target.begin(), target.end());
            const auto key = std::make_pair(policy, AtomSets{source, target});
            const auto inserted = fragments.try_emplace(key, static_cast<int>(fragments.size()));
            ids.push_back(inserted.first->second);
            if (inserted.second) {
                py::tuple witness(source.size());
                for (size_t i = 0; i < source.size(); ++i)
                    witness[i] = py::make_tuple(source[i], mapping.at(source[i]));
                records.append(py::make_tuple(policy, atom_tuple(source), atom_tuple(target), witness));
            }
        }
        compositions.append(py::make_tuple(terminal, atom_tuple(ids)));
    }
    return py::make_tuple(records, compositions);
}

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
    Images fixed;

    static std::shared_ptr<PairAction> from_permutation(const Images& images) {
        check_permutation(images);
        auto action = std::make_shared<PairAction>(images.size());
        action->fixed = images;
        action->fixed.push_back(images.size());
        return action;
    }

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
    struct Predecessor { size_t pair; int generator; };
    struct Projection {
        std::vector<size_t> pairs;
        std::vector<std::vector<Predecessor>> history;
        size_t visited = 0;
    };

    Projection project_states(Bond pair, bool trace) const {
        for (int atom : {pair.first, pair.second})
            if (atom < -1 || atom >= degree) throw std::invalid_argument("endpoint outside frame");
        Projection result;
        const size_t width = static_cast<size_t>(degree)+1;
        std::vector<unsigned char> seen(width*width, false);
        result.pairs = {static_cast<size_t>(pair.first < 0 ? degree : pair.first)*width
                       + (pair.second < 0 ? degree : pair.second)};
        for (const auto& action : actions) {
            std::fill(seen.begin(), seen.end(), false);
            std::vector<size_t> next;
            if (trace) result.history.emplace_back(width*width);
            auto add = [&](int a, int b, size_t previous, int generator = -1) {
                const size_t key = static_cast<size_t>(a)*width+b;
                if (!seen[key]) {
                    seen[key] = true;
                    next.push_back(key);
                    if (trace) result.history.back()[key] = {previous, generator};
                }
            };
            if (!action->fixed.empty()) {
                for (size_t key : result.pairs)
                    add(action->fixed[key/width], action->fixed[key%width], key);
            } else if (!action->pool.empty()) {
                bool both_expanded[2] = {false, false};
                std::vector<unsigned char> rows(width, false), columns(width, false);
                for (size_t key : result.pairs) {
                    const int a = key/width, b = key%width;
                    const bool left = action->member[a], right = action->member[b];
                    if (left && right) {
                        const bool diagonal = a == b;
                        if (both_expanded[diagonal]) continue;
                        both_expanded[diagonal] = true;
                        for (int x : action->pool) for (int y : action->pool)
                            if ((x == y) == diagonal) add(x, y, key);
                    } else if (left) {
                        if (columns[b]) continue;
                        columns[b] = true;
                        for (int x : action->pool) add(x, b, key);
                    } else if (right) {
                        if (rows[a]) continue;
                        rows[a] = true;
                        for (int y : action->pool) add(a, y, key);
                    } else {
                        add(a, b, key);
                    }
                }
            } else {
                for (size_t key : result.pairs) add(key/width, key%width, key);
                for (size_t i = 0; i < next.size(); ++i) {
                    const size_t key = next[i];
                    const int a = key/width, b = key%width;
                    for (size_t j = 0; j < action->generators.size(); ++j) {
                        const auto& g = action->generators[j];
                        add(g[a], g[b], key, j);
                    }
                }
            }
            result.visited += next.size();
            result.pairs.swap(next);
        }
        return result;
    }

    std::vector<Images> realize(size_t key, const Projection& projection) const {
        const size_t width = static_cast<size_t>(degree)+1;
        std::vector<Images> choices(actions.size(), Images(degree));
        for (size_t layer = actions.size(); layer-- > 0;) {
            const auto& action = actions[layer];
            auto& choice = choices[layer];
            std::iota(choice.begin(), choice.end(), 0);
            const auto& history = projection.history[layer];
            if (!action->fixed.empty()) {
                choice.assign(action->fixed.begin(), action->fixed.begin()+degree);
                key = history[key].pair;
            } else if (!action->pool.empty()) {
                const size_t previous = history[key].pair;
                std::map<int, int> selected;
                for (const auto& [a, b] : {Bond(previous/width, key/width), Bond(previous%width, key%width)})
                    if (action->member[a]) selected[a] = b;
                std::set<int> unused(action->pool.begin(), action->pool.end());
                for (const auto& [a, b] : selected) { choice[a] = b; unused.erase(b); }
                for (int a : action->pool) if (!selected.count(a)) {
                    choice[a] = *unused.begin();
                    unused.erase(unused.begin());
                }
                key = previous;
            } else {
                Images word;
                while (history[key].generator >= 0) {
                    word.push_back(history[key].generator);
                    key = history[key].pair;
                }
                key = history[key].pair;
                for (auto i = word.rbegin(); i != word.rend(); ++i)
                    for (int& image : choice) image = action->generators[*i][image];
            }
        }
        return choices;
    }
public:
    EndpointPairProgram(int n, std::vector<std::shared_ptr<PairAction>> program)
        : degree(n), actions(std::move(program)) {
        if (degree < 0) throw std::invalid_argument("negative permutation degree");
        for (const auto& action : actions)
            if (action->degree != degree) throw std::invalid_argument("action frame differs");
    }

    py::tuple project(Bond pair) const {
        std::vector<Bond> pairs;
        {
            py::gil_scoped_release release;
            const size_t width = static_cast<size_t>(degree)+1;
            const auto projection = project_states(pair, false);
            pairs.reserve(projection.pairs.size());
            for (size_t key : projection.pairs) {
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

    py::tuple minimum(Bond pair, py::array_t<int, py::array::c_style | py::array::forcecast> costs) const {
        if (pair.first < 0 || pair.second < 0)
            throw std::invalid_argument("minimum requires two mapped endpoints");
        if (costs.ndim() != 2 || costs.shape(0) != degree || costs.shape(1) != degree)
            throw std::invalid_argument("pair costs differ from the action frame");
        int best;
        size_t visited;
        std::vector<Images> choices;
        {
            py::gil_scoped_release release;
            const auto projection = project_states(pair, true);
            const size_t width = static_cast<size_t>(degree)+1;
            const int* values = costs.data();
            auto value = [&](size_t key) { return values[(key/width)*degree + key%width]; };
            const size_t selected = *std::min_element(projection.pairs.begin(), projection.pairs.end(),
                [&](size_t a, size_t b) { return value(a) < value(b); });
            best = value(selected);
            visited = projection.visited;
            choices = realize(selected, projection);
        }
        return py::make_tuple(best, choices, visited);
    }
};

PYBIND11_MODULE(_group_ops, m) {
    m.def("fragment_compositions", &fragment_compositions);
    py::class_<permutation_native::SetOrbit>(m, "SetOrbit")
        .def(py::init<int, const permutation_native::Generators&>())
        .def("contains", &permutation_native::SetOrbit::contains)
        .def("stabilizer", &permutation_native::SetOrbit::stabilizer);
    py::class_<PairAction, std::shared_ptr<PairAction>>(m, "PairAction")
        .def_static("pool", &PairAction::from_pool)
        .def_static("permutation", &PairAction::from_permutation)
        .def_static("group", &PairAction::from_group);
    py::class_<EndpointPairProgram>(m, "EndpointPairProgram")
        .def(py::init<int, std::vector<std::shared_ptr<PairAction>>>())
        .def("project", &EndpointPairProgram::project)
        .def("minimum", &EndpointPairProgram::minimum);

    m.def("occupation_orbit", &occupation_orbit,
          py::arg("witness"), py::arg("degree"), py::arg("stages"),
          py::arg("attachments"), py::arg("fragments"), py::arg("bonds"),
          py::arg("limit"), py::arg("observed_atoms") = py::none(),
          py::call_guard<py::gil_scoped_release>());
    m.def("project_generators", &project_generators);
    m.def("conjugate_generators", &conjugate_generators);
    py::register_exception<std::length_error>(m, "OccupationLimitExceeded");
}
