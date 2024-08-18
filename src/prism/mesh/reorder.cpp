#include "reorder.h"

#include <algorithm>
#include <map>
#include <queue>

namespace prism::mesh {
CuthillMckee::CuthillMckee(PMesh& mesh) noexcept : _mesh(mesh) {
    auto n_cells = mesh.cells().size();
    _nodes.reserve(n_cells);

    // fill the adjacency graph
    for (const auto& cell : mesh.cells()) {
        auto cell_id = cell.id();
        auto node = Node(cell_id, 0);

        for (const auto& face_id : cell.facesIds()) {
            const auto& face = mesh.faces()[face_id];

            if (face.isBoundary()) {
                continue;
            }

            bool is_owner = face.owner() == cell_id;
            auto neighbor_cell_id = is_owner ? face.neighbor().value() : face.owner();
            node.neighbors().push_back(neighbor_cell_id);
            node.degree()++;
        }

        _nodes.push_back(node);
    }

    // sort nodes by degree, ascending
    std::sort(_nodes.begin(), _nodes.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.degree() < rhs.degree();
    });
}

auto CuthillMckee::permute() -> std::vector<std::size_t> {
    auto n = _mesh.cells().size();

    // initialize the permutation vector. This is the output of the algorithm
    std::vector<std::size_t> permutations;

    // The queue of nodes to visit
    std::deque<Node> queue;

    // initialize the set of not-visited nodes
    std::deque<Node> not_visited;
    for (const auto& node : _nodes) {
        not_visited.push_back(node);
    }

    while (!not_visited.empty()) {
        // add first node (that is not visited) to the queue,
        // and remove it from the set of not-visited nodes
        queue.push_back(not_visited.front());
        not_visited.pop_front();


        while (!queue.empty()) {
            // vector of neighbors of the current node, sorted by degree
            // and not visited yet
            std::vector<Node> not_visited_neighbors;

            for (const auto& neighbor_id : queue.front().neighbors()) {
                auto it = std::find_if(
                    not_visited.begin(), not_visited.end(), [neighbor_id](const auto& node) {
                        return node.id() == neighbor_id;
                    });

                if (it != not_visited.end()) {
                    not_visited_neighbors.push_back(*it);

                    // remove the neighbor from the set of not-visited nodes
                    not_visited.erase(it);
                }
            }

            // sort neighbors by degree, ascending
            std::sort(
                not_visited_neighbors.begin(),
                not_visited_neighbors.end(),
                [](const auto& lhs, const auto& rhs) { return lhs.degree() < rhs.degree(); });

            queue.insert(queue.end(), not_visited_neighbors.begin(), not_visited_neighbors.end());
            permutations.push_back(queue.front().id());

            queue.pop_front();
        }
    }

    return permutations;
}

void CuthillMckee::reorder(bool reverse) {
    auto perms = permute();
    auto n_cells = perms.size();

    // reverse CutHill-McKee permutation if needed
    if (reverse) {
        std::reverse(perms.begin(), perms.end());
    }

    std::map<std::size_t, std::size_t> old_to_new;
    for (auto i = 0; i < n_cells; ++i) {
        old_to_new[perms[i]] = i;
    }

    // update order of cells vector using the permutation map
    auto& cells = _mesh.cells();

    std::vector<Cell> new_cells = cells;

    // TODO: replace this with std::transform
    for (auto& cell : cells) {
        auto new_id = old_to_new[cell.id()];
        cell.id() = new_id;
        new_cells[new_id] = cell;
    }

    // swap the old cells vector with the new one
    cells.swap(new_cells);

    // update face owners and neighbors
    for (auto& face : _mesh.faces()) {
        face.setOwner(old_to_new[face.owner()]);

        if (face.isInterior()) {
            face.setNeighbor(old_to_new[face.neighbor().value()]);
        }
    }
}
} // namespace prism::mesh