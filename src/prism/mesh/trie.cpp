#include "trie.h"

#include <fmt/core.h>

#include <algorithm>

namespace prism::mesh {
FacesLookupTrie::FacesLookupTrie(std::size_t n_vertices) {
    _root_nodes.resize(n_vertices);
    for (std::size_t i = 0; i < n_vertices; ++i) {
        _root_nodes[i].v_id = i;
    }
}

void FacesLookupTrie::insert(const std::vector<std::size_t>& face_labels, std::size_t face_id) {
    auto root_node_id = face_labels[0];
    auto& root_node = _root_nodes[root_node_id];

    auto* current_node = &root_node;

    for (std::size_t i = 1; i < face_labels.size(); ++i) {
        auto v_id = face_labels[i];
        auto it = std::find_if(current_node->children.begin(),
                               current_node->children.end(),
                               [v_id](const auto& child) {
                                   return static_cast<VertexNode*>(child.get())->v_id == v_id;
                               });

        if (it == current_node->children.end()) {
            auto new_node = std::make_unique<VertexNode>();
            new_node->v_id = v_id;
            current_node->children.push_back(std::move(new_node));
            current_node = static_cast<VertexNode*>(current_node->children.back().get());
        }

        else {
            current_node = static_cast<VertexNode*>(it->get());
        }
    }

    auto tail_node = std::make_unique<TailNode>();
    tail_node->face_id = face_id;
    current_node->children.push_back(std::move(tail_node));
    current_node->is_end = true;
}

auto FacesLookupTrie::find(const std::vector<std::size_t>& face_labels) const
    -> std::optional<std::size_t> {
    auto root_node_id = face_labels[0];
    const auto& root_node = _root_nodes[root_node_id];

    const auto* current_node = &root_node;

    for (std::size_t i = 1; i < face_labels.size(); ++i) {
        auto v_id = face_labels[i];
        auto it = std::find_if(current_node->children.begin(),
                               current_node->children.end(),
                               [v_id](const auto& child) {
                                   return static_cast<VertexNode*>(child.get())->v_id == v_id;
                               });

        if (it == current_node->children.end()) {
            return std::nullopt;
        }

        else {
            current_node = static_cast<VertexNode*>(it->get());
        }
    }

    if (current_node->is_end) {
        return static_cast<TailNode*>(current_node->children[0].get())->face_id;
    }

    return std::nullopt;
}
} // namespace prism::mesh