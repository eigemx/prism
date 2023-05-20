#include "trie.h"

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
    auto face_labels_size = face_labels.size();
    auto& root_node = _root_nodes[root_node_id];

    auto* current_node = &root_node;

    for (std::size_t i = 1; i < face_labels_size; ++i) {
        auto v_id = face_labels[i];
        auto it = std::find_if(current_node->children.begin(),
                               current_node->children.end(),
                               [v_id](const auto& child) {
                                   return std::get_if<VertexNode>(&child) != nullptr &&
                                          std::get<VertexNode>(child).v_id == v_id;
                               });

        if (it == current_node->children.end()) {
            current_node->children.emplace_back(VertexNode {v_id, {}, false});
            current_node = &std::get<VertexNode>(current_node->children.back());
        }

        else {
            current_node = &std::get<VertexNode>(*it);
        }
    }

    auto tail_node = TailNode {face_id};
    current_node->children.emplace_back(tail_node);
    current_node->is_end = true;
}

auto FacesLookupTrie::find(const std::vector<std::size_t>& face_labels) const
    -> std::optional<std::size_t> {
    auto root_node_id = face_labels[0];
    auto face_labels_size = face_labels.size();
    const auto& root_node = _root_nodes[root_node_id];

    const auto* current_node = &root_node;

    for (std::size_t i = 1; i < face_labels_size; ++i) {
        auto v_id = face_labels[i];
        auto it = std::find_if(current_node->children.begin(),
                               current_node->children.end(),
                               [v_id](const auto& child) {
                                   return std::get_if<VertexNode>(&child) != nullptr &&
                                          std::get<VertexNode>(child).v_id == v_id;
                               });

        if (it == current_node->children.end()) {
            return std::nullopt;
        }

        current_node = &std::get<VertexNode>(*it);
    }

    if (current_node->is_end) {
        return std::get<TailNode>(current_node->children[0]).face_id;
    }

    return std::nullopt;
}
} // namespace prism::mesh