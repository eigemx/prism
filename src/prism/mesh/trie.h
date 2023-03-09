#pragma once

#include <memory>
#include <optional>
#include <vector>

namespace prism::mesh {
class FacesLookupTrie {
  public:
    FacesLookupTrie(std::size_t n_vertices);

    void add(const std::vector<std::size_t>& face_labels, std::size_t face_id);
    auto contains(const std::vector<std::size_t>& face_labels) const
        -> std::optional<std::size_t>;

  private:
    struct Node {};

    struct VertexNode : Node {
        std::size_t v_id;
        std::vector<std::unique_ptr<Node>> children;
        bool is_end {false};
    };

    struct TailNode : Node {
        std::size_t face_id;
    };

    std::vector<VertexNode> _root_nodes;
};
} // namespace prism::mesh