#pragma once

#include <optional>
#include <variant>
#include <vector>

namespace prism::mesh {
class FacesLookupTrie {
  public:
    FacesLookupTrie(size_t n_vertices);

    void insert(const std::vector<size_t>& face_labels, size_t face_id);
    auto find(const std::vector<size_t>& face_labels) const -> std::optional<size_t>;

  private:
    // forward declarations
    struct VertexNode;
    struct TailNode;

    using Node = std::variant<VertexNode, TailNode>;

    struct VertexNode {
        size_t v_id;
        std::vector<Node> children;
        bool is_end {false};
    };

    struct TailNode {
        size_t face_id;
    };

    std::vector<VertexNode> _root_nodes;
};
} // namespace prism::mesh