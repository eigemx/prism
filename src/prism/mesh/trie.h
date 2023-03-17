#pragma once

#include <memory>
#include <optional>
#include <variant>
#include <vector>

namespace prism::mesh {
class FacesLookupTrie {
  public:
    FacesLookupTrie(std::size_t n_vertices);

    void insert(const std::vector<std::size_t>& face_labels, std::size_t face_id);
    auto find(const std::vector<std::size_t>& face_labels) const -> std::optional<std::size_t>;

  private:
    // forward declaration
    struct VertexNode;
    struct TailNode;
    using Node = std::variant<VertexNode, TailNode>;

    struct VertexNode {
        std::size_t v_id;
        std::vector<Node> children;
        bool is_end {false};
    };

    struct TailNode {
        std::size_t face_id;
    };

    std::vector<VertexNode> _root_nodes;
};
} // namespace prism::mesh