#pragma once


#include "pmesh.h"
#include "prism/types.h"

namespace prism::mesh {

class CuthillMckee {
  public:
    CuthillMckee(const SharedPtr<PMesh>& mesh) noexcept;
    auto reorder(bool reverse = false) -> void;

  private:
    auto permute() -> std::vector<size_t>;

    class Node {
      public:
        Node(size_t id, size_t degree) noexcept : _id(id), _degree(degree) {}

        auto inline id() const -> size_t { return _id; }
        auto inline id() -> size_t& { return _id; }

        auto inline degree() const -> size_t { return _degree; }
        auto inline degree() -> size_t& { return _degree; }

        auto inline neighbors() const -> const std::vector<size_t>& { return _neighbors; }
        auto inline neighbors() -> std::vector<size_t>& { return _neighbors; }

      private:
        size_t _id {};
        size_t _degree {};
        std::vector<size_t> _neighbors;
    };

    SharedPtr<PMesh> _mesh;
    std::vector<Node> _nodes;
};
} // namespace prism::mesh
