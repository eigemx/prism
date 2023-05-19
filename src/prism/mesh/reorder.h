#pragma once

#include <map>

#include "../types.h"
#include "pmesh.h"

namespace prism::mesh {

class CuthillMckee {
  public:
    CuthillMckee(PMesh& mesh);
    auto permute() -> std::vector<std::size_t>;
    auto reorder() -> void;

  private:
    class Node {
      public:
        Node(std::size_t id, std::size_t degree) : _id(id), _degree(degree) {}
        auto inline id() const -> std::size_t { return _id; }
        auto inline id() -> std::size_t& { return _id; }
        auto inline degree() const -> std::size_t { return _degree; }
        auto inline degree() -> std::size_t& { return _degree; }
        auto inline neighbors() const -> const std::vector<std::size_t>& { return _neighbors; }
        auto inline neighbors() -> std::vector<std::size_t>& { return _neighbors; }

      private:
        std::size_t _id {};
        std::size_t _degree {};
        std::vector<std::size_t> _neighbors {};
    };


    PMesh& _mesh;
    std::vector<Node> _nodes;
};
} // namespace prism::mesh