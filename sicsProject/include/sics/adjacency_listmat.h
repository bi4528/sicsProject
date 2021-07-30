#ifndef SICS_ADJACENCY_LISTMAT_H_
#define SICS_ADJACENCY_LISTMAT_H_

#include <vector>

#include <boost/range/iterator_range.hpp>

#include "graph_traits.h"
#include "adjacency_list.h"

namespace sics {

template <
    typename Index,
    typename DirectedCategory,
    typename VertexLabel = void,
    typename EdgeLabel = void>
class adjacency_listmat;

// TODO implement EdgeLabel != void

template <
    typename Index,
    typename DirectedCategory,
    typename VertexLabel>
class adjacency_listmat<Index, DirectedCategory, VertexLabel, void>
    : public adjacency_list<Index, DirectedCategory, VertexLabel, void> {
 private:
  using base = adjacency_list<Index, DirectedCategory, VertexLabel, void>;

 protected:
  typename base::index_type n;
  std::vector<bool> m_mat;

 public:
  template <
      typename Temp = typename base::directed_category,  // needed to make enable_if work
      typename std::enable_if_t<std::is_base_of<undirected_tag, Temp>::value, int> = 0>
  explicit adjacency_listmat(typename base::index_type n)
      : base(n),
        n{n},
        m_mat((n * (n+1)) / 2, false) {
  }

  template <
      typename G,
      typename std::enable_if_t<!std::is_integral<G>::value, int> = 0,
      typename Temp = typename base::directed_category,  // needed to make enable_if work
      typename std::enable_if_t<std::is_base_of<undirected_tag, Temp>::value, int> = 0>
  explicit adjacency_listmat(G const & g)
      : base(g),
        n{base::num_vertices()},
        m_mat((n * (n+1)) / 2, false) {
    for (typename base::index_type u=0; u<n; ++u) {
      for (auto e : g.edges(u)) {
        auto v = e.target;
        if (u >= v) {
          m_mat[(u * (u+1)) / 2 + v] = true;
        }
        else {
          m_mat[(v * (v+1)) / 2 + u] = true;
        }
      }
    }
  }

  template <
      typename Temp = typename base::directed_category,  // needed to make enable_if work
      typename std::enable_if_t<std::is_base_of<directed_tag, Temp>::value, int> = 0>
  explicit adjacency_listmat(typename base::index_type n)
      : base(n),
        n{n},
        m_mat(n * n, false) {
  }

  template <
      typename G,
      typename std::enable_if_t<!std::is_integral<G>::value, int> = 0,
      typename Temp = typename base::directed_category,  // needed to make enable_if work
      typename std::enable_if_t<std::is_base_of<directed_tag, Temp>::value, int> = 0>
  explicit adjacency_listmat(G const & g)
      : base(g),
        n{base::num_vertices()},
        m_mat(n * n, false) {
    for (typename base::index_type u=0; u<n; ++u) {
      for (auto oe : g.out_edges(u)) {
        m_mat[u*n + oe.target] = true;
      }
    }
  }

  void add_edge(typename base::index_type u, typename base::index_type v) {
    base::add_edge(u, v);
    if constexpr (std::is_base_of<undirected_tag, typename base::directed_category>::value) {
      if (u >= v) {
        m_mat[(u * (u+1)) / 2 + v] = true;
      } else {
        m_mat[(v * (v+1)) / 2 + u] = true;
      }
    } else { // directed
      m_mat[u*n + v] = true;
    }
  }

  bool edge(typename base::index_type u, typename base::index_type v) const {
    if constexpr (std::is_base_of<undirected_tag, typename base::directed_category>::value) {
      if (u >= v) {
        return m_mat[(u * (u+1)) / 2 + v];
      } else {
        return m_mat[(v * (v+1)) / 2 + u];
      }
    } else { // directed
      return m_mat[u*n + v];
    }
  }
};

}  // namespace sics

#endif  // SICS_ADJACENCY_LISTMAT_H_
