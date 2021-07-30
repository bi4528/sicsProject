#ifndef SICS_FORWARDCHECKING_DEGREESEQUENCEPRUNE_IND_H_
#define SICS_FORWARDCHECKING_DEGREESEQUENCEPRUNE_IND_H_

#include <iterator>
#include <utility>
#include <vector>

#include "adjacency_degreesortedlistmat.h"
#include "graph_traits.h"
#include "label_equivalence.h"
#include "multi_stack.h"

#include "stats.h"

namespace sics {

template <
    typename G,
    typename H,
    typename Callback,
    typename IndexOrderG,
    typename VertexEquiv = default_vertex_label_equiv<G, H>,
    typename EdgeEquiv = default_edge_label_equiv<G, H>>
void forwardchecking_degreesequenceprune_ind(
    G const & g,
    H const & h,
    Callback const & callback,
    IndexOrderG const & index_order_g,
    VertexEquiv const & vertex_equiv = VertexEquiv(),
    EdgeEquiv const & edge_equiv = EdgeEquiv()) {

  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;

  struct explorer {

    adjacency_degreesortedlistmat<
        IndexG,
        typename G::directed_category,
        typename G::vertex_label_type,
        typename G::edge_label_type> g;
    adjacency_degreesortedlistmat<
        IndexH,
        typename H::directed_category,
        typename H::vertex_label_type,
        typename H::edge_label_type> h;

    Callback callback;

    IndexOrderG const & index_order_g;

    vertex_equiv_helper<VertexEquiv> vertex_equiv;
    edge_equiv_helper<EdgeEquiv> edge_equiv;

    IndexG m;
    IndexH n;

    IndexG level;

    std::vector<IndexH> map;

    std::vector<char> M;
    bool M_get(IndexG u, IndexH v) {
      return M[u*n + v];
    }
    void M_set(IndexG u, IndexH v) {
      M[u*n + v] = true;
    }
    void M_unset(IndexG u, IndexH v) {
      M[u*n + v] = false;
    }
    void build_M() {
      for (IndexG u=0; u<m; ++u) {
        for (IndexH v=0; v<n; ++v) {
          if (vertex_equiv(g, u, h, v) &&
              degree_condition(g, u, h, v) &&
              degree_sequence_condition(g, u, h, v)) {
            M_set(u, v);
          }
        }
      }
    }
    multi_stack<std::pair<IndexG,IndexH>> M_mst;

    explorer(
        G const & g,
        H const & h,
        Callback const & callback,
        IndexOrderG const & index_order_g,
        VertexEquiv const & vertex_equiv,
        EdgeEquiv const & edge_equiv)
        : g{g},
          h{h},
          callback{callback},
          index_order_g{index_order_g},
          vertex_equiv{vertex_equiv},
          edge_equiv{edge_equiv},

          m{g.num_vertices()},
          n{h.num_vertices()},
          level{0},
          map(m, n),
          M(m * n, false),
          M_mst(m*n, m) {
      build_M();
    }

    bool explore() {
      SICS_STATS_STATE;
      if (level == m) {
        return callback(std::as_const(*this));
      } else {
        auto x = index_order_g[level];
        bool proceed = true;
        for (IndexH y=0; y<n; ++y) {
          if (M_get(x, y)) {
            M_mst.push_level();
            if (forward_check(y)) {
              map[x] = y;
              ++level;
              proceed = explore();
              --level;
              map[x] = n;
            }
            revert_M();
            M_mst.pop_level();
            if (!proceed) {
              break;
            }
          }
        }
        return proceed;
      }
    }

    bool forward_check(IndexH y) {
      auto x = index_order_g[level];

      bool not_empty = true;
      for (IndexG i=level+1; i<m && not_empty; ++i) {
        auto u = index_order_g[i];
        not_empty = false;

        if constexpr (is_directed_v<G>) {
          bool x_out = g.edge(x, u);
          bool x_in = g.edge(u, x);
          for (IndexH v=0; v<n; ++v) {
            if (M_get(u, v)) {
              if (v == y ||
                  (x_out != h.edge(y, v) || (x_out && !edge_equiv(g, x, u, h, y, v))) ||
                  (x_in != h.edge(v, y) || (x_in && !edge_equiv(g, u, x, h, v, y)))) {
                M_unset(u, v);
                M_mst.push({u, v});
              } else {
                not_empty = true;
              }
            }
          }
        } else {
          bool x_edge = g.edge(x, u);
          for (IndexH v=0; v<n; ++v) {
            if (M_get(u, v)) {
              if (v == y ||
                  (x_edge != h.edge(y, v) || (x_edge && !edge_equiv(g, x, u, h, y, v)))) {
                M_unset(u, v);
                M_mst.push({u, v});
              } else {
                not_empty = true;
              }
            }
          }
        }
      }
      return true;//not_empty;
    }

    void revert_M() {
      while (!M_mst.level_empty()) {
        auto [u, v] = M_mst.top();
        M_mst.pop();
        M_set(u, v);
      }
    }
  } e(g, h, callback, index_order_g, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_FORWARDCHECKING_DEGREESEQUENCEPRUNE_IND_H_
