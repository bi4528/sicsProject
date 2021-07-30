#ifndef SICS_FORWARDCHECKING_MRV_DEGREEPRUNE_IND_H_
#define SICS_FORWARDCHECKING_MRV_DEGREEPRUNE_IND_H_

#include <iterator>
#include <utility>
#include <tuple>
#include <numeric>
#include <vector>

#include "graph_traits.h"
#include "label_equivalence.h"
#include "consistency_utilities.h"
#include "multi_stack.h"

#include "stats.h"

namespace sics {

template <
    typename G,
    typename H,
    typename Callback,
    typename VertexEquiv = default_vertex_label_equiv<G, H>,
    typename EdgeEquiv = default_edge_label_equiv<G, H>>
void forwardchecking_mrv_degreeprune_ind(
    G const & g,
    H const & h,
    Callback const & callback,
    VertexEquiv const & vertex_equiv = VertexEquiv(),
    EdgeEquiv const & edge_equiv = EdgeEquiv()) {

  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;

  struct explorer {

    G const & g;
    H const & h;
    Callback callback;

    vertex_equiv_helper<VertexEquiv> vertex_equiv;
    edge_equiv_helper<EdgeEquiv> edge_equiv;

    IndexG m;
    IndexH n;

    IndexG level;

    std::vector<IndexG> index_order_g;

    std::vector<IndexH> map;

    std::vector<IndexH> num_candidates;
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
              degree_condition(g, u, h, v)) {
            M_set(u, v);
            ++num_candidates[u];
          }
        }
      }
    }
    multi_stack<std::pair<IndexG,IndexH>> M_mst;

    explorer(
        G const & g,
        H const & h,
        Callback const & callback,
        VertexEquiv const & vertex_equiv,
        EdgeEquiv const & edge_equiv)
        : g{g},
          h{h},
          callback{callback},
          vertex_equiv{vertex_equiv},
          edge_equiv{edge_equiv},

          m{g.num_vertices()},
          n{h.num_vertices()},
          level{0},
          index_order_g(m),
          map(m, n),
          num_candidates(m, 0),
          M(m * n, false),
          M_mst(m*n, m) {
      std::iota(index_order_g.begin(), index_order_g.end(), 0);
      build_M();
    }

    bool explore() {
      SICS_STATS_STATE;
      if (level == m) {
        return callback(std::as_const(*this));
      } else {
        auto it = std::min_element(
            std::next(index_order_g.begin(), level),
            index_order_g.end(),
            [this](auto a, auto b) {
              return std::forward_as_tuple(num_candidates[a], g.degree(a), a) < std::forward_as_tuple(num_candidates[b], g.degree(b), b);
            });
        std::swap(index_order_g[level], *it);
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
                --num_candidates[u];
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
                --num_candidates[u];
                M_mst.push({u, v});
              } else {
                not_empty = true;
              }
            }
          }
        }
      }
      return not_empty;
    }

    void revert_M() {
      while (!M_mst.level_empty()) {
        auto [u, v] = M_mst.top();
        M_mst.pop();
        M_set(u, v);
        ++num_candidates[u];
      }
    }
  } e(g, h, callback, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_FORWARDCHECKING_MRV_DEGREEPRUNE_IND_H_
