#ifndef SICS_LAZYFORWARDCHECKING_LOW_IND_H_
#define SICS_LAZYFORWARDCHECKING_LOW_IND_H_

#include <iterator>
#include <utility>
#include <vector>
#include <stack>

#include "graph_traits.h"
#include "label_equivalence.h"

#include "stats.h"

namespace sics {

template <
    typename G,
    typename H,
    typename Callback,
    typename IndexOrderG,
    typename VertexEquiv = default_vertex_label_equiv<G, H>,
    typename EdgeEquiv = default_edge_label_equiv<G, H>>
void lazyforwardchecking_low_ind(
    G const & g,
    H const & h,
    Callback const & callback,
    IndexOrderG const & index_order_g,
    VertexEquiv const & vertex_equiv = VertexEquiv(),
    EdgeEquiv const & edge_equiv = EdgeEquiv()) {

  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;

  struct explorer {

    G const & g;
    H const & h;
    Callback callback;

    IndexOrderG const & index_order_g;

    vertex_equiv_helper<VertexEquiv> vertex_equiv;
    edge_equiv_helper<EdgeEquiv> edge_equiv;

    IndexG m;
    IndexH n;

    IndexG level;

    std::vector<IndexH> map;

    std::vector<IndexG> low;
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
          if (vertex_equiv(g, u, h, v)) {
            M_set(u, v);
          }
        }
      }
    }
    std::vector<std::stack<std::pair<IndexG,IndexH>>> M_sts;

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
          low(m, 0),
          M(m * n, false),
          M_sts(m) {
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
          if (M_get(x, y) &&
              consistency(y)) {
            for (IndexG i=level+1; i<m && level<low[i]; ++i) {
              low[i] = level;
            }
            map[x] = y;
            ++level;
            proceed = explore();
            --level;
            map[x] = n;
            revert_M();
            if (!proceed) {
              break;
            }
          }
        }
        low[level] = level;
        return proceed;
      }
    }

    bool consistency(IndexH y) {
      auto x = index_order_g[level];
      for (IndexG i=low[level]; i<level; ++i) {
        auto u = index_order_g[i];
        auto v = map[u];
        if (v == y) {
          M_unset(x, y);
          M_sts[i].emplace(x, y);
          return false;
        }
        auto x_out = g.edge(x, u);
        if (x_out != h.edge(y, v) || (x_out && !edge_equiv(g, x, u, h, y, v))) {
          M_unset(x, y);
          M_sts[i].emplace(x, y);
          return false;
        }
        if constexpr (is_directed_v<G>) {
          auto x_in = g.edge(u, x);
          if (x_in != h.edge(v, y) || (x_in && !edge_equiv(g, u, x, h, v, y))) {
            M_unset(x, y);
            M_sts[i].emplace(x, y);
            return false;
          }
        }
      }
      return true;
    }

    void revert_M() {
      while (!M_sts[level].empty()) {
        IndexG u;
        IndexH v;
        std::tie(u, v) = M_sts[level].top();
        M_sts[level].pop();
        M_set(u, v);
      }
    }
  } e(g, h, callback, index_order_g, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_LAZYFORWARDCHECKING_LOW_IND_H_
