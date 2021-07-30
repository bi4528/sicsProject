#ifndef SICS_BACKJUMPING_DEGREEPRUNE_IND_H_
#define SICS_BACKJUMPING_DEGREEPRUNE_IND_H_

#include <iterator>
#include <utility>
#include <vector>

#include "graph_traits.h"
#include "label_equivalence.h"
#include "consistency_utilities.h"

#include "stats.h"

namespace sics {

template <
    typename G,
    typename H,
    typename Callback,
    typename IndexOrderG,
    typename VertexEquiv = default_vertex_label_equiv<G, H>,
    typename EdgeEquiv = default_edge_label_equiv<G, H>>
void backjumping_degreeprune_ind(
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

    using x_it_type = typename IndexOrderG::const_iterator;
    x_it_type x_it;
    x_it_type backjump_it;

    std::vector<IndexH> map;

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
          x_it(std::cbegin(index_order_g)),
          backjump_it{x_it},
          map(m, n) {
    }

    bool explore() {
      SICS_STATS_STATE;
      if (x_it == std::cend(index_order_g)) {
        return callback(std::as_const(*this));
      } else {
        auto x = *x_it;
        bool proceed = true;
        backjump_it = std::next(x_it);
        x_it_type latest_it = std::cbegin(index_order_g);
        for (IndexH y=0; y<n; ++y) {
          x_it_type culprit_it = consistency(y);
          if (culprit_it > x_it) {
            map[x] = y;
            ++x_it;
            proceed = explore();
            --x_it;
            map[x] = n;
            if (!proceed || backjump_it <= x_it) {
              break;
            }
          }
          if (culprit_it > latest_it) {
            latest_it = culprit_it;
          }
        }
        if (backjump_it > x_it && latest_it <= x_it) {
          backjump_it = latest_it;
        }
        return proceed;
      }
    }

    x_it_type consistency(IndexH y) {
      auto x = *x_it;

      if (!vertex_equiv(g, x, h, y) ||
          !degree_condition(g, x, h, y)) {
        return std::cbegin(index_order_g);
      }

      x_it_type it;
      for (it=std::cbegin(index_order_g); it!=x_it; ++it) {
        auto u = *it;
        auto v = map[u];
        if (v == y) {
          break;
        }
        auto x_out = g.edge(x, u);
        if (x_out != h.edge(y, v) || (x_out && !edge_equiv(g, x, u, h, y, v))) {
          break;
        }
        if constexpr (is_directed_v<G>) {
          auto x_in = g.edge(u, x);
          if (x_in != h.edge(v, y) || (x_in && !edge_equiv(g, u, x, h, v, y))) {
            break;
          }
        }
      }
      return std::next(it);
    }
  } e(g, h, callback, index_order_g, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_BACKJUMPING_DEGREEPRUNE_IND_H_
