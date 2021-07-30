#ifndef SICS_LAZYFORWARDCHECKING_LOW_BITSET_DEGREEPRUNE_IND_H_
#define SICS_LAZYFORWARDCHECKING_LOW_BITSET_DEGREEPRUNE_IND_H_

#include <iterator>
#include <utility>
#include <vector>
#include <stack>

#include <boost/dynamic_bitset.hpp>

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
void lazyforwardchecking_low_bitset_degreeprune_ind(
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

    using bits_type = std::conditional_t<
        is_directed_v<H>,
        std::tuple<boost::dynamic_bitset<>, boost::dynamic_bitset<>>,
        std::tuple<boost::dynamic_bitset<>>>;

    std::vector<bits_type> h_bits;
    std::vector<bits_type> h_c_bits;
    void build_h_bits() {
      for (IndexH i=0; i<n; ++i) {
        std::get<0>(h_bits[i]).resize(n);
        std::get<0>(h_c_bits[i]).resize(n);
        if constexpr (is_directed_v<H>) {
          std::get<1>(h_bits[i]).resize(n);
          std::get<1>(h_c_bits[i]).resize(n);
        }
      }

      for (IndexH i=0; i<n; ++i) {
        for (IndexH j=0; j<n; ++j) {
          if (h.edge(i, j)) {
            std::get<0>(h_bits[i]).set(j);
            if constexpr (is_directed_v<H>) {
              std::get<1>(h_bits[j]).set(i);
            }
          } else {
            std::get<0>(h_c_bits[i]).set(j);
            if constexpr (is_directed_v<H>) {
              std::get<1>(h_c_bits[j]).set(i);
            }
          }
        }
      }
    }

    IndexG level;

    std::vector<IndexH> map;

    std::vector<IndexG> low;
    std::vector<boost::dynamic_bitset<>> Ms;
    boost::dynamic_bitset<> & M_get(IndexG level, IndexG i) {
      return Ms[(level * (2*m - level + 1)) / 2 + i - level];
    }
    void build_M() {
      for (IndexG i=0; i<m; ++i) {
        auto u = index_order_g[i];
        for (IndexH v=0; v<n; ++v) {
          if (vertex_equiv(g, u, h, v) &&
              degree_condition(g, u, h, v)) {
            M_get(0, i).set(v);
          }
        }
      }
    }

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
          h_bits(n),
          h_c_bits(n),
          level{0},
          map(m, n),
          low(m, 0),
          Ms((m*(m+1))/2, boost::dynamic_bitset<>(n)) {
      build_h_bits();
      build_M();
    }

    bool explore() {
      SICS_STATS_STATE;
      if (level == m) {
        return callback(std::as_const(*this));
      } else {
        auto x = index_order_g[level];
        lazy_forward_check();
        bool proceed = true;
        for (auto y=M_get(level, level).find_first(); y!=boost::dynamic_bitset<>::npos; y=M_get(level, level).find_next(y)) {
          for (IndexG i=level+1; i<m && level<low[i]; ++i) {
            low[i] = level;
          }
          map[x] = y;
          ++level;
          proceed = explore();
          --level;
          map[x] = n;
          if (!proceed) {
            break;
          }
        }
        low[level] = level;
        return proceed;
      }
    }

    void lazy_forward_check() {
      auto x = index_order_g[level];
      for (IndexG i=low[level]; i<level; ++i) {
        auto u = index_order_g[i];
        auto v = map[u];
        M_get(i+1, level) = M_get(i, level);
        M_get(i+1, level).reset(v);
        if (g.edge(u, x)) {
          M_get(i+1, level) &= std::get<0>(h_bits[v]);
        } else {
          M_get(i+1, level) &= std::get<0>(h_c_bits[v]);
        }
        if constexpr (is_directed_v<G>) {
          if (g.edge(x, u)) {
            M_get(i+1, level) &= std::get<1>(h_bits[v]);
          } else {
            M_get(i+1, level) &= std::get<1>(h_c_bits[v]);
          }
        }
      }
    }
  } e(g, h, callback, index_order_g, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_LAZYFORWARDCHECKING_LOW_BITSET_DEGREEPRUNE_IND_H_
