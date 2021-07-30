#ifndef SICS_FORWARDCHECKING_BITSET_MRV_DEGREEPRUNE_IND_H_
#define SICS_FORWARDCHECKING_BITSET_MRV_DEGREEPRUNE_IND_H_

#include <iterator>
#include <utility>
#include <tuple>
#include <numeric>
#include <vector>

#include <boost/dynamic_bitset.hpp>

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
void forwardchecking_bitset_mrv_degreeprune_ind(
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

    std::vector<IndexG> index_order_g;

    std::vector<IndexH> map;

    std::vector<boost::dynamic_bitset<>> M;
    void build_M() {
      for (IndexG u=0; u<m; ++u) {
        for (IndexH v=0; v<n; ++v) {
          if (vertex_equiv(g, u, h, v) &&
              degree_condition(g, u, h, v)) {
            M[u].set(v);
          }
        }
      }
    }
    multi_stack<std::tuple<IndexG, boost::dynamic_bitset<>>> M_mst;

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
          h_bits(n),
          h_c_bits(n),
          level{0},
          index_order_g(m),
          map(m, n),
          M(m, boost::dynamic_bitset<>(n)),
          M_mst(m*n, m) {
      build_h_bits();
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
              return std::forward_as_tuple(M[a].count(), g.degree(a), a) < std::forward_as_tuple(M[b].count(), g.degree(b), b);
            });
        std::swap(index_order_g[level], *it);
        auto x = index_order_g[level];
        bool proceed = true;
        for (auto y=M[x].find_first(); y!=boost::dynamic_bitset<>::npos; y=M[x].find_next(y)) {
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
        return proceed;
      }
    }

    bool forward_check(IndexH y) {
      auto x = index_order_g[level];

      bool not_empty = true;
      for (IndexG i=level+1; i<m && not_empty; ++i) {
        auto u = index_order_g[i];

        M_mst.push({u, M[u]});

        M[u].reset(y);
        if (g.edge(x, u)) {
          M[u] &= std::get<0>(h_bits[y]);
        } else {
          M[u] &= std::get<0>(h_c_bits[y]);
        }

        if constexpr (is_directed_v<G>) {
          if (g.edge(u, x)) {
            M[u] &= std::get<1>(h_bits[y]);
          } else {
            M[u] &= std::get<1>(h_c_bits[y]);
          }
        }

        not_empty = M[u].any();
      }
      return not_empty;
    }

    void revert_M() {
      while (!M_mst.level_empty()) {
        auto & [u, row] = M_mst.top();
        M[u] = row;
        M_mst.pop();
      }
    }
  } e(g, h, callback, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_FORWARDCHECKING_BITSET_MRV_DEGREEPRUNE_IND_H_
