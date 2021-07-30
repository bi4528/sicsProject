#ifndef SICS_FORWARDCHECKING_BITSET_DEGREEPRUNE_AC1_IND_H_
#define SICS_FORWARDCHECKING_BITSET_DEGREEPRUNE_AC1_IND_H_

#include <iterator>
#include <utility>
#include <numeric>
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
void forwardchecking_bitset_degreeprune_ac1_ind(
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

    std::vector<IndexG> index_order_g;

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
      ac1();
    }
    std::stack<std::vector<boost::dynamic_bitset<>>> M_st;

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
          M(m, boost::dynamic_bitset<>(n)) {
      build_h_bits();
      build_M();
    }

    bool explore() {
      SICS_STATS_STATE;
      if (level == m) {
        return callback(std::as_const(*this));
      } else {
        auto x = index_order_g[level];
        bool proceed = true;
        for (auto y=M[x].find_first(); y!=boost::dynamic_bitset<>::npos; y=M[x].find_next(y)) {
          M_st.push(M);
          if (forward_check(y) && ac1()) {
            map[x] = y;
            ++level;
            proceed = explore();
            --level;
            map[x] = n;
          }
          M = M_st.top();
          M_st.pop();
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

    bool ac1() {
      bool change;
      do {
        change = false;
        for (IndexG i0=level+1; i0<m; ++i0) {
          auto u0 = index_order_g[i0];
          for (auto v0=M[u0].find_first(); v0!=boost::dynamic_bitset<>::npos; v0=M[u0].find_next(v0)) {
            for (IndexG u1=0; u1<m; ++u1) {
              if (u1 != u0) {
                bool exists = false;
                for (auto v1=M[u1].find_first(); v1!=boost::dynamic_bitset<>::npos; v1=M[u1].find_next(v1)) {
                  if constexpr (is_directed_v<G>) {
                    if (v0 != v1 &&
                        (g.edge(u0, u1) == h.edge(v0, v1)) && (!g.edge(u0, u1) || edge_equiv(g, u0, u1, h, v0, v1)) &&
                        (g.edge(u1, u0) == h.edge(v1, v0)) && (!g.edge(u1, u0) || edge_equiv(g, u1, u0, g, v1, v0))) {
                      exists = true;
                      break;
                    }
                  } else {
                    if (v0 != v1 &&
                        g.edge(u0, u1) == h.edge(v0, v1) && (!g.edge(u0, u1) || edge_equiv(g, u0, u1, h, v0, v1))) {
                      exists = true;
                      break;
                    }
                  }
                }
                if (!exists) {
                  M[u0].reset(v0);
                  change = true;
                  break;
                }
              }
            }
          }
          if (!M[u0].any()) {
            return false;
          }
        }
      } while (change);
      return true;
    }
  } e(g, h, callback, index_order_g, vertex_equiv, edge_equiv);

  e.explore();
}

}  // namespace sics

#endif  // SICS_FORWARDCHECKING_BITSET_DEGREEPRUNE_AC1_IND_H_
