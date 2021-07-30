#ifndef SICS_VERTEX_ORDER_H_
#define SICS_VERTEX_ORDER_H_

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>
#include <set>

#include "graph_traits.h"
#include "graph_utilities.h"

namespace sics {

template <typename G>
std::vector<typename G::index_type> vertex_order_DEG(G const & g) {
  std::vector<typename G::index_type> order(g.num_vertices());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&g](auto u, auto v) {
    return g.degree(u) > g.degree(v);
  });
  return order;
}

template <typename G>
std::vector<typename G::index_type> vertex_order_RDEG(G const & g) {
  using Index = typename G::index_type;

  auto n = g.num_vertices();

  std::vector<Index> vertex_order(n);

  std::vector<bool> avail(n, true);
  for (Index idx=0; idx<n; ++idx) {
    auto bestn = n;
    int bestv = -1;
    for (Index i=0; i<n; ++i) {
      if (avail[i]) {
        auto i_out_edges = edges_or_out_edges(g, i);
        int rdeg = std::count_if(std::begin(i_out_edges), std::end(i_out_edges), [&avail](auto oe) {
          return !avail[oe.target];
        });
        if constexpr (is_directed_v<G>) {
          auto i_in_edges = g.in_edges(i);
          rdeg += std::count_if(std::begin(i_in_edges), std::end(i_in_edges), [&avail](auto ie) {
            return !avail[ie.target];
          });
        }

        if (bestn == n || std::forward_as_tuple(rdeg, g.degree(i), i) > std::forward_as_tuple(bestv, g.degree(bestn), bestn)) {
          bestn = i;
          bestv = rdeg;
        }
      }
    }
    avail[bestn] = false;
    vertex_order[idx] = bestn;
  }
  return vertex_order;
}

template <typename G>
std::vector<typename G::index_type> vertex_order_GreatestConstraintFirst(G const & g) {
  using Index = typename G::index_type;

  auto n = g.num_vertices();

  std::vector<Index> vertex_order(n);
  std::iota(std::begin(vertex_order), std::end(vertex_order), 0);

  enum struct Flag {
    vis,
    neigh,
    unv
  };

  std::vector<Flag> flags(n, Flag::unv);
  std::vector<std::tuple<Index,Index,Index,Index>> ranks(n);

  for (Index i=0; i<n; ++i) {
    std::get<2>(ranks[i]) = g.degree(i);
    std::get<3>(ranks[i]) = i;
  }

  for (Index m=0; m<n; ++m) {
    auto u_it = std::max_element(
        std::next(std::begin(vertex_order), m),
        std::end(vertex_order),
        [&ranks](auto u, auto v) {
          return ranks[u] < ranks[v];
        });
    Index u = *u_it;

    if (flags[u] == Flag::unv) {
      for (auto oe : edges_or_out_edges(g, u)) {
        --std::get<2>(ranks[oe.target]);
      }
      if constexpr (is_directed_v<G>) {
        for (auto ie : g.in_edges(u)) {
          --std::get<2>(ranks[ie.target]);
        }
      }
    } else if (flags[u] == Flag::neigh) {
      for (auto oe : edges_or_out_edges(g, u)) {
        --std::get<1>(ranks[oe.target]);
      }
      if constexpr (is_directed_v<G>) {
        for (auto ie : g.in_edges(u)) {
          --std::get<1>(ranks[ie.target]);
        }
      }
    }

    std::swap(vertex_order[m], *u_it);
    flags[u] = Flag::vis;

    for (auto oe : edges_or_out_edges(g, u)) {
      auto v = oe.target;
      ++std::get<0>(ranks[v]);
      if (flags[v] == Flag::unv) {
        flags[v] = Flag::neigh;
        for (auto v_oe : edges_or_out_edges(g, v)) {
          ++std::get<1>(ranks[v_oe.target]);
        }
        if constexpr (is_directed_v<G>) {
          for (auto v_ie : g.in_edges(v)) {
            ++std::get<1>(ranks[v_ie.target]);
          }
        }
      }
    }
    if constexpr (is_directed_v<G>) {
      for (auto ie : g.in_edges(u)) {
        auto v = ie.target;
        ++std::get<0>(ranks[v]);
        if (flags[v] == Flag::unv) {
          flags[v] = Flag::neigh;
          for (auto v_oe : g.out_edges(v)) {
            ++std::get<1>(ranks[v_oe.target]);
          }
          for (auto v_ie : g.in_edges(v)) {
            ++std::get<1>(ranks[v_ie.target]);
          }
        }
      }
    }
  }
  return vertex_order;
}

}  // namespace sics

#endif  // SICS_VERTEX_ORDER_H_
