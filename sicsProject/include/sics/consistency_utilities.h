#ifndef SICS_CONSISTENCY_UTILITIES_H_
#define SICS_CONSISTENCY_UTILITIES_H_

#include <optional>
#include <tuple>

#include "graph_traits.h"
#include "graph_utilities.h"

namespace sics {

template <
    typename G,
    typename H,
    typename Map,
    typename Inv,
    typename EdgeEquiv>
bool adjacent_consistency_mono(
    G const & g,
    typename G::index_type u,
    H const & h,
    typename H::index_type v,
    Map const & map,
    Inv const & inv,
    EdgeEquiv const & edge_equiv) {
  auto n = h.num_vertices();
  for (auto oe : edges_or_out_edges(g, u)) {
    auto i = oe.target;
    auto j = map[i];
    if (j != n) {
      if (!h.edge(v, j) || !edge_equiv(g, u, oe, h, v, j)) {
        return false;
      }
    }
  }
  if constexpr (is_directed_v<G>) {
    for (auto ie : g.in_edges(u)) {
      auto i = ie.target;
      auto j = map[i];
      if (j != n) {
        if (!h.edge(j, v) || !edge_equiv(g, ie, u, h, j, v)) {
          return false;
        }
      }
    }
  }
  return true;
}

template <
    typename G,
    typename H,
    typename Map,
    typename Inv,
    typename EdgeEquiv>
bool adjacent_consistency_ind(
    G const & g,
    typename G::index_type u,
    H const & h,
    typename H::index_type v,
    Map const & map,
    Inv const & inv,
    EdgeEquiv const & edge_equiv) {
  if (!adjacent_consistency_mono(g, u, h, v, map, inv, edge_equiv)) {
    return false;
  }

  auto m = g.num_vertices();
  for (auto oe : edges_or_out_edges(h, v)) {
    auto j = oe.target;
    auto i = inv[j];
    if (i != m) {
      if (!g.edge(u, i)) {
        return false;
      }
    }
  }
  if constexpr (is_directed_v<H>) {
    for (auto ie : h.in_edges(v)) {
      auto j = ie.target;
      auto i = inv[j];
      if (i != m) {
        if (!g.edge(i, u)) {
          return false;
        }
      }
    }
  }
  return true;
}

template <
    typename G,
    typename H,
    typename Map,
    typename Inv,
    typename EdgeEquiv>
auto h_adjacent_consistency_mono(
    G const & g,
    typename G::index_type u,
    H const & h,
    typename H::index_type v,
    Map const & map,
    Inv const & inv,
    EdgeEquiv const & edge_equiv) {
  using h_count_type = std::conditional_t<
      is_directed_v<H>,
      std::tuple<typename H::index_type, typename H::index_type>,
      typename H::index_type>;
  h_count_type v_count{};

  auto m = g.num_vertices();
  for (auto oe : edges_or_out_edges(h, v)) {
    auto j = oe.target;
    auto i = inv[j];
    if (i != m) {
      if (!g.edge(u, i) || !edge_equiv(g, u, i, h, v, oe)) {
        return std::optional<h_count_type>{};
      }
      if constexpr (is_directed_v<H>) {
        ++std::get<0>(v_count);
      } else {
        ++v_count;
      }
    }
  }
  if constexpr (is_directed_v<H>) {
    for (auto ie : h.in_edges(v)) {
      auto j = ie.target;
      auto i = inv[j];
      if (i != m) {
        if (!g.edge(i, u) || !edge_equiv(g, i, u, h, ie, v)) {
          return std::optional<h_count_type>{};
        }
        ++std::get<1>(v_count);
      }
    }
  }
  return std::optional{v_count};
}

template <typename G, typename H>
bool degree_condition(
    G const & g,
    typename G::index_type u,
    H const & h,
    typename H::index_type v) {
  if constexpr (is_directed_v<G>) {
    return g.out_degree(u) <= h.out_degree(v) && g.in_degree(u) <= h.in_degree(v);
  } else {
    return g.degree(u) <= h.degree(v);
  }
}

template <typename G, typename H>
bool degree_sequence_condition(
    G const & g,
    typename G::index_type u,
    H const & h,
    typename H::index_type v) {
  if constexpr (is_directed_v<G>) {
    auto u_out = g.out_edges(u);
    auto v_out = h.out_edges(v);
    for (typename G::index_type i=0; i<g.out_degree(u); ++i) {
      auto ui = u_out[i].target;
      auto vi = v_out[i].target;
      if (g.degree(ui) > h.degree(vi)) {
        //return false;
      }
    }
    auto u_in = g.in_edges(u);
    auto v_in = h.in_edges(v);
    for (typename G::index_type i=0; i<g.in_degree(u); ++i) {
      auto ui = u_in[i].target;
      auto vi = v_in[i].target;
      if (g.degree(ui) > h.degree(vi)) {
        return false;
      }
    }
    return true;
  } else {
    auto u_edges = g.edges(u);
    auto v_edges = h.edges(v);
    for (typename G::index_type i=0; i<g.degree(u); ++i) {
      auto ui = u_edges[i].target;
      auto vi = v_edges[i].target;
      if (g.degree(ui) > h.degree(vi)) {
        return false;
      }
    }
    return true;
  }
}

}  // namespace sics

#endif  // SICS_CONSISTENCY_UTILITIES_H_
