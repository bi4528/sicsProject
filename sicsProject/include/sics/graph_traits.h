#ifndef SICS_GRAPH_TRAITS_H_
#define SICS_GRAPH_TRAITS_H_

#include <type_traits>

namespace sics {

struct directed_category_tag {};
struct directed_tag : public directed_category_tag {};
struct undirected_tag : public directed_category_tag {};
struct bidirectional_tag : public directed_tag {};

template <
    typename DirectedCategoryTag,
    typename G>
struct is_directed_category_of {
  using G_directed_category = typename G::directed_category;
  static_assert(std::is_base_of<directed_category_tag, G_directed_category>::value);
  static constexpr bool value = std::is_base_of<DirectedCategoryTag, G_directed_category>::value;
};
template <
    typename DirectedCategoryTag,
    typename G>
inline constexpr bool is_directed_category_of_v = is_directed_category_of<DirectedCategoryTag, G>::value;


template <typename G> struct is_directed : public is_directed_category_of<directed_tag, G> {};
template <typename G> inline constexpr bool is_directed_v = is_directed<G>::value;
template <typename G> struct is_undirected : public is_directed_category_of<undirected_tag, G> {};
template <typename G> inline constexpr bool is_undirected_v = is_undirected<G>::value;
template <typename G> struct is_bidirectional : public is_directed_category_of<bidirectional_tag, G> {};
template <typename G> inline constexpr bool is_bidirectional_v = is_bidirectional<G>::value;

template <typename G>
struct is_vertex_labelled {
  using vertex_label_type = typename G::vertex_label_type;
  static constexpr bool value = ! std::is_same<void, vertex_label_type>::value;
};
template <typename G> inline constexpr bool is_vertex_labelled_v = is_vertex_labelled<G>::value;

template <typename G>
struct is_edge_labelled {
  using edge_label_type = typename G::edge_label_type;
  static constexpr bool value = ! std::is_same<void, edge_label_type>::value;
};
template <typename G> inline constexpr bool is_edge_labelled_v = is_edge_labelled<G>::value;

}  // namespace sics

#endif  // SICS_GRAPH_TRAITS_H_
