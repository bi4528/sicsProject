#ifndef SICS_LABEL_EQUIVALENCE_H_
#define SICS_LABEL_EQUIVALENCE_H_

#include "graph_traits.h"

namespace sics {

template <
    typename G,
    typename H,
    typename LabelEquiv = void,
    typename SFINAE = void>
struct default_vertex_label_equiv;

template <
    typename G,
    typename H,
    typename LabelEquiv>
struct default_vertex_label_equiv<
    G,
    H,
    LabelEquiv,
    std::enable_if_t<is_vertex_labelled_v<G> && is_vertex_labelled_v<H>>> {

  LabelEquiv label_equiv;

  default_vertex_label_equiv(LabelEquiv const & label_equiv)
      : label_equiv{label_equiv} {
  }

  bool operator()(
      typename G::index_type u,
      typename G::vertex_label_type u_label,
      typename H::index_type v,
      typename H::vertex_label_type v_label) const {
    return label_equiv(u_label, v_label);
  }
};
template <
    typename G,
    typename H>
struct default_vertex_label_equiv<
    G,
    H,
    void,
    std::enable_if_t<is_vertex_labelled_v<G> && is_vertex_labelled_v<H>>> {
  bool operator()(
      typename G::index_type u,
      typename G::vertex_label_type u_label,
      typename H::index_type v,
      typename H::vertex_label_type v_label) const {
    return u_label == v_label;
  }
};
template <
    typename G,
    typename H>
struct default_vertex_label_equiv<
    G,
    H,
    void,
    std::enable_if_t<!is_vertex_labelled_v<G> && !is_vertex_labelled_v<H>>> {
  bool operator()(
      typename G::index_type u,
      typename H::index_type v) const {
    return true;
  }
};


template <
    typename G,
    typename H,
    typename LabelEquiv = void,
    typename SFINAE = void>
struct default_edge_label_equiv;

template <
    typename G,
    typename H,
    typename LabelEquiv>
struct default_edge_label_equiv<
    G,
    H,
    LabelEquiv,
    std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>>> {

  LabelEquiv label_equiv;

  default_edge_label_equiv(LabelEquiv const & label_equiv)
      : label_equiv{label_equiv} {
  }

  bool operator()(
      typename G::index_type u0,
      typename G::index_type u1,
      typename G::vertex_label_type u_label,
      typename H::index_type v0,
      typename H::index_type v1,
      typename H::vertex_label_type v_label) const {
    return label_equiv(u_label, v_label);
  }
};
template <
    typename G,
    typename H>
struct default_edge_label_equiv<
    G,
    H,
    void,
    std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>>> {
  bool operator()(
      typename G::index_type u0,
      typename G::index_type u1,
      typename G::vertex_label_type u_label,
      typename H::index_type v0,
      typename H::index_type v1,
      typename H::vertex_label_type v_label) const {
    return u_label == v_label;
  }
};
template <
    typename G,
    typename H>
struct default_edge_label_equiv<
    G,
    H,
    void,
    std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>>> {
  bool operator()(
      typename G::index_type u0,
      typename G::index_type u1,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return true;
  }
};


template <typename VertexEquiv>
struct vertex_equiv_helper {
  VertexEquiv vertex_equiv;

  vertex_equiv_helper(VertexEquiv const & vertex_equiv)
      : vertex_equiv{vertex_equiv} {
  }

  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_vertex_labelled_v<G> && is_vertex_labelled_v<H>, int> = 0>
  bool operator()(G const & g, typename G::index_type u, H const & h, typename H::index_type v) const {
    return vertex_equiv(u, g.get_vertex_label(u), v, h.get_vertex_label(v));
  }

  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_vertex_labelled_v<G> && !is_vertex_labelled_v<H>, int> = 0>
  bool operator()(G const & g, typename G::index_type u, H const & h, typename H::index_type v) const {
    return vertex_equiv(u, v);
  }
};

template <typename EdgeEquiv>
struct edge_equiv_helper {
  EdgeEquiv edge_equiv;

  edge_equiv_helper(EdgeEquiv const & edge_equiv)
      : edge_equiv{edge_equiv} {
  }

  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::index_type u1,
      H const & h,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return edge_equiv(
        u0, u1, g.get_edge_label(u0, u1),
        v0, v1, h.get_edge_label(v0, v1));
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::index_type u1,
      H const & h,
      typename H::index_type v0,
      typename H::half_edge_type hf) const {
    return edge_equiv(
        u0, u1, g.get_edge_label(u0, u1),
        v0, hf.target, hf.label);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::half_edge_type he,
      H const & h,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return edge_equiv(
        u0, he.target, he.label,
        v0, v1, h.get_edge_label(v0, v1));
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::half_edge_type he,
      H const & h,
      typename H::index_type v0,
      typename H::half_edge_type hf) const {
    return edge_equiv(
        u0, he.target, he.label,
        v0, hf.target, hf.label);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::index_type u1,
      H const & h,
      typename H::half_edge_type hf,
      typename H::index_type v1) const {
    return edge_equiv(
        u0, u1, g.get_edge_label(u0, u1),
        hf.target, v1, hf.label);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::half_edge_type he,
      typename G::index_type u1,
      H const & h,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return edge_equiv(
        he.target, u1, he.label,
        v0, v1, h.get_edge_label(v0, v1));
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<is_edge_labelled_v<G> && is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::half_edge_type he,
      typename G::index_type u1,
      H const & h,
      typename H::half_edge_type hf,
      typename H::index_type v1) const {
    return edge_equiv(
        he.target, u1, he.label,
        hf.target, v1, hf.label);
  }

  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::index_type u1,
      H const & h,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return edge_equiv(
        u0, u1,
        v0, v1);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::index_type u1,
      H const & h,
      typename H::index_type v0,
      typename H::half_edge_type hf) const {
    return edge_equiv(
        u0, u1,
        v0, hf.target);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::half_edge_type he,
      H const & h,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return edge_equiv(
        u0, he.target,
        v0, v1);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::half_edge_type he,
      H const & h,
      typename H::index_type v0,
      typename H::half_edge_type hf) const {
    return edge_equiv(
        u0, he.target,
        v0, hf.target);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::index_type u0,
      typename G::index_type u1,
      H const & h,
      typename H::half_edge_type hf,
      typename H::index_type v1) const {
    return edge_equiv(
        u0, u1,
        hf.target, v1);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::half_edge_type he,
      typename G::index_type u1,
      H const & h,
      typename H::index_type v0,
      typename H::index_type v1) const {
    return edge_equiv(
        he.target, u1,
        v0, v1);
  }
  template <
      typename G,
      typename H,
      typename std::enable_if_t<!is_edge_labelled_v<G> && !is_edge_labelled_v<H>, int> = 0>
  bool operator()(
      G const & g,
      typename G::half_edge_type he,
      typename G::index_type u1,
      H const & h,
      typename H::half_edge_type hf,
      typename H::index_type v1) const {
    return edge_equiv(
        he.target, u1,
        hf.target, v1);
  }
};

}  // namespace sics

#endif  // SICS_LABEL_EQUIVALENCE_H_
