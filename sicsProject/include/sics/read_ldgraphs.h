#ifndef SICS_READ_LDGRAPHS_H_
#define SICS_READ_LDGRAPHS_H_

#include <cstdint>

#include <istream>
#include <string>

namespace sics {

template <typename G>
G read_ldgraphs_unl(std::istream & in) {
  using index_type = typename G::index_type;

  index_type n;
  in >> n;

  G g(n);

  for (index_type i=0; i<n; ++i) {
    index_type u;
    std::string ignore;
    in >> u >> ignore;
  }

  for (index_type i=0; i<n; ++i) {
    int cnt;
    in >> cnt;
    for (decltype(cnt) j=0; j<cnt; ++j) {
      index_type u, v;
      in >> u >> v;
      g.add_edge(u, v);
    }
  }
  return g;
}

template <
    typename G,
    typename VertexLabel = typename G::vertex_label_type>
G read_ldgraphs_lab(std::istream & in) {
  using index_type = typename G::index_type;

  index_type n;
  in >> n;

  G g(n);

  for (index_type i=0; i<n; ++i) {
    index_type u;
    VertexLabel label;
    in >> u >> label;
    g.set_vertex_label(u, label);
  }

  for (index_type i=0; i<n; ++i) {
    int cnt;
    in >> cnt;
    for (decltype(cnt) j=0; j<cnt; ++j) {
      index_type u, v;
      in >> u >> v;
      g.add_edge(u, v);
    }
  }
  return g;
}

}  // namespace sics

#endif  // SICS_READ_LDGRAPHS_H_
