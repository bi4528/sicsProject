#ifndef SICS_READ_GF_H_
#define SICS_READ_GF_H_

#include <cstddef>

#include <istream>
#include <string>

namespace sics {

template <typename G>
G read_gf(std::istream & in) {
  using index_type = typename G::index_type;
  using vertex_label_type = typename G::vertex_label_type;

  std::string name;
  in >> name;

  index_type n;
  in >> n;

  G g(n);

  for (index_type u=0; u<n; ++u) {
    vertex_label_type vertex_label;
    in >> vertex_label;
    g.set_vertex_label(u, vertex_label);
  }

  std::size_t edge_count;
  in >> edge_count;
  for (std::size_t i=0; i<edge_count; ++i) {
    index_type u;
    index_type v;
    in >> u >> v;
    g.add_edge(u, v);
  }

  return g;
}

}  // namespace sics

#endif  // SICS_READ_GF_H_
