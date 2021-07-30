#ifndef SICS_READ_GALV_H_
#define SICS_READ_GALV_H_

#include <istream>

namespace sics {

template <typename G>
G read_galv(std::istream & in) {
  using index_type = typename G::index_type;
  using vertex_label_type = typename G::vertex_label_type;

  index_type n;
  in >> n;

  G g(n);

  for (index_type u=0; u<n; ++u) {
    vertex_label_type vertex_label;
    in >> vertex_label;
    g.set_vertex_label(u, vertex_label);
    index_type cnt;
    in >> cnt;
    for (index_type j=0; j<cnt; ++j) {
      index_type v;
      in >> v;
      g.add_edge(u, v);
    }
  }
  return g;
}

}  // namespace sics

#endif  // SICS_READ_GALV_H_
