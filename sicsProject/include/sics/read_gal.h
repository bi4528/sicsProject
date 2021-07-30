#ifndef SICS_READ_GAL_H_
#define SICS_READ_GAL_H_

#include <istream>

namespace sics {

template <typename G>
G read_gal(std::istream & in) {
  using index_type = typename G::index_type;

  index_type n;
  in >> n;

  G g(n);

  for (index_type u=0; u<n; ++u) {
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

#endif  // SICS_READ_GAL_H_
