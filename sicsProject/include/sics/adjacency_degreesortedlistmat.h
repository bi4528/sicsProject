#ifndef SICS_ADJACENCY_DEGREESORTEDLISTMAT_H_
#define SICS_ADJACENCY_DEGREESORTEDLISTMAT_H_

#include <type_traits>
#include <tuple>
#include <algorithm>

#include "adjacency_listmat.h"
#include "graph_traits.h"

namespace sics {

template <
    typename Index,
    typename DirectedCategory,
    typename VertexLabel = void,
    typename EdgeLabel = void>
class adjacency_degreesortedlistmat
    : public adjacency_listmat<Index, DirectedCategory, VertexLabel, EdgeLabel> {
 private:
  using base = adjacency_listmat<Index, DirectedCategory, VertexLabel, EdgeLabel>;

 public:

  template <
      typename G,
      typename std::enable_if_t<!std::is_integral<G>::value, int> = 0>
  explicit adjacency_degreesortedlistmat(G const & g)
      : base(g) {
    auto n = base::num_vertices();
    for (Index u=0; u<n; ++u) {
      if constexpr (std::is_base_of_v<directed_tag, DirectedCategory>) {
        std::sort(base::m_out_edges[u].begin(), base::m_out_edges[u].end(), [this](auto ha, auto hb) {
          auto a = ha.target;
          auto b = hb.target;
          return std::forward_as_tuple(base::degree(a), a) > std::forward_as_tuple(base::degree(b), b);
        });
        std::sort(base::m_in_edges[u].begin(), base::m_in_edges[u].end(), [this](auto ha, auto hb) {
          auto a = ha.target;
          auto b = hb.target;
          return std::forward_as_tuple(base::degree(a), a) > std::forward_as_tuple(base::degree(b), b);
        });
      } else {
        std::sort(base::m_edges[u].begin(), base::m_edges[u].end(), [this](auto ha, auto hb) {
          auto a = ha.target;
          auto b = hb.target;
          return std::forward_as_tuple(base::degree(a), a) > std::forward_as_tuple(base::degree(b), b);
        });
      }
    }
  }
};

}  // namespace sics


#endif  // SICS_ADJACENCY_DEGREESORTEDLISTMAT_H_
