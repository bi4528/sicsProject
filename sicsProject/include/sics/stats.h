#ifndef SICS_STATS_H_
#define SICS_STATS_H_

#ifdef SICS_DO_STATS

#include <cstddef>
#include <ostream>

namespace sics {

struct sics_stats_t {
  std::size_t state;
  std::size_t vequiv;
  std::size_t unique;
  std::size_t top;
  std::size_t eequiv;
  std::size_t misc;

  friend std::ostream & operator<<(std::ostream & out, sics_stats_t const & stats) {
    return out <<
        stats.state/* << "," <<
        stats.vequiv << "," <<
        stats.unique << "," <<
        stats.top << "," <<
        stats.eequiv << "," <<
        stats.misc*/;
  }
} sics_stats{};

}  // namespace sics

#define SICS_STATS_STATE ++sics::sics_stats.state
#define SICS_STATS_STATE_E(x) (SICS_STATS_STATE, x)

#define SICS_STATS_VEQUIV ++sics::sics_stats.vequiv
#define SICS_STATS_VEQUIV_E(x) (SICS_STATS_VEQUIV, x)

#define SICS_STATS_UNIQUE ++sics::sics_stats.unique
#define SICS_STATS_UNIQUE_E(x) (SICS_STATS_UNIQUE, x)

#define SICS_STATS_TOP ++sics::sics_stats.top
#define SICS_STATS_TOP_E(x) (SICS_STATS_TOP, x)

#define SICS_STATS_EEQUIV ++sics::sics_stats.eequiv
#define SICS_STATS_EEQUIV_E(x) (SICS_STATS_EEQUIV, x)

#define SICS_STATS_MISC ++sics::sics_stats.misc
#define SICS_STATS_MISC_E(x) (SICS_STATS_MISC, x)

#define SICS_STATS_PRINT(out) ((out) << sics::sics_stats << std::endl)

#else  // SICS_DO_STATS

#define SICS_STATS_STATE
#define SICS_STATS_STATE_E(x) (x)

#define SICS_STATS_VEQUIV
#define SICS_STATS_VEQUIV_E(x) (x)

#define SICS_STATS_UNIQUE
#define SICS_STATS_UNIQUE_E(x) (x)

#define SICS_STATS_TOP
#define SICS_STATS_TOP_E(x) (x)

#define SICS_STATS_EEQUIV
#define SICS_STATS_EEQUIV_E(x) (x)

#define SICS_STATS_MISC
#define SICS_STATS_MISC_E(x) (x)

#define SICS_STATS_PRINT(out)

#endif  // SICS_DO_STATS

#endif  // SICS_STATS_H_
