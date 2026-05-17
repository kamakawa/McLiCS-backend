#ifndef RNG_POOL_H_
#define RNG_POOL_H_

#include <vector>
#include <gsl/gsl_rng.h>
#include <omp.h>

// RAII wrapper for a per-thread array of GSL RNG instances.
// Allocates one RNG per OpenMP thread and frees them automatically
// when the object goes out of scope (even on exceptions or early exit).
class RngPool {
 public:
  explicit RngPool(int n_threads = -1) {
    if (n_threads < 0) n_threads = omp_get_max_threads();
    gsl_rng_env_setup();
    rngs_.resize(n_threads);
    for (int i = 0; i < n_threads; i++) {
      rngs_[i] = gsl_rng_alloc(gsl_rng_ranlxs0);
      gsl_rng_set(rngs_[i], static_cast<unsigned long>(i));
    }
  }

  ~RngPool() {
    for (auto r : rngs_) gsl_rng_free(r);
  }

  // Non-copyable, movable
  RngPool(const RngPool&)            = delete;
  RngPool& operator=(const RngPool&) = delete;
  RngPool(RngPool&&)                 = default;
  RngPool& operator=(RngPool&&)      = default;

  gsl_rng** data()  { return rngs_.data(); }
  int       size()  const { return static_cast<int>(rngs_.size()); }
  gsl_rng*  operator[](int i) { return rngs_[i]; }

 private:
  std::vector<gsl_rng*> rngs_;
};

#endif
