//-----------------------------------------------------------------------------
// Per-call pseudorandom state for CheapTrick() and D4C().
//
// The 2017 WORLD sources these files derive from drew their safeguard noise
// from randn() in matlabfunctions.cc, whose xorshift state is a process-wide
// global that is never reseeded. Results therefore depended on everything
// called earlier in the R session, and the OpenMP frame loop in D4C() raced
// on that state. Current upstream WORLD (mmorise/World) seeds a RandnState at
// the start of every CheapTrick()/D4C() call; this header does the same with
// the identical generator, so a call's output depends only on its inputs.
//-----------------------------------------------------------------------------
#ifndef SUPERASSP_WORLD_RANDN_STATE_H_
#define SUPERASSP_WORLD_RANDN_STATE_H_

#include <stdint.h>

namespace sptk {
namespace world {

struct RandnState {
  uint32_t x, y, z, w;
};

inline void randn_reseed(RandnState *s) {
  s->x = 123456789;
  s->y = 362436069;
  s->z = 521288629;
  s->w = 88675123;
}

inline uint32_t randn_step(RandnState *s) {
  uint32_t t = s->x ^ (s->x << 11);
  s->x = s->y;
  s->y = s->z;
  s->z = s->w;
  s->w = (s->w ^ (s->w >> 19)) ^ (t ^ (t >> 8));
  return s->w;
}

// Same arithmetic as randn() in matlabfunctions.cc (sum of 12 uniforms).
inline double randn(RandnState *s) {
  uint32_t tmp = randn_step(s) >> 4;
  for (int i = 0; i < 11; ++i) tmp += randn_step(s) >> 4;
  return tmp / 268435456.0 - 6.0;
}

// Advance the state past n randn() draws without computing them.
inline void randn_skip(RandnState *s, long n) {
  for (long i = 0; i < 12 * n; ++i) randn_step(s);
}

}  // namespace world
}  // namespace sptk

#endif  // SUPERASSP_WORLD_RANDN_STATE_H_
