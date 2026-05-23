//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// Adapted for superassp/SPTK integration:
// All symbols wrapped in namespace sptk::world:: to match SPTK convention.
//-----------------------------------------------------------------------------
#ifndef WORLD_STONEMASK_H_
#define WORLD_STONEMASK_H_

#include "world/macrodefinitions.h"

#if 1
namespace sptk {
namespace world {
#endif

WORLD_BEGIN_C_DECLS

void StoneMask(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    double *refined_f0);

WORLD_END_C_DECLS

#if 1
}  // namespace world
}  // namespace sptk
#endif

#endif  // WORLD_STONEMASK_H_
