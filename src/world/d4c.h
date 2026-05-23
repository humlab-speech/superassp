//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// Adapted for superassp/SPTK integration:
// All symbols wrapped in namespace sptk::world:: to match SPTK convention.
//-----------------------------------------------------------------------------
#ifndef WORLD_D4C_H_
#define WORLD_D4C_H_

#include "world/macrodefinitions.h"

#if 1
namespace sptk {
namespace world {
#endif

WORLD_BEGIN_C_DECLS

typedef struct {
  double threshold;
} D4COption;

void D4C(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    int fft_size, const D4COption *option, double **aperiodicity);

void InitializeD4COption(D4COption *option);

WORLD_END_C_DECLS

#if 1
}  // namespace world
}  // namespace sptk
#endif

#endif  // WORLD_D4C_H_
