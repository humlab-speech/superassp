//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// Adapted for superassp/SPTK integration:
// All symbols wrapped in namespace sptk::world:: to match SPTK's existing
// WORLD namespace convention (same as world/common.h, world/dio.h etc.).
//-----------------------------------------------------------------------------
#ifndef WORLD_CHEAPTRICK_H_
#define WORLD_CHEAPTRICK_H_

#include "world/macrodefinitions.h"

#if 1
namespace sptk {
namespace world {
#endif

WORLD_BEGIN_C_DECLS

typedef struct {
  double q1;
  double f0_floor;
  int fft_size;
} CheapTrickOption;

void CheapTrick(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    const CheapTrickOption *option, double **spectrogram);

void InitializeCheapTrickOption(int fs, CheapTrickOption *option);

int GetFFTSizeForCheapTrick(int fs, const CheapTrickOption *option);

double GetF0FloorForCheapTrick(int fs, int fft_size);

WORLD_END_C_DECLS

#if 1
}  // namespace world
}  // namespace sptk
#endif

#endif  // WORLD_CHEAPTRICK_H_
