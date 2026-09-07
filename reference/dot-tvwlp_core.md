# TVWLP core pipeline

Pre-emphasize, estimate pitch + GCI (for weighted methods), solve
time-varying LP per frame, extract formants from roots, downsample,
median filter. Expects audio pre-resampled to 8 kHz.

## Usage

``` r
.tvwlp_core(s, fs, lptype, p, q, npeaks, preemp, fint)
```

## Value

List with Fi (n_frames x npeaks), Bw (n_frames x npeaks), n_frames,
frame_rate.
