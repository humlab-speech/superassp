# TVWLP Rcpp vs MATLAB Validation Report

Generated: 2026-05-15 14:29:19
Test audio: `test-dr8-mjln0-sx279.wav` (16 kHz, 3.59 s TIMIT sentence)

| LP method | Formant | n_valid | MATLAB mean Hz | R mean Hz | corr | RMSE Hz | MAE Hz | bias Hz | within 50 Hz |
|---|---|---|---|---|---|---|---|---|---|
| tvwlp_l2 (weighted) | F1 | 358/358 | 566.9 | 548.9 | 0.4626 | 276.1 | 131.6 | -18.0 | 61.5% |
| tvwlp_l2 (weighted) | F2 | 358/358 | 1530.9 | 1556.4 | 0.7512 | 287.0 | 130.2 | +25.5 | 53.1% |
| tvwlp_l2 (weighted) | F3 | 358/358 | 2597.6 | 2612.7 | 0.6942 | 261.7 | 141.2 | +15.1 | 53.9% |
| tvlp_l2 (unweighted) | F1 | 358/358 | 611.5 | 567.1 | 0.3374 | 263.0 | 122.2 | -44.4 | 63.4% |
| tvlp_l2 (unweighted) | F2 | 358/358 | 1529.3 | 1521.0 | 0.8731 | 199.8 | 100.3 | -8.3 | 58.4% |
| tvlp_l2 (unweighted) | F3 | 358/358 | 2586.4 | 2615.3 | 0.6556 | 302.7 | 162.9 | +28.8 | 50.8% |

## Interpretation

Post-fix state (May 2026). Three bugs fixed against MATLAB reference:
1. **Levinson-Durbin sign error** in `levinson_durbin()` (tvwlp_core.cpp):
   reflection coefficient + predictor update had wrong signs. LPC residual
   correlation jumped from 0.16 to 0.998 against MATLAB `lpc()`.
2. **Wrong Blackman window** in SRH and SEDREAMS: `av::blackman()` is the
   periodic DFT-even variant; MATLAB uses the symmetric form. GCI count
   dropped from 581 (231 spurious) to 349, matching MATLAB's 350.
3. **Wrong Hann window** convention in get_lpc_residual_cpp (minor).

Remaining gap vs Python port (0.85-0.95 corr): tvwlp_l2 F1=0.46/F2=0.75/F3=0.69.
Probable cause is Armadillo `arma::solve` behaviour on near-singular per-frame
LP systems (the runtime emits `solve(): system is singular` warnings). F1 in
tvlp_l2 (unweighted) stays at 0.34 even though signal correlation is 0.9999,
suggesting the residual gap is in the LP solver, not the preprocessing.

## Caveats

- MATLAB does not output bandwidths; superassp's `bw` track cannot be ground-truthed here.
- R produces 358 frames vs MATLAB 359 (one-frame boundary difference).
- FFmpeg vs MATLAB resampling: signal correlation 0.9999, negligible contribution.
