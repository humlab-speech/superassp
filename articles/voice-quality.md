# Voice quality and creak

Voice quality analysis in `superassp` spans three layers: glottal-source
extraction (`trk_*` functions returning a per-sample or per-frame
track), creak/GCI detection (`trk_*` functions returning event tracks),
and whole-recording summary scores (`lst_*` functions).

## Glottal flow: Iterative Adaptive Inverse Filtering

[`trk_covarep_iaif()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md)
separates the glottal source from the vocal tract, returning the glottal
flow and its derivative (MFDR) at the native audio sample rate — the
standard input to jitter/shimmer/HNR-style measures.

``` r

wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

iaif <- trk_covarep_iaif(wav, toFile = FALSE, verbose = FALSE)
track_names(iaif)
#> [1] "glottal_flow"       "glottal_derivative"
```

## Detecting creaky voice

Two creak detectors are available with the same output schema
(`creak_pp`: posterior probability; `creak_bin`: binary decision), so
they are interchangeable in downstream code:

- [`trk_creak_vat()`](https://humlab-speech.github.io/superassp/reference/trk_creak_vat.md)
  — a bit-faithful Rcpp port of the Kane-Drugman + Ishi ANN pipeline
  from the MATLAB Voice Analysis Toolkit. Prefer this when parity with
  the original MATLAB VAT output matters.
- [`trk_covarep_creak()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)
  — a partial R-side reimplementation; may differ from
  [`trk_creak_vat()`](https://humlab-speech.github.io/superassp/reference/trk_creak_vat.md)
  on borderline frames.

``` r

creak <- trk_creak_vat(wav, toFile = FALSE)
#> Applying `trk_creak_vat()` to 1 recording
head(as.data.frame(creak))
#> Warning: Package 'units' not available. Skipping unit assignment.
#>   frame_time    creak_pp creak_bin
#> 1       0.00 0.018787402         0
#> 2       0.01 0.028024589         0
#> 3       0.02 0.031346901         0
#> 4       0.03 0.031346901         0
#> 5       0.04 0.005675122         0
#> 6       0.05 0.005072619         0
```

Glottal closure instants (GCIs) — the time points used to anchor many
voice-quality measures — are available separately via
[`trk_gci_vat()`](https://humlab-speech.github.io/superassp/reference/trk_gci_vat.md).

## Whole-recording voice quality summaries

[`lst_vq()`](https://humlab-speech.github.io/superassp/reference/lst_vq.md)
runs the Praat-based `VQ_measurements_V2.praat` pipeline (two-pass
adaptive pitch detection followed by 36 voice-quality parameters:
jitter, shimmer, multi-band HNR, spectral energy measures, and
glottal-to-noise excitation ratio) and returns one row per file.

``` r

vq <- lst_vq(wav, toFile = FALSE, verbose = FALSE)
str(vq)
#> 'data.frame':    1 obs. of  39 variables:
#>  $ file_name            : chr "a1"
#>  $ start                : num 0
#>  $ end                  : num 4.04
#>  $ duration             : num 4.04
#>  $ duration_msec        : num 4035
#>  $ mean_period          : num 0.00832
#>  $ sd_period            : num 0.000156
#>  $ jitter_local_percent : num 0.0051
#>  $ jitter_local_abs_db  : num 4.24e-05
#>  $ jitter_rap_percent   : num 0.00264
#>  $ jitter_ppq5_percent  : num 0.00241
#>  $ jitter_ddp_percent   : num 0.00791
#>  $ shimmer_local_percent: num 0.0419
#>  $ shimmer_local_db     : num 0.338
#>  $ shimmer_apq3_percent : num 0.0179
#>  $ shimmer_apq5_percent : num 0.0251
#>  $ shimmer_apq11_percent: num 0.0418
#>  $ shimmer_dda_percent  : num 0.0537
#>  $ hnr_mean_full_db     : num 23.1
#>  $ hnr_sd_full_db       : num 5.1
#>  $ hnr_mean_500_db      : num 26
#>  $ hnr_sd_500_db        : num 6.08
#>  $ hnr_mean_1500_db     : num 23.3
#>  $ hnr_sd_1500_db       : num 5.03
#>  $ hnr_mean_2500_db     : num 23.2
#>  $ hnr_sd_2500_db       : num 5.13
#>  $ hnr_mean_3500_db     : num 23.1
#>  $ hnr_sd_3500_db       : num 4.97
#>  $ energy_1000_db       : num -10.4
#>  $ energy_2000_db       : num -28.9
#>  $ energy_4000_db       : num -10.4
#>  $ energy_6000_db       : num -49.2
#>  $ hammarberg_index_db  : num -29
#>  $ slope_db             : num -1.25
#>  $ tilt_db              : num -2.49
#>  $ bed_db               : num 37.4
#>  $ gne_3500             : num 0.898
#>  $ gne_4500             : num 0.898
#>  $ cpp_db               : num 17.7
```

Restrict the initial pitch search range for atypical voices
(e.g. children, high-pitched pathological voices) via
`minPitchInitial`/`maxPitchInitial`:

``` r

lst_vq("child.wav", minPitchInitial = 150, maxPitchInitial = 600)
```

## Which function should I use?

| Goal | Function | Requires |
|----|----|----|
| Glottal flow / MFDR at sample rate | [`trk_covarep_iaif()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md) | — |
| Creak probability, MATLAB-VAT parity | [`trk_creak_vat()`](https://humlab-speech.github.io/superassp/reference/trk_creak_vat.md) | — |
| Creak probability, lighter R-side path | [`trk_covarep_creak()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md) | — |
| Glottal closure instants | [`trk_gci_vat()`](https://humlab-speech.github.io/superassp/reference/trk_gci_vat.md) | — |
| 36-parameter voice-quality summary | [`lst_vq()`](https://humlab-speech.github.io/superassp/reference/lst_vq.md) | `pladdrr` |
