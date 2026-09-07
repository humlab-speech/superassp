# Getting started with superassp

`superassp` is a self-contained R package for speech signal processing.
It bundles a large collection of DSP routines (pitch, formants, voice
quality, prosody, …) behind one consistent, `wrassp`-like interface, and
produces output that plugs directly into `emuR`. No external Praat or
`wrassp` installation is required.

## The two function families

Almost everything you call falls into one of two prefixes:

| Prefix | Produces | Output class | Example |
|----|----|----|----|
| `trk_*` | A **time-series track** following the audio (F0, formants, RMS, …) | `AsspDataObj` | [`trk_pitch_rapt()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_rapt.md) |
| `lst_*` | **Summary statistics** (scalars/vectors: jitter, voice-quality scores, …) | `JsonTrackObj` | [`lst_voice_report()`](https://humlab-speech.github.io/superassp/reference/lst_voice_report.md) |

Two more families round things out: `ucnv_*` for unit conversions
(Hz↔︎Bark/Mel/ERB/semitone, dB↔︎phon/sone) and `read_*`/`write_*` for I/O
(`read_audio`, `read_ssff`, `read_jstf`, and their `write_`
counterparts).

## Output modes: in memory vs. to file

Every `trk_*` and `lst_*` function takes a `toFile` argument:

- `toFile = FALSE` — return the result **in memory** (an `AsspDataObj`
  or `JsonTrackObj`) for immediate use in R. Single file only.
- `toFile = TRUE` (the default for most functions) — **write** an SSFF
  track file (`trk_*`) or JSTF file (`lst_*`) next to each input, and
  return a count. This is the batch mode; pass a vector of paths and it
  parallelises across files.

``` r

wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# In-memory analysis
rms <- trk_rms(wav, toFile = FALSE, verbose = FALSE)
rms
#> In-memory Assp Data Object
#> Format: SSFF (binary)
#> 805 records at 199.547511312217 Hz
#> Duration: 4.034127 s
#> Number of tracks: 1 
#>   RMS[dB] (1 fields)
```

## Inspecting an `AsspDataObj`

Use the accessor generics rather than digging into attributes:

``` r

track_names(rms)       # which tracks are present
#> [1] "RMS[dB]"
sample_rate(rms)       # frame rate in Hz
#> [1] 199.5475
n_records(rms)         # number of analysis frames
#> [1] 805
signal_duration(rms)   # seconds
#> [1] 4.034127
```

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (or
`as_tibble()`) flattens all tracks into one time-indexed table — the
usual bridge to `dplyr` and plotting:

``` r

df <- as.data.frame(rms)
head(df)
#>    frame_time   RMS_dB
#> 1 0.002505669 26.59310
#> 2 0.007517007 30.23827
#> 3 0.012528345 31.27690
#> 4 0.017539683 30.04947
#> 5 0.022551020 26.95653
#> 6 0.027562358 23.09678
```

Individual track matrices are reached by name, e.g. `rms[["RMS[dB]"]]`.

## Inspecting a `JsonTrackObj`

Summary (`lst_*`) functions return a `JsonTrackObj`: a self-describing
container of measures with their schema and provenance.

``` r

vq <- lst_voice_report(wav, toFile = FALSE)
vq                       # compact summary
names(vq$field_schema)   # available measures
vq$slices[[1]]$values    # values for the first slice
```

## Which pitch tracker should I use?

All pitch trackers emit a fundamental-frequency track (0 marks unvoiced
frames). They differ in method, speed, and robustness:

| Function | Method / origin | Notes |
|----|----|----|
| `trk_pitch_rapt` | RAPT autocorrelation (SPTK) | Fast, robust general-purpose default |
| `trk_pitch_swipe` | SWIPE′ (SPTK) | Accurate, good on noisy speech |
| `trk_pitch_reaper` | REAPER (Google) | Also yields voicing / GCI information |
| `trk_pitch_yin` | YIN | Classic difference-function tracker; +probability |
| `trk_pitch_pyin` | probabilistic YIN | Smoother voicing decisions than YIN |
| `trk_pitch_crepe` | CREPE deep neural net (ONNX) | Most accurate; heavier (downloads a model) |
| `trk_pitch_srh` | Summation of Residual Harmonics | Strong under additive noise |
| `trk_pitch_dio` / `trk_pitch_harvest` | DIO / Harvest (WORLD) | Fast (DIO) vs. accurate (Harvest) |
| `trk_pitch_ac` / `trk_pitch_cc` / `trk_pitch_shs` / `trk_pitch_spinet` | Praat (autocorrelation / cross-correlation / subharmonic / SPINET) | Praat-faithful behaviour |
| `trk_pitch_ksv` / `trk_pitch_mhs` | ASSP (K. Schäfer-Vincent / Michel) | Bundled C library, no external deps |

Start with `trk_pitch_rapt` (fast, reliable). Reach for
`trk_pitch_swipe` or `trk_pitch_crepe` when you need maximum accuracy,
and `trk_pitch_srh` for noisy recordings.

## Which formant tracker should I use?

| Function | Method / origin | Notes |
|----|----|----|
| `trk_formant_burg` | Praat Burg LPC (via pladdrr) | Praat-faithful default |
| `trk_formant_forest` | ASSP forest tracker | Bundled C, no external deps |
| `trk_formant_snack` | Snack ESPS-style LPC | Familiar to Snack/WaveSurfer users |
| `trk_formant_cgdzp` | Complex-group-delay (COVAREP) | Robust bandwidth estimates |
| `trk_formant_tvwlp` | Time-varying weighted LP | Better on rapid transitions |
| `trk_formant_deepformants` / `trk_formant_formantnet` | Deep neural nets (ONNX) | Highest accuracy; heavier |

Start with `trk_formant_burg` or `trk_formant_forest`. Use the neural
trackers (`trk_formant_deepformants`, `trk_formant_formantnet`) when
accuracy matters more than runtime.

## Batch processing to files

``` r

files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")

# Writes one .rap SSFF file next to each input, parallelised across files
trk_pitch_rapt(files, toFile = TRUE)
```

The written SSFF/JSTF files load straight back with
[`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md)
/
[`read_jstf()`](https://humlab-speech.github.io/superassp/reference/read_jstf.md)
and are ready for use with `emuR`.
