# Formant frequencies and bandwidths via Praat's Burg method

Tracks formant frequencies (F1–F5), bandwidths (B1–B5), and optionally
spectral intensities (L1–L5) using Praat's Burg LPC algorithm via
pladdrr. Optional HMM-based formant tracking smooths trajectories across
time. Prefer this over `trk_formant_forest` when Praat-compatible Burg
estimates are needed.

## Usage

``` r
trk_formant_burg(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- timeStep:

  Numeric. Frame shift in seconds; sets output frame rate (1 / timeStep
  Hz). Default 0.005 s (200 Hz).

- number_of_formants:

  Integer. Number of formant candidates to find (up to 5). Default 5.

- maxHzFormant:

  Numeric. Formant ceiling in Hz; typically 5500 for female, 5000 for
  male speakers. Default 5500 Hz.

- windowLength:

  Numeric. LPC analysis window length in seconds. Default 0.025 s.

- pre_emphasis:

  Numeric. Pre-emphasis onset frequency in Hz. Default 50 Hz.

- track_formants:

  Logical. Apply HMM-based formant tracking to smooth trajectories.
  Default `FALSE`.

- number_of_tracks:

  Integer. Number of tracks to retain when `track_formants = TRUE`.
  Default 3.

- reference_F1:

  Numeric. Reference F1 frequency for HMM tracking in Hz. Default 550
  Hz.

- reference_F2:

  Numeric. Reference F2 frequency for HMM tracking in Hz. Default 1650
  Hz.

- reference_F3:

  Numeric. Reference F3 frequency for HMM tracking in Hz. Default 2750
  Hz.

- reference_F4:

  Numeric. Reference F4 frequency for HMM tracking in Hz. Default 3850
  Hz.

- reference_F5:

  Numeric. Reference F5 frequency for HMM tracking in Hz. Default 4950
  Hz.

- frequency_cost:

  Numeric. HMM cost weight for deviation from reference frequencies.
  Default 1.0.

- bandwidth_cost:

  Numeric. HMM cost weight for formant bandwidth. Default 1.0.

- transition_cost:

  Numeric. HMM cost weight for frame-to-frame transitions. Default 1.0.

- windowShape:

  Character. Window shape applied before analysis. Default
  `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the analysis window. Default 1.0.

- include_intensity:

  Logical. Extract spectral intensity at each formant frequency (L1–L5
  tracks via spectrogram lookup). Default `TRUE`.

- spectrogram_resolution:

  Numeric. Frequency resolution of the spectrogram used for intensity
  extraction in Hz. Default 40 Hz.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"pfm"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `fm1`–`fm5`:

  REAL32, Hz, n_frames x 1. Formant frequencies F1–F5. 0 encodes an
  undefined (unvoiced) frame.

- `bw1`–`bw5`:

  REAL32, Hz, n_frames x 1. Formant bandwidths B1–B5. 0 encodes an
  undefined frame.

- `L1`–`L5`:

  REAL32, Pa^2/Hz, n_frames x 1. Spectral power at each formant
  frequency (only when `include_intensity = TRUE`).

Frame rate: `1 / timeStep` Hz (default 200 Hz). If `toFile = TRUE`:
integer count of files written, returned invisibly.

## Details

Uses Praat's Burg all-pole LPC estimator. Spectral intensity tracks
(L1–L5) are extracted from a separate spectrogram and require pladdrr
\>= 4.8.20; earlier versions may crash with `include_intensity = TRUE`.
HMM tracking (`track_formants = TRUE`) may not be available in all
pladdrr builds; a warning is issued and untracked formants are returned.

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze formants for a single file
formants <- trk_formant_burg("speech.wav", toFile = FALSE)

# Batch process multiple files with tracking
files <- c("speech1.wav", "speech2.wav")
trk_formant_burg(files, track_formants = TRUE, toFile = TRUE)

# With custom formant ceiling (e.g., male speaker)
formants <- trk_formant_burg("speech.wav",
                         maxHzFormant = 5000,
                         number_of_formants = 4,
                         toFile = FALSE)
} # }
```
