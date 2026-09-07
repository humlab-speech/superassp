# Comprehensive voice quality feature set via PraatSauce

Extracts 37 time-series voice quality measures — F0, formants, harmonic
amplitudes (corrected and uncorrected), HNR in four bands, and CPP —
replicating the VoiceSauce feature set via Praat's algorithms in
pladdrr. Formant corrections follow Iseli & Alwan (2004); bandwidth can
optionally be estimated with Hawks & Miller (1995).

## Usage

``` r
trk_praatsauce(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5 ms (200 Hz).

- windowSize:

  Numeric. Analysis window length in milliseconds. Default 25 ms.

- minF:

  Numeric. Lower F0 bound in Hz. Default 50 Hz.

- maxF:

  Numeric. Upper F0 bound (ceiling) in Hz. Default 300 Hz.

- formantTracking:

  Logical. Attempt HMM formant tracking. Currently unsupported in
  pladdrr — a warning is issued and untracked formants are used. Default
  `TRUE`.

- numFormants:

  Integer. Number of formants to extract (up to 5). Default 5.

- maxFormantHz:

  Numeric. Formant ceiling in Hz. Default 5000 Hz.

- nominalF1:

  Numeric. Reference F1 for formant tracking in Hz. Default 500 Hz.

- nominalF2:

  Numeric. Reference F2 for formant tracking in Hz. Default 1500 Hz.

- nominalF3:

  Numeric. Reference F3 for formant tracking in Hz. Default 2500 Hz.

- preEmphFrom:

  Numeric. Pre-emphasis onset frequency in Hz. Default 50 Hz.

- useBandwidthFormula:

  Logical. If `TRUE`, estimate bandwidths with the Hawks & Miller (1995)
  formula instead of Praat's Burg estimates. Default `FALSE`.

- channel:

  Integer. Audio channel to use for multi-channel files. Default 1.

- resample_to_16k:

  Logical. Resample to 16 kHz before analysis. Default `TRUE`.

- windowShape:

  Character. Window shape for audio extraction. Default `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the extraction window. Default 1.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"psa"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with 37 REAL32 tracks, n_frames x
1 each. Frame rate: `1000 / windowShift` Hz (default 200 Hz).

- `f0`:

  Hz. Fundamental frequency (AC method; 0 = unvoiced).

- `F1`, `F2`, `F3`:

  Hz. Formant frequencies.

- `B1`, `B2`, `B3`:

  Hz. Formant bandwidths.

- `H1u`, `H2u`, `H4u`:

  dB. Uncorrected harmonic amplitudes at 1×, 2×, 4× F0.

- `H2Ku`, `H5Ku`:

  dB. Uncorrected spectral amplitude near 2 kHz and 5 kHz.

- `A1u`, `A2u`, `A3u`:

  dB. Uncorrected amplitude at F1, F2, F3.

- `H1H2u`, `H2H4u`, `H1A1u`, `H1A2u`, `H1A3u`, `H2KH5Ku`:

  dB. Uncorrected spectral slope differences.

- `H1c`, `H2c`, `H4c`:

  dB. Formant-corrected harmonic amplitudes (Iseli & Alwan 2004).

- `A1c`, `A2c`, `A3c`:

  dB. Formant-corrected amplitudes at F1, F2, F3.

- `H1H2c`, `H2H4c`, `H1A1c`, `H1A2c`, `H1A3c`:

  dB. Corrected spectral slope differences.

- `CPP`:

  dB. Cepstral Peak Prominence.

- `HNR05`, `HNR15`, `HNR25`, `HNR35`:

  dB. Harmonics-to-noise ratio in bands 0–500 Hz, 0–1500 Hz, 0–2500 Hz,
  0–3500 Hz.

If `toFile = TRUE`: character vector of output file paths, returned
invisibly.

## Details

Audio is optionally resampled to 16 kHz (`resample_to_16k = TRUE`)
before analysis. Harmonic amplitudes are searched in windows of ±10%
around the expected harmonic frequency. Formant corrections use the
Iseli & Alwan (2004) formula applied to F1 and F2 for H1/H2/H4, and
additionally F3 for A3.

## References

(Iseli and Alwan 2004)

(Hawks and Miller 1995)

(Shue et al. 2011)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with 5ms frame shift
result <- trk_praatsauce("speech.wav", windowShift = 5, toFile = FALSE)

# Access corrected H1-H2 (breathiness measure)
plot(result$H1H2c, type = "l", main = "H1-H2 Corrected")

# Process with custom F0 range and bandwidth formula
result <- trk_praatsauce(
  "speech.wav",
  minF = 75,
  maxF = 300,
  useBandwidthFormula = TRUE,
  toFile = FALSE
)

# Batch process multiple files
trk_praatsauce(c("f1.wav", "f2.wav", "f3.wav"), toFile = TRUE)
} # }
```
