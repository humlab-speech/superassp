# Voice Quality Measurements using pladdrr

Extract comprehensive voice quality measures from audio signals
including jitter, shimmer, harmonics-to-noise ratio (HNR) at multiple
frequency bands, spectral energy measures, glottal-to-noise excitation
ratio (GNE), and cepstral peak prominence (CPP).

## Usage

``` r
lst_vq(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector with path(s) to audio file(s)

- beginTime:

  Numeric. Start time in seconds (default 0, or NULL = 0)

- endTime:

  Numeric. End time in seconds (0 or NULL = end of file)

- minPitchInitial:

  Numeric. Initial minimum pitch for two-pass detection in Hz (default
  50)

- maxPitchInitial:

  Numeric. Initial maximum pitch for two-pass detection in Hz (default
  800)

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- explicitExt:

  Character. File extension for output. Default "vq".

- outputDirectory:

  Character. Output directory path. Default NULL (use input directory).

- verbose:

  Logical. Print progress messages (default TRUE)

## Value

If `toFile=FALSE`, returns data frame (single file) or list of data
frames (multiple files) with columns:

- file_name:

  File name without extension

- start, end:

  Analysis time bounds (seconds)

- duration, duration_msec:

  Duration in seconds and milliseconds

- mean_period, sd_period:

  Mean and SD of pitch period

- jitter_local_percent, jitter_local_abs_db:

  Local jitter measures

- jitter_rap_percent, jitter_ppq5_percent, jitter_ddp_percent:

  Relative jitter measures

- shimmer_local_percent, shimmer_local_db:

  Local shimmer measures

- shimmer_apq3_percent, shimmer_apq5_percent, shimmer_apq11_percent,
  shimmer_dda_percent:

  Amplitude perturbation

- hnr_mean_full_db, hnr_sd_full_db:

  Full-spectrum HNR mean and SD

- hnr_mean_500_db, hnr_sd_500_db:

  HNR 0-500 Hz

- hnr_mean_1500_db, hnr_sd_1500_db:

  HNR 0-1500 Hz

- hnr_mean_2500_db, hnr_sd_2500_db:

  HNR 0-2500 Hz

- hnr_mean_3500_db, hnr_sd_3500_db:

  HNR 0-3500 Hz

- energy_1000_db, energy_2000_db, energy_4000_db, energy_6000_db:

  Band energy differences

- hammarberg_index_db:

  Hammarberg index (0-2kHz vs 2-5kHz)

- slope_db, tilt_db:

  LTAS slope and tilt

- bed_db:

  Band Energy Difference (low vs high)

- gne_3500, gne_4500:

  Glottal-to-Noise Excitation ratio

- cpp_db:

  Cepstral Peak Prominence

If `toFile=TRUE`, invisibly returns path(s) to JSTF file(s).

## Details

This function implements the algorithm from VQ_measurements_V2.praat,
using two-pass adaptive pitch detection for speaker-specific F0 range
estimation, followed by extraction of 36 voice quality parameters.

The function performs a two-pass pitch detection:

1.  Initial pass with wide range (50-800 Hz) to estimate speaker F0

2.  Adaptive range calculated from quartiles (Q1×0.75, Q3×1.5)

3.  Second pass with adaptive range for accurate F0 tracking

Voice quality measures include:

- **Jitter**: Period-to-period variation (local, RAP, PPQ5, DDP)

- **Shimmer**: Amplitude variation (local, APQ3, APQ5, APQ11, DDA)

- **HNR**: Harmonics-to-Noise Ratio at 5 frequency ranges

- **Spectral energy**: Band energy differences and Hammarberg index

- **LTAS**: Long-term average spectrum slope and tilt

- **GNE**: Glottal-to-Noise Excitation ratio

- **CPP**: Cepstral Peak Prominence

All measures are calculated from voiced segments only.

## References

( )

(Boersma and Weenink 2023)

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze single file
vq <- lst_vq("speech.wav")
print(vq)

# Check voice quality
cat(sprintf("Jitter: %.2f%%\n", vq$jitter_local_percent))
cat(sprintf("Shimmer: %.2f%%\n", vq$shimmer_local_percent))
cat(sprintf("HNR: %.1f dB\n", vq$hnr_mean_full_db))
cat(sprintf("CPP: %.1f dB\n", vq$cpp_db))

# Batch analysis with file output
lst_vq(c("f1.wav", "f2.wav", "f3.wav"), toFile = TRUE)

# Custom pitch range for child speech
vq <- lst_vq("child.wav", minPitchInitial = 150, maxPitchInitial = 600)
} # }
```
