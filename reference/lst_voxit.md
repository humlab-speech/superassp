# Extract Voxit prosodic complexity features from audio files

Computes voice and articulation complexity measures from audio files
using word-level alignments (optional) and WORLD-based pitch/spectral
analysis. Features include speaking rate (WPM), pause statistics, F0
dynamics, pitch complexity, and intensity.

Uses superassp's C++ implementations of WORLD vocoder analysis (Harvest
F0, intensity, LZ complexity, SG smoothing) — no external dependencies
beyond those already in superassp.

## Usage

``` r
lst_voxit(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths

- alignmentFiles:

  Character vector of alignment CSV paths (optional). CSV must have
  columns: `word`, `start`, `end` (in seconds). If NULL, word-based
  metrics (WPM, pause stats, rhythmic complexity) are set to NA;
  signal-based metrics are always computed.

- beginTime:

  Numeric. Start of analysis window in seconds (default: 0)

- endTime:

  Numeric. End of analysis window in seconds (default: 0)

- minF:

  Numeric. Minimum F0 in Hz for Harvest (default: 60)

- maxF:

  Numeric. Maximum F0 in Hz for Harvest (default: 600)

- verbose:

  Logical. Show progress (default: TRUE)

- parallel:

  Logical. Use parallel processing (default: TRUE)

- n_cores:

  Integer. Number of cores (default: auto-detect)

- toFile:

  Logical. Write JSTF files (default: FALSE)

- explicitExt:

  Character. Output extension (default: "vxt")

- outputDirectory:

  Character. Output directory (default: NULL = input dir)

## Value

If `toFile=FALSE` (default): for single file, a named list with up to 18
prosodic features. For multiple files, a list where each element is a
file's feature list. If `toFile=TRUE`, invisibly returns the output file
path(s).

## Features

**Word-level (requires alignment):**

- `WPM`: Words per minute (speaking rate)

- `pause_count`: Number of pauses (100–3000 ms)

- `long_pause_count`: Pauses \> 3 s

- `average_pause_length`: Mean pause duration (s)

- `average_pause_rate`: Pauses per second

- `rhythmic_complexity_of_pauses`: LZ complexity of pause/speech
  sequence

**Signal-level (WORLD-based):**

- `average_pitch`: Mean F0 (Hz) — arithmetic mean of geometric mean per
  voiced segment

- `pitch_range`: F0 range in octaves

- `pitch_speed`: F0 velocity (octaves/s)

- `pitch_acceleration`: F0 acceleration (octaves/s²)

- `pitch_entropy`: Shannon entropy of F0 histogram (bits)

- `f0_geometric_mean_hz`: Geometric mean F0 (Hz)

- `voicing_percent`: Percentage of voiced frames

- `intensity_mean_db`: Mean intensity (dB)

- `lz_complexity_voiced`: LZ complexity of VUV sequence

- `f0_velocity_mean_abs`: Mean absolute F0 velocity (octaves/s)

- `f0_accel_mean_abs`: Mean absolute F0 acceleration (octaves/s²)

- `dynamism`: Composite measure: \|f0MeanAbsVel\| \* f0Entropy + LZ \*
  0.439

## Examples

``` r
if (FALSE) { # \dontrun{
# With alignment (all metrics)
features <- lst_voxit("speech.wav", alignmentFiles = "speech_align.csv")
print(features$WPM)
print(features$f0_geometric_mean_hz)

# Without alignment (signal metrics only)
features <- lst_voxit("speech.wav")
print(features$f0_geometric_mean_hz)  # OK
print(features$WPM)  # NA
} # }
```
