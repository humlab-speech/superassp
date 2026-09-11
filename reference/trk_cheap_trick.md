# CheapTrick Spectral Envelope Estimation (WORLD vocoder, C++ implementation)

Estimate the spectral envelope using the CheapTrick algorithm from the
WORLD vocoder (Morise 2015). CheapTrick uses a Hanning-windowed,
F0-adaptive spectral analysis to produce a smooth, high-quality power
spectral envelope per frame.

Internally: F0 is extracted via Harvest (WORLD), then fed to CheapTrick.
The output is a multi-column SSFF track (`sp`) with `fft_size/2 + 1`
spectral coefficients per frame, where `fft_size` is determined by `fs`
and `f0_floor` following the WORLD formula.

All input media formats supported by av are accepted.

## Usage

``` r
trk_cheap_trick(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 5,
  minF = 60,
  maxF = 400,
  voicing_threshold = 0.1,
  q1 = -0.15,
  f0_floor = 71,
  toFile = TRUE,
  explicitExt = "sp",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0.

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds. Default 5.0 ms.

- minF:

  Numeric. Minimum F0 in Hz for the Harvest pitch extractor. Default
  60.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz for the Harvest pitch extractor. Default
  400.0 Hz.

- voicing_threshold:

  Numeric. Voicing threshold for Harvest (default 0.1).

- q1:

  Numeric. CheapTrick spectral regularization parameter (default -0.15).

- f0_floor:

  Numeric. Lower F0 bound used to determine FFT size (default 71.0). For
  16 kHz audio, 71 Hz gives fft_size = 2048 (sp_length = 1025).

- toFile:

  Logical. If `TRUE`, write SSFF output and return count invisibly. If
  `FALSE`, return `AsspDataObj`. Default `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"sp"`.

- outputDirectory:

  Character. Output directory. `NULL` (default) writes alongside the
  input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `sp`:

  REAL32, power spectral envelope, n_frames x (fft_size/2+1).

Frame rate: `1000 / windowShift` Hz. If `toFile = TRUE`: integer count
of files written, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Estimate spectral envelope
trk_cheap_trick("recording.wav")

# Return data without writing file
sp <- trk_cheap_trick("recording.wav", toFile = FALSE)
dim(sp$sp)  # n_frames x 1025 for 16 kHz audio

# Process multiple files
trk_cheap_trick(c("file1.wav", "file2.wav"))
} # }
```
