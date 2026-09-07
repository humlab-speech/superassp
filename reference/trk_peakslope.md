# Track spectral tilt using D'Alessandro PeakSlope (Morlet wavelet)

Returns a per-frame spectral tilt estimate related to breathiness and
voice quality (Henrich et al. 2004) . Captures similar information to
CPP and H1-H2 but via multi-scale wavelet analysis.

## Usage

``` r
trk_peakslope(listOfFiles, beginTime = 0, endTime = 0, toFile = FALSE, explicitExt = "psl", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"psl"`.

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `peakslope`:

  FLOAT, spectral tilt slope (log10(magnitude) per scale unit), n_frames
  × 1. Negative = energy at low scales (bright/tense voice); positive =
  energy at high scales (breathy/relaxed voice).

Frame rate: 100 Hz (10 ms hop; 40 ms analysis window). If
`toFile = TRUE`: character vector of output file paths, returned
invisibly.

## Details

Spectral tilt is estimated by fitting a linear slope across peak
magnitudes at seven Morlet wavelet scales. Lower (more negative) values
indicate a steeper spectral slope, associated with creakier phonation;
higher values indicate breathier or modal phonation.

Seven D'Alessandro Morlet wavelets at octave-spaced scales (2^0 to 2^6)
are convolved with the signal. The maximum magnitude in each 40 ms frame
is taken per scale, log10-transformed, and a linear regression slope is
computed across the seven scale-magnitude pairs.

## References

Henrich N, d'Alessandro C, Castellengo M, Doval B (2004). “On the use of
the derivative of electroglottographic signal for pitch detection.”
*Journal of the Acoustical Society of America*, **115**, 3040–3058.
[doi:10.1121/1.1738025](https://doi.org/10.1121/1.1738025) . Spectral
tilt and energy distribution measures via Morlet wavelet.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single file, return object
peakslope <- trk_peakslope("speech.wav", toFile = FALSE)

# Batch process multiple files
files <- c("file1.wav", "file2.wav")
trk_peakslope(files, toFile = TRUE, outputDirectory = "output/")
} # }
```
