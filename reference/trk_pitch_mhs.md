# Track pitch using the Modified Harmonic Sieve algorithm

Estimates pitch in Hz using Michel Scheffers' Modified Harmonic Sieve
(MHS) algorithm implemented in the *libassp* C library (Scheffers 2012)
. MHS operates in the frequency domain and is robust to noise; it is a
complementary alternative to the waveform-based `trk_pitch_ksv`.

## Usage

``` r
trk_pitch_mhs(
  listOfFiles,
  beginTime = 0,
  centerTime = FALSE,
  endTime = 0,
  windowShift = 5,
  gender = "u",
  maxF = 600,
  minF = 50,
  minAmp = 50,
  minAC1 = 0.25,
  minRMS = 18,
  maxZCR = 3000,
  minProb = 0.52,
  plainSpectrum = FALSE,
  toFile = FALSE,
  explicitExt = "pit",
  outputDirectory = NULL,
  assertLossless = NULL,
  logToFile = FALSE,
  convertOverwrites = FALSE,
  keepConverted = FALSE,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- centerTime:

  Numeric or logical. Single-frame analysis time point in seconds;
  overrides `beginTime`, `endTime`, and `windowShift`. Default `FALSE`.

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5 ms.

- gender:

  Character. Gender-specific pitch search range: `"f"` (female), `"m"`
  (male), `"u"` (unknown, default).

- maxF:

  Numeric. Maximum pitch in Hz. Default 600.0.

- minF:

  Numeric. Minimum pitch in Hz. Default 50.0.

- minAmp:

  Numeric. Minimum signal amplitude threshold. Default 50.0.

- minAC1:

  Numeric. Minimum first autocorrelation coefficient. Default 0.25.

- minRMS:

  Numeric. Minimum RMS amplitude in dB for voiced detection. Default
  18.0.

- maxZCR:

  Numeric. Maximum zero-crossing rate in Hz for voiced detection.
  Default 3000.0.

- minProb:

  Numeric. Minimum harmonic sieve fit quality (0–1) for a frame to be
  considered voiced. Default 0.52.

- plainSpectrum:

  Logical. Use plain (non-pre-emphasised) spectrum. Default `FALSE`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"pit"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- assertLossless:

  Character vector of additional file extensions to treat as losslessly
  encoded.

- logToFile:

  Logical. Write processing log to a file in `outputDirectory` rather
  than the console. Default `FALSE`.

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- verbose:

  Logical. Print per-file progress. Default `FALSE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `pitch[Hz]`:

  REAL32, Hz, n_frames x 1 column. Estimated pitch frequency; 0
  indicates unvoiced frames.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

Voicing is determined by joint thresholds on `minAmp`, `minAC1`,
`minRMS`, `maxZCR`, and `minProb`. Increase `minProb` to reduce false
voiced decisions in noisy conditions.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

[wrassp::mhsF0](https://rdrr.io/pkg/wrassp/man/mhsF0.html)

[trk_pitch_ksv](https://humlab-speech.github.io/superassp/reference/trk_pitch_ksv.md)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(
   system.file("samples", "sustained", package = "superassp"),
   pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate short-term autocorrelation
res <- trk_pitch_mhs(path2wav, toFile=FALSE)
#> Applying `method(trk_pitch_mhs, class_character)()` to 1 recording

# plot fundamental frequency contour
plot(seq(0, n_records(res) - 1) / sample_rate(res) +
       attr(res, 'startTime'),
     res[["pitch[Hz]"]],
     type='l',
     xlab='time (s)',
     ylab="Pitch (Hz)")

```
