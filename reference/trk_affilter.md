# Apply a digital filter to audio signals

Filters audio waveforms using FIR or Butterworth IIR designs implemented
in the *libassp* C library (Scheffers 2012) . Supports high-pass,
low-pass, and band-pass configurations. At least one of `highPass` or
`lowPass` must be specified.

## Usage

``` r
trk_affilter(listOfFiles, highPass = NULL, lowPass = NULL, stopBand = 96, transition = 250, useIIR = FALSE, numIIRsections = 4L, beginTime = 0, endTime = 0, toFile = TRUE, explicitExt = "flt", outputDirectory = NULL, assertLossless = NULL, logToFile = FALSE, keepConverted = FALSE, convertOverwrites = FALSE, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- highPass:

  Numeric. High-pass cutoff frequency in Hz. `NULL` disables high-pass
  filtering. Default `NULL`.

- lowPass:

  Numeric. Low-pass cutoff frequency in Hz. `NULL` disables low-pass
  filtering. Default `NULL`.

- stopBand:

  Numeric. FIR stop-band attenuation in dB (Kaiser window design).
  Default 96.0.

- transition:

  Numeric. FIR transition band width in Hz. Default 250.0.

- useIIR:

  Logical. Use Butterworth IIR filter instead of FIR. Default `FALSE`.

- numIIRsections:

  Integer. Number of 2nd-order IIR sections (filter order). Default
  `4L`.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"flt"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- assertLossless:

  Character vector of additional file extensions to treat as losslessly
  encoded.

- logToFile:

  Logical. Write processing log to a file in `outputDirectory` rather
  than the console. Default `FALSE`.

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track name preserved from
libassp output (the filtered waveform, same label as the input channel),
containing INT16 or REAL32 sample values. If `toFile = TRUE`: integer
count of files written, returned invisibly.

## Details

Filter mode is determined by `highPass`/`lowPass`: supply only
`highPass` for a high-pass filter, only `lowPass` for a low-pass filter,
or both for a band-pass filter. `stopBand` and `transition` govern the
FIR Kaiser-window design; set `useIIR = TRUE` to switch to a Butterworth
IIR design instead.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## Author

Fredrik Nylén

## Examples

``` r
path2wav <- list.files(system.file("samples", "sustained", package = "superassp"),
                       pattern = glob2rx("a1.wav"), full.names = TRUE)
# High-pass filter above 100 Hz
res <- trk_affilter(path2wav, highPass = 100, toFile = FALSE)
#> Applying `method(trk_affilter, class_character)()` to 1 recording
```
