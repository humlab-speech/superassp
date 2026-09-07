# Track short-term zero-crossing rate

Computes the average of the short-term positive and negative
zero-crossing rates of audio signals using the *libassp* C library
(Scheffers 2012) . ZCR is a simple, fast measure correlated with
spectral centroid and useful for voicing detection and fricative
classification.

## Usage

``` r
trk_zcr(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowSize:

  Numeric. Analysis window size in milliseconds. Default 25 ms.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `ZCR[Hz]`:

  REAL32, Hz, n_frames x 1 column. Average zero-crossing rate (positive
  and negative crossings combined) per frame, expressed as a rate in Hz.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

The ZCR is reported in Hz (crossings per second), averaged over the
positive and negative zero-crossing rates within each analysis window.

## See also

wrassp::zcrana

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate zcr values
res <- trk_zcr(path2wav, toFile=FALSE)
#> Applying `method(trk_zcr, class_character)()` to 1 recording

# plot zcr values
plot(seq(0, n_records(res) - 1) / sample_rate(res) +
      attr(res, 'startTime'),
    res[["ZCR[Hz]"]],
    type='l',
    xlab='time (s)',
    ylab='Zero Crossing Rates (Hz)')

```
