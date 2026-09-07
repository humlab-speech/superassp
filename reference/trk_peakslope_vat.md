# Peak slope via voiceanalysis Daless wavelet bank

Computes the peak slope acoustic parameter (Kane and Gobl 2011) using
the bit-faithful Daless wavelet bank in voiceanalysis. Faster and more
numerically stable than the pure-R
[`trk_peakslope`](https://humlab-speech.github.io/superassp/reference/trk_peakslope.md)
(which approximates Daless via db4 wavelets) for the same algorithm.

## Usage

``` r
trk_peakslope_vat(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

## Value

If `toFile = FALSE`: an `AsspDataObj` with one track:

- `peak_slope`:

  REAL32, peak-slope coefficient per frame, n_frames × 1.

Frame rate: 100 Hz (10 ms hop).

## References

(Kane and Gobl 2011)

## See also

[`trk_peakslope`](https://humlab-speech.github.io/superassp/reference/trk_peakslope.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_peakslope_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
