# Track Maxima Dispersion Quotient (MDQ) for breathy/tense voice discrimination

Returns MDQ (Kane and Gobl 2013) per frame at 100 Hz — higher values
indicate more dispersed wavelet maxima (breathier voice).

## Usage

``` r
trk_mdq_vat(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

## Value

If `toFile = FALSE`: an `AsspDataObj` with one track:

- `mdq`:

  REAL32, MDQ value, n_frames × 1. Higher values indicate more dispersed
  wavelet maxima (breathier voice).

Frame rate: 100 Hz (10 ms hop).

## Details

MDQ is computed per-GCI using the bit-faithful Daless wavelet bank from
voiceanalysis, then resampled to a fixed 100 Hz frame grid (10 ms hop)
so it can sit alongside other `trk_*` tracks.

## References

(Kane and Gobl 2013)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_mdq_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
