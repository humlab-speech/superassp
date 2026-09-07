# Track short-term RMS amplitude

Computes the short-term root-mean-square (RMS) amplitude of audio
signals using the *libassp* C library (Scheffers 2012) . By default,
output is expressed in dB (short-term power). A fast, low-memory energy
measure suitable as a voicing or loudness feature.

## Usage

``` r
trk_rms(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- linear:

  Logical. If `TRUE`, return linear RMS amplitude instead of dB. Default
  `FALSE` (dB scale).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `RMS[dB]`:

  REAL32, dB (or linear amplitude if `linear = TRUE`), n_frames x 1
  column. Short-term RMS amplitude per frame.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

Window function defaults to HAMMING (not BLACKMAN as in other libassp
functions). Set `linear = TRUE` for raw amplitude values without log
conversion.

## See also

wrassp::rmsana
