# Track cepstrally-smoothed spectrum

Computes a short-term power spectrum smoothed by liftering in the
cepstral domain, using the *libassp* C library (Scheffers 2012) .
Cepstral smoothing suppresses spectral fine structure (harmonics),
leaving the vocal tract envelope; prefer `trk_css_spectrum` over
`trk_dft_spectrum` when a smooth spectral envelope is needed without LP
assumptions.

## Usage

``` r
trk_css_spectrum(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- numCeps:

  Integer. Number of cepstral coefficients used for liftering. Default 0
  sets this to sample rate in kHz + 1 (minimum 2).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `CSS[dB]`:

  REAL32, dB, n_frames x (FFT_length/2 + 1) columns. Cepstrally-smoothed
  power spectral amplitude from 0 Hz to the Nyquist rate.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

`numCeps` controls the number of cepstral terms retained before
back-transforming; lower values give a smoother envelope. The FFT size
is governed by `resolution` (or overridden by `fftLength`).

## See also

wrassp::cssSpectrum

[AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)

[av::av_audio_convert](https://docs.ropensci.org/av//reference/encoding.html)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate cepstrally smoothed spectrum
res <- trk_css_spectrum(path2wav, toFile=FALSE)
#> Applying `method(trk_css_spectrum, class_character)()` to 1 recording
resolution <- attr(res,"origFreq") / ncol(res[[1]])

# plot spectral values at midpoint of signal
plot(y=res[["CSS[dB]"]][400,],
    x=seq(1,ncol(res[[1]]),1)* resolution,
    type='l',
    xlab='Frequency (Hz)',
    ylab='Amplitude (dB)')

```
