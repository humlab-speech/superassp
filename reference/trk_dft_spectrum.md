# Track short-term DFT power spectrum

Computes a short-term power spectrum via the Fast Fourier Transform
using the *libassp* C library (Scheffers 2012) . Produces an unsmoothed
narrow-band spectrum from 0 Hz to the Nyquist rate. Prefer this function
when raw spectral detail is needed; use `trk_css_spectrum` or
`trk_lps_spectrum` for smoothed spectral envelopes.

## Usage

``` r
trk_dft_spectrum(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- bandwidth:

  Numeric. Effective analysis bandwidth in Hz. Default 0 yields the
  minimum bandwidth determined by the FFT length.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `DFT[dB]`:

  REAL32, dB power, n_frames x (FFT_length/2 + 1) columns. Power
  spectral amplitude from 0 Hz to the Nyquist rate.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

The FFT length is determined by `resolution` unless overridden by
`fftLength`. `bandwidth` widens the effective analysis window, trading
spectral resolution for reduced side-lobe leakage.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

wrassp::dftSpectrum

[AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)

[av::av_audio_convert](https://docs.ropensci.org/av//reference/encoding.html)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(system.file("extdata", package = "wrassp"),
                       pattern = glob2rx("*.wav"),
                       full.names = TRUE)[1]

# calculate dft spectrum
res <- trk_dft_spectrum(path2wav, toFile=FALSE)
#> Applying `method(trk_dft_spectrum, class_character)()` to 1 recording
#> Warning: path[1]="NA": No such file or directory
#> Warning: ! Found 1 recording in lossy format
#> ℹ Lossy compression may affect `spectrum()` accuracy
#> ✖ For accurate DSP, use lossless formats: "wav", "au", "kay", "nist", and "nsp"
#> Error in av_to_asspDataObj(file_path, start_time = bt, end_time = if (et ==     0) NULL else et, target_sample_rate = NULL): path[1]="NA": No such file or directory

# plot spectral values at midpoint of signal
plot(res$dft[dim(res$dft)[1]/2,],
     type='l',
     xlab='spectral value index',
     ylab='spectral value')
#> Error: object 'res' not found
```
