# Track short-term cepstral coefficients

Computes the short-term real cepstrum of audio signals via FFT using the
*libassp* C library (Scheffers 2012; Oppenheim and Schafer 2004;
Childers et al. 1977) . Output coefficients index quefrency from 0 to
FFT_length/2 in steps of 0.5 ms (at the default 40 Hz resolution).
Useful for pitch-period detection and spectral tilt estimation.

## Usage

``` r
trk_cepstrum(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- resolution:

  Numeric. Target FFT frequency resolution in Hz; the FFT length is set
  to the smallest power-of-2 meeting this target. Default 40.0.

- fftLength:

  Integer. Explicit FFT length in points; overrides `resolution`.
  Default 0 (use `resolution`).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `C[dB]`:

  REAL32, dB amplitude, n_frames x (FFT_length/2 + 1) columns. Each
  column corresponds to a quefrency of col_index × 0.5 ms.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

The number of coefficients per frame equals FFT_length / 2 + 1
(one-sided, not mirrored). Use `fftLength` to fix the transform size;
otherwise `resolution` governs it. Pre-emphasis is not applied by this
function; apply `trk_afdiff` first if needed.

## References

Childers DG, Skinner DP, Kemerait RC (1977). “The cepstrum: A guide to
processing.” *Proceedings of the IEEE*, **65**(10), 1428–1443. ISSN
0018-9219.
[doi:10.1109/proc.1977.10747](https://doi.org/10.1109/proc.1977.10747)
.  
  
Oppenheim AV, Schafer RW (2004). “From frequency to quefrency: a history
of the cepstrum.” *IEEE Signal Processing Magazine*, **21**(5), 95–106.
ISSN 1053-5888.
[doi:10.1109/msp.2004.1328092](https://doi.org/10.1109/msp.2004.1328092)
.  
  
Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

wrassp::cepstrum

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

# calulate cepstrum
res <- trk_cepstrum(path2wav, toFile=FALSE)
#> Applying `method(trk_cepstrum, class_character)()` to 1 recording

# plot cepstral values at midpoint of signal
plot(y=res[["C[dB]"]][400,],
    x=seq(1,ncol(res[["C[dB]"]])),
    type='l',
    xlab='Quefrency (ms)',
    ylab='Amplitude (dB)')

```
