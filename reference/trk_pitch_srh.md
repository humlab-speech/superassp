# Track fundamental frequency using the Summation of Residual Harmonics (SRH)

Extracts F0 and a voiced/unvoiced decision using SRH (Drugman & Alwan
2011), a two-pass harmonic-summation pitch estimator operating on the
LPC residual. SRH is robust in noisy conditions and produces an
integrated VAD decision. Audio is resampled to 16 kHz internally. The
fixed 10 ms hop differs from RAPT/SWIPE, which honour the `windowShift`
parameter.

## Usage

``` r
trk_pitch_srh(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `f0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

- `vad`:

  REAL32, voiced/unvoiced decision (0 = unvoiced, 1 = voiced), n_frames
  × 1.

Frame rate: 100 Hz (fixed 10 ms hop). If `toFile = TRUE`: integer count
of files written, returned invisibly.

## References

(Drugman and Alwan 2011)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract F0 using SRH
trk_pitch_srh("recording.wav")

# Process with custom F0 range
trk_pitch_srh("speech.wav", minF = 80, maxF = 300)
} # }
```
