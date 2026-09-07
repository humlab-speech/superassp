# Snack Pitch Extraction — full 4-track output (C++ implementation)

Normalized cross-correlation + dynamic-programming pitch tracker from
the Snack Sound Toolkit (Talkin, 1995). Returns F0, voicing probability,
RMS energy and autocorrelation-peak per frame.

## Usage

``` r
snackp_cpp(
  audio_obj,
  minF = 50,
  maxF = 550,
  windowShift = 10,
  voiceBias = 0,
  verbose = FALSE
)
```

## Arguments

- audio_obj:

  An AsspDataObj containing audio data

- minF:

  Minimum F0 in Hz (default 50)

- maxF:

  Maximum F0 in Hz (default 550)

- windowShift:

  Frame shift in milliseconds (default 10)

- voiceBias:

  Bias toward voiced hypothesis (default 0.0)

- verbose:

  Print processing info (default FALSE)

## Value

List with f0, voicing, rms, acpeak (matrices), times, sample_rate,
n_frames
