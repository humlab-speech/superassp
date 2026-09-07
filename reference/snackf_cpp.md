# Snack Formant Extraction (C++ implementation)

LPC + dynamic-programming formant tracker from the Snack Sound Toolkit
(Talkin / AT&T / KTH). Returns formant frequencies and bandwidths.

## Usage

``` r
snackf_cpp(
  audio_obj,
  numFormants = 4L,
  lpcOrder = 12L,
  windowLength = 0.049,
  windowShift = 10,
  preEmphasis = 0.7,
  dsFreq = 10000,
  nomF1 = -10,
  lpcType = 0L,
  windowType = 2L,
  verbose = FALSE
)
```

## Arguments

- audio_obj:

  An AsspDataObj containing audio data

- numFormants:

  Number of formants to track (default 4, max 7)

- lpcOrder:

  LPC order (default 12)

- windowLength:

  Analysis window duration in seconds (default 0.049)

- windowShift:

  Frame shift in milliseconds (default 10)

- preEmphasis:

  Pre-emphasis factor (default 0.7)

- dsFreq:

  Downsample target frequency in Hz (default 10000)

- nomF1:

  Nominal F1 for DP cost (default -10 = use defaults)

- lpcType:

  LPC method: 0=autocorrelation, 1=stabilized covariance, 2=covariance
  (default 0)

- windowType:

  Window type: 0=rectangular, 1=Hamming, 2=cos^4, 3=Hanning (default 2)

- verbose:

  Print processing info (default FALSE)

## Value

List with fm (frequency matrix), bw (bandwidth matrix), times,
sample_rate, n_frames
