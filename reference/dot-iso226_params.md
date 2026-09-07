# ISO 226:2023 Parameter Table

Parameters for calculating equal-loudness-level contours according to
ISO 226:2023 Table 1. These are used in the formulas for converting
between sound pressure level and loudness level.

## Usage

``` r
.iso226_params
```

## Format

A data frame with 29 rows (one per 1/3-octave frequency) and 4 columns:

- freq_hz:

  Frequency in Hz (20 to 12500 Hz)

- alpha_f:

  Exponent for loudness perception at frequency f

- L_U:

  Magnitude of linear transfer function normalized at 1000 Hz (dB)

- T_f:

  Threshold of hearing at frequency f (dB)

## Source

ISO 226:2023, Table 1
