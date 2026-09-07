# Get ISO 226 Parameters for a Frequency

Retrieves or interpolates the ISO 226:2023 parameters (alpha_f, L_U,
T_f) for a given frequency. Uses linear interpolation (on log-frequency
scale) for frequencies between the standard 1/3-octave values.

## Usage

``` r
.get_iso226_params(freq_hz)
```

## Arguments

- freq_hz:

  Numeric; frequency in Hz (20 to 12500 Hz)

## Value

Named list with elements: alpha_f, L_U, T_f, T_r (reference threshold)
