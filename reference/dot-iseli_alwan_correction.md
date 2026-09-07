# Iseli & Alwan (2004) Formant Correction

Internal helper function implementing Iseli & Alwan formant correction.
Calculates how much a formant affects harmonic amplitude measurement.

## Usage

``` r
.iseli_alwan_correction(f, fx, bx, fs)
```

## Arguments

- f:

  Frequency of harmonic being measured (Hz)

- fx:

  Formant frequency (Hz)

- bx:

  Formant bandwidth (Hz)

- fs:

  Sampling frequency (Hz)

## Value

Correction in dB to subtract from measured harmonic
