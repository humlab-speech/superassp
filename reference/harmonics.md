# Compute the harmonic frequency structure from f0 measurements

This function takes a pre-computed f0 track and derive `n` harmonic
tracks from it so that the vector of f0 values are now a matrix with `n`
columns. Each column then encode the `n`th harmonic values.

## Usage

``` r
harmonics(track, column = "f0", n = 5, explicitExt = "har", toFile = TRUE)
```

## Arguments

- track:

  An f0 track, either as an SSFF object or as the name of an SSFF
  formatted file. It is recommended that the

- column:

  The name of the column to use as an f0 track.

- n:

  The number of harmonics to compute.

- explicitExt:

  The output file extension.

- toFile:

  boolean;Should the SSFF track be returned or stored on disc?

## Value

The SSFF track object, if required. `NULL` otherwise.

## Details

The stored harmonic frequencies are simply multiples of the fundamental
frequency (f0) track, and not derived independently from the speech
signal. Therefore, errors in the frequency tracking of the f0 signal
will be carried over to these tracks. The primary use case for the track
is to have have estimates of the harmonic frequencies to visualize
harmonic frequency (n\*f~0~) against harmonic amplitude (L~1-n~) .
