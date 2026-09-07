# Detect if track name has placeholder 'i'

Checks if a track name template contains the placeholder 'i' using the
uniform pattern: uppercase letter + 'i' + (bracket or end-of-string).

## Usage

``` r
.has_placeholder(name)
```

## Arguments

- name:

  Character. Track name template to check.

## Value

Logical. TRUE if name contains placeholder, FALSE otherwise.

## Details

The pattern `[A-Z]i(\\[|$)` matches:

- Fi\[Hz\] → TRUE (formant frequency template)

- Bi\[Hz\] → TRUE (bandwidth template)

- LPCi → TRUE (LP coefficient template)

- ARFi → TRUE (ARF coefficient template)

- Hi\[dB\] → TRUE (harmonic template)

- Ai\[dB\] → TRUE (amplitude template)

But excludes:

- foi\[Hz\] → FALSE (lowercase before i)

- intensity → FALSE (lowercase, i not at end)

- pitch\[Hz\] → FALSE (no i before bracket)

- gain\[dB\] → FALSE (no i)

## Examples

``` r
if (FALSE) { # \dontrun{
.has_placeholder("Fi[Hz]")     # TRUE
.has_placeholder("LPCi")       # TRUE
.has_placeholder("fo[Hz]")     # FALSE
.has_placeholder("intensity")  # FALSE
} # }
```
