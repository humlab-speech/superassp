# Synthesise an LF model glottal pulse via voiceanalysis

Convenience wrapper around `.vat_rd2r()` and `.vat_lf_cont()` for
generating Liljencrants-Fant (Fant et al. 1985) glottal flow derivative
pulses from a Rd shape descriptor. Useful for source-modelling
experiments and synthesis demos.

## Usage

``` r
lst_lf_vat_synthesis(Rd, F0, fs, EE = 1)
```

## Arguments

- Rd:

  Numeric. LF shape descriptor (typical range 0.3–2.5; smaller =
  tenser).

- F0:

  Numeric. Fundamental frequency in Hz.

- fs:

  Numeric. Sampling frequency in Hz.

- EE:

  Numeric. Excitation strength (default 1.0).

## Value

Named list:

- `pulse`:

  Numeric vector — LF glottal flow derivative sampled at `fs`.

- `Rd`:

  The Rd that was used.

- `Ra`, `Rk`, `Rg`:

  Derived LF R-parameters.

- `F0`, `fs`, `EE`:

  Inputs (for reproducibility).

## Details

This is a pure synthesis utility — it does not analyse an input file.
For per-cycle Rd estimation on real speech, a future `lst_lf_vat_fit()`
will land once `dyProg_LF` is ported.

## References

(Fant et al. 1985)
