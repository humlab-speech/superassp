# Abort with the standard "pladdrr is missing" error

Every pladdrr-backed function calls this when
[`pladdrr_available()`](https://humlab-speech.github.io/superassp/reference/pladdrr_available.md)
is FALSE, so that a missing optional dependency produces one consistent,
actionable message instead of a per-function variant.

## Usage

``` r
pladdrr_unavailable()
```
