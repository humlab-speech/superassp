# Infer track label from column name

Intelligently infers a display label from a cleaned column name.

## Usage

``` r
.infer_track_label(col, full = FALSE)
```

## Arguments

- col:

  Character. Column name.

- full:

  Logical. If TRUE, generate full descriptive label.

## Value

Character. Display label.

## Details

First checks the predefined mapping from
[`.get_track_label_mapping()`](https://humlab-speech.github.io/superassp/reference/dot-get_track_label_mapping.md).
If not found, applies intelligent inference rules:

- Formants: `F1_Hz` → "F1 \[Hz\]" or "First formant frequency \[Hz\]"

- Bandwidths: `B1_Hz` → "B1 \[Hz\]" or "First formant bandwidth \[Hz\]"

- Generic: `param_unit` → "param unit"
