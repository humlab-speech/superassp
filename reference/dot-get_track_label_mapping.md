# Get track label mapping

Returns a list mapping column names to display labels.

## Usage

``` r
.get_track_label_mapping()
```

## Value

Named list with 'short' and 'full' label mappings.

## Details

Provides two levels of labels:

- **short**: For plot axes (e.g., "fo \[Hz\]", "F1 \[Hz\]")

- **full**: For documentation/papers (e.g., "Frequency of oscillation
  \[Hz\]")
