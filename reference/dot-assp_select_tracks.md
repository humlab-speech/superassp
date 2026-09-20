# Select track columns by name

Accepts the cleaned column name (`"RMS_dB"`), the track template
(`"RMS[dB]"`) or the object's track name.

## Usage

``` r
.assp_select_tracks(col, band, tracks, available)
```

## Arguments

- col, band:

  Character vectors. Candidate column and band names.

- tracks:

  Character vector. Requested names.

- available:

  Character vector. All candidate names, for error messages.

## Value

Logical vector selecting `col`.
