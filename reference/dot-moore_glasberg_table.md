# Moore-Glasberg Lookup Table (ISO 532-2)

Reference table for sone-phon conversion using the Moore-Glasberg
method. Contains 23 calibration points from 0 to 120 phon.

## Usage

``` r
.moore_glasberg_table
```

## Format

A data frame with 23 rows and 2 columns:

- phon:

  Loudness level in phons (0 to 120)

- sone:

  Loudness in sones (0.001 to 337.6)

## References

MATLAB Audio Toolbox documentation for sone2phon (ISO 532-2 reference
table). (The MathWorks, Inc. 2024)
