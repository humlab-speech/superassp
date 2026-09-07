# Parse unit from column name

Extracts the unit suffix from a cleaned column name.

## Usage

``` r
.parse_unit_from_colname(col_name)
```

## Arguments

- col_name:

  Character. Column name to parse.

## Value

Character. Unit name (e.g., "Hz", "dB") or NA if no unit found.

## Details

Handles both track-name conventions used across the package:

- Bracket form (raw track names): `fo[Hz]` → "Hz", `F1[Hz]` → "Hz"

- Cleaned underscore form: `fo_Hz` → "Hz", `H1_H2c_dB` → "dB"

For the underscore form only *known* units are recognised, so ordinary
underscored names are not mistaken for units:

- `frame_time` → NA (`time` is not a unit)

- `fm_1` → NA

- `intensity` → NA (no unit)

## Examples

``` r
if (FALSE) { # \dontrun{
.parse_unit_from_colname("fo[Hz]")       # "Hz"
.parse_unit_from_colname("fo_Hz")        # "Hz"
.parse_unit_from_colname("CPP_dB")       # "dB"
.parse_unit_from_colname("frame_time")   # NA
} # }
```
