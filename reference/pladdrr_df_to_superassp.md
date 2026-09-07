# Convert pladdrr data frame to superassp format

Converts data frames from pladdrr (long format with
time/frequency/value) to superassp's wide format expected for
AsspDataObj tracks.

## Usage

``` r
pladdrr_df_to_superassp(
  df,
  type = c("pitch", "formant", "intensity"),
  n_formants = 5
)
```

## Arguments

- df:

  Data frame from pladdrr R6 object's `as_data_frame()` method

- type:

  Character indicating data type: "pitch", "formant", "intensity"

- n_formants:

  Integer number of formants (default: 5, for formant data only)

## Value

Data frame in wide format suitable for AsspDataObj

## Details

Transformations by type:

**Pitch**: pladdrr returns (time, frequency) → superassp expects (time,
F0)

**Formant**: pladdrr returns (time, formant_number, frequency,
bandwidth) → superassp expects (time, fm1, fm2, ..., fm5, bw1, bw2, ...,
bw5)

**Intensity**: pladdrr returns (time, intensity) → superassp expects
same

## Examples

``` r
if (FALSE) { # \dontrun{
sound <- pladdrr::Sound("audio.wav")
pitch <- sound$to_pitch_cc()
pitch_df <- pitch$as_data_frame()
wide_df <- pladdrr_df_to_superassp(pitch_df, type = "pitch")
} # }
```
