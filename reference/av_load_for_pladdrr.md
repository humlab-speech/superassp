# Load audio file as pladdrr Sound object

Simple wrapper to load audio files directly with pladdrr. Handles time
windowing via pladdrr's native `extract_part()` method.

## Usage

``` r
av_load_for_pladdrr(
  file_path,
  start_time = 0,
  end_time = 0,
  window_type = "Gaussian1",
  relative_width = 1
)
```

## Arguments

- file_path:

  Character path to audio file (WAV, AIFF, FLAC, MP3, etc.)

- start_time:

  Numeric start time in seconds (default: 0 = beginning)

- end_time:

  Numeric end time in seconds (default: 0 = end of file)

- window_type:

  Character window shape for part extraction (default: "Gaussian1")

- relative_width:

  Numeric relative width of window (default: 1.0)

## Value

A pladdrr Sound R6 object

## Details

pladdrr reads audio files directly via its C++ backend (no av needed).
Supported formats: WAV, AIFF, FLAC, MP3, NIST (via native readers), or
any format via av package fallback.

If time windowing is requested (start_time \> 0 or end_time \> 0), the
function uses pladdrr's `extract_part()` method with appropriate
windowing.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load entire file
sound <- av_load_for_pladdrr("speech.wav")

# Load with time windowing (1.0 to 3.0 seconds)
sound <- av_load_for_pladdrr("speech.wav", start_time = 1.0, end_time = 3.0)

# Use with pladdrr functions
pitch <- sound$to_pitch_cc(time_step = 0.01, pitch_floor = 75, pitch_ceiling = 600)
} # }
```
