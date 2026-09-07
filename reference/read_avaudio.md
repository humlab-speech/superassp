# Create AVAudio Object from File

Read an audio file and create an AVAudio object.

## Usage

``` r
read_avaudio(
  file_path,
  format = "wav",
  sample_rate = NULL,
  channels = NULL,
  start_time = NULL,
  end_time = NULL,
  ...
)
```

## Arguments

- file_path:

  Character; path to audio file

- format:

  Character; output format (default: "wav")

- sample_rate:

  Integer; target sample rate in Hz (default: NULL, keep original)

- channels:

  Integer; number of channels (default: NULL, keep original)

- start_time:

  Numeric; start time in seconds (default: NULL)

- end_time:

  Numeric; end time in seconds (default: NULL)

- ...:

  Additional arguments passed to
  [`prep_recode()`](https://humlab-speech.github.io/superassp/reference/prep_recode.md)

## Value

AVAudio object

## Examples

``` r
if (FALSE) { # \dontrun{
# Read entire file
audio <- read_avaudio("speech.wav")

# Read with time windowing
audio <- read_avaudio("speech.wav", start_time = 1.0, end_time = 3.0)

# Read with resampling
audio <- read_avaudio("speech.wav", sample_rate = 16000)
} # }
```
