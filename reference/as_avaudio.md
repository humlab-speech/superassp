# Convert to AVAudio Object

Convert audio data (from av::read_audio_bin or prep_recode) to AVAudio
object.

## Usage

``` r
as_avaudio(x, file_path = NA_character_)
```

## Arguments

- x:

  Integer vector with audio samples (must have channels and sample_rate
  attributes)

- file_path:

  Character; optional source file path

## Value

AVAudio object

## Examples

``` r
if (FALSE) { # \dontrun{
# From av::read_audio_bin
audio_data <- av::read_audio_bin("speech.wav")
audio <- as_avaudio(audio_data)

# From prep_recode
audio_data <- prep_recode("speech.wav", format = "wav")
audio <- as_avaudio(audio_data, file_path = "speech.wav")
} # }
```
