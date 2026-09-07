# Convert AVAudio to av::read_audio_bin Format

Convert AVAudio object back to the format returned by
av::read_audio_bin(), which is an integer vector with channels and
sample_rate attributes.

## Usage

``` r
avaudio_to_av(audio)
```

## Arguments

- audio:

  AVAudio object

## Value

Integer vector with attributes (channels, sample_rate)

## Examples

``` r
if (FALSE) { # \dontrun{
audio <- read_avaudio("speech.wav")
audio_vec <- avaudio_to_av(audio)

# audio_vec is now compatible with av::read_audio_bin output
} # }
```
