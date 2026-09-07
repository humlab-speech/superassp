# Perform RMS analysis on AsspDataObj in memory

This function performs RMS analysis on an AsspDataObj without writing
intermediate files. It's useful when processing audio data read from
video files or other non-WAV formats using the av package.

## Usage

``` r
rmsana_memory(audio_obj, ...)
```

## Arguments

- audio_obj:

  AsspDataObj containing audio data

- ...:

  Additional parameters passed to rmsana

## Value

AsspDataObj with RMS analysis results

## Examples

``` r
if (FALSE) { # \dontrun{
audio_obj <- av_to_asspDataObj("myfile.mp4")
rms_result <- rmsana_memory(audio_obj, windowShift = 5)
} # }
```
