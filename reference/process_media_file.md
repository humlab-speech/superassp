# Process audio from any media file format

This function provides a convenient interface to process audio from any
media file format supported by av (including video files), perform
superassp analysis, and return results without creating intermediate
files.

## Usage

``` r
process_media_file(
  file_path,
  analysis_function = "trk_rms",
  start_time = 0,
  end_time = NULL,
  target_sample_rate = NULL,
  ...
)
```

## Arguments

- file_path:

  Path to the media file

- analysis_function:

  Name of analysis function (e.g., "trk_rms", "trk_formant_forest")

- start_time:

  Start time in seconds (default 0)

- end_time:

  End time in seconds (default NULL for end)

- target_sample_rate:

  Target sample rate (default NULL for original)

- ...:

  Additional parameters for the analysis function

## Value

Result from the analysis function

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze RMS from a video file
rms <- process_media_file("video.mp4", "trk_rms", windowShift = 5)

# Analyze formants from audio segment
formants <- process_media_file("audio.m4a", "trk_formant_forest",
                                start_time = 10, end_time = 20)
} # }
```
