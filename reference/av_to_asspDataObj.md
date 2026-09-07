# Convert audio file to AsspDataObj

Reads audio into an AsspDataObj for in-memory DSP processing. Native
formats (wav, au, kay, nist, nsp) are read via the C-level
[`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md)
for maximum speed. All other formats (mp3, mp4, aac, flac, ogg, …) use
[`av::read_audio_bin()`](https://docs.ropensci.org/av//reference/read_audio.html)
as the fallback. Resampling (via `target_sample_rate`) always uses the
av path.

## Usage

``` r
av_to_asspDataObj(
  file_path,
  start_time = 0,
  end_time = NULL,
  target_sample_rate = NULL
)
```

## Arguments

- file_path:

  Path to the audio/video file

- start_time:

  Start time in seconds (default 0)

- end_time:

  End time in seconds (default NULL for end of file)

- target_sample_rate:

  Target sample rate (default NULL to keep original).

## Value

An AsspDataObj containing the audio data

## Examples

``` r
if (FALSE) { # \dontrun{
audio_obj <- av_to_asspDataObj("myfile.wav")
audio_obj <- av_to_asspDataObj("myfile.mp3")
audio_obj <- av_to_asspDataObj("myfile.wav", start_time = 1.0, end_time = 3.0)
} # }
```
