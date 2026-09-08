# Re-encode Media File with Custom Parameters

Re-encodes any media file (audio/video) supported by the av package into
16-bit PCM WAV, with optional resampling, channel remixing, and time
windowing. Returns the audio data in the same format as
[`av::read_audio_bin`](https://docs.ropensci.org/av//reference/read_audio.html).

## Usage

``` r
prep_recode(
  listOfFiles,
  codec,
  sample_rate = NULL,
  bit_rate = NULL,
  start_time = NULL,
  end_time = NULL,
  channels = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- listOfFiles:

  Character vector of file paths to media files

- codec:

  Either `"none"` (read the file as-is, no re-encoding) or `"pcm_s16le"`
  (re-encode to 16-bit PCM WAV). Required. These are the only two
  re-encoding needs superassp has internally; for anything else call
  [`av::av_audio_convert()`](https://docs.ropensci.org/av//reference/encoding.html)
  directly.

- sample_rate:

  Target sample rate in Hz (default: NULL keeps original)

- bit_rate:

  Ignored for `"pcm_s16le"` (lossless); kept for interface symmetry with
  [`av::av_audio_convert()`](https://docs.ropensci.org/av//reference/encoding.html).

- start_time:

  Start time in seconds (default: NULL = start of file)

- end_time:

  End time in seconds (default: NULL = end of file)

- channels:

  Number of output channels: 1 (mono), 2 (stereo), or NULL (keep
  original)

- verbose:

  Logical; show progress messages (default: TRUE)

- ...:

  Additional arguments passed to
  [`av::av_audio_convert`](https://docs.ropensci.org/av//reference/encoding.html)

## Value

For single file: Integer vector with audio samples in s32le format
(32-bit signed integers), with attributes:

- `channels`: Number of audio channels (integer)

- `sample_rate`: Sample rate in Hz (integer)

For multiple files: List of integer vectors, one per file

This matches the format returned by
[`av::read_audio_bin()`](https://docs.ropensci.org/av//reference/read_audio.html).

## Details

Re-encoding goes through a temporary WAV file
([`av::av_audio_convert()`](https://docs.ropensci.org/av//reference/encoding.html)
followed by
[`av::read_audio_bin()`](https://docs.ropensci.org/av//reference/read_audio.html));
the temp file is removed on exit. It's useful for:

- Converting sample rates for analysis

- Extracting audio from video files

- Time-windowing large files

- Remixing channel count across a corpus

## References

Ooms J (2024). *av: Working with Audio and Video in R*. rOpenSci. R
package, <https://docs.ropensci.org/av/>.

FFmpeg Developers (2024). *FFmpeg Codecs Documentation*. FFmpeg project.
<https://ffmpeg.org/ffmpeg-codecs.html>.

## See also

[`av_audio_convert`](https://docs.ropensci.org/av//reference/encoding.html),
[`read_audio_bin`](https://docs.ropensci.org/av//reference/read_audio.html),
[`av_to_asspDataObj`](https://humlab-speech.github.io/superassp/reference/av_to_asspDataObj.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read as-is
audio <- prep_recode("speech.wav", codec = "none")

# Extract segment from 1-3 seconds
audio_segment <- prep_recode("long.wav",
                              codec = "pcm_s16le",
                              start_time = 1.0,
                              end_time = 3.0)

# Downsample to 16 kHz
audio_16k <- prep_recode("high_res.wav",
                         codec = "pcm_s16le",
                         sample_rate = 16000)

# Convert to mono
audio_mono <- prep_recode("stereo.wav",
                          codec = "pcm_s16le",
                          channels = 1)

# Batch processing
files <- c("file1.mp4", "file2.wav", "file3.flac")
audio_list <- prep_recode(files,
                          codec = "pcm_s16le",
                          sample_rate = 44100,
                          channels = 1)

# Access audio data (same as av::read_audio_bin)
audio <- prep_recode("test.wav", codec = "pcm_s16le")
cat("Channels:", attr(audio, "channels"), "\n")
cat("Sample rate:", attr(audio, "sample_rate"), "\n")
cat("Duration:", length(audio) / attr(audio, "channels") / attr(audio, "sample_rate"), "s\n")
} # }
```
