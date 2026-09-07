# Re-encode Media File with Custom Parameters

Re-encodes any media file (audio/video) supported by the av package into
a specified format with custom codec, sample rate, bit rate, and
optional time windowing. Returns the audio data in the same format as
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

  Output codec (e.g., "pcm_s16le", "mp3", "flac", "vorbis"). Required.
  See
  [`av::av_encoders()`](https://docs.ropensci.org/av//reference/formats.html)
  for available codecs.

- sample_rate:

  Target sample rate in Hz (default: NULL keeps original)

- bit_rate:

  Target bit rate for lossy codecs (default: NULL uses codec default).
  Specify as integer (bits/second), e.g., 128000, 192000, 320000

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

This function performs in-memory transcoding using
[`av::av_audio_transcode()`](https://docs.ropensci.org/av//reference/encoding.html),
avoiding intermediate files on disk. It's useful for:

- Converting sample rates for analysis

- Extracting audio from video files

- Time-windowing large files

- Normalizing formats across a corpus

- Testing codec-specific effects

**Supported Formats:**

The av package supports a wide range of formats through FFmpeg:

- **Lossless:** wav, flac, alac, ape, wv

- **Lossy:** mp3, ogg, aac, opus, wma

- **Video:** mp4, mkv, avi, mov, webm (extracts audio)

**Common Codec Examples:**

- **WAV:** "pcm_s16le" (16-bit), "pcm_s24le" (24-bit), "pcm_f32le"
  (32-bit float)

- **MP3:** "libmp3lame"

- **FLAC:** "flac"

- **OGG:** "libvorbis"

- **AAC:** "aac"

- **OPUS:** "libopus"

**Processing Strategy:**

1.  If no re-encoding needed (no codec/sample_rate/channels change, no
    windowing):

    - Returns
      [`av::read_audio_bin()`](https://docs.ropensci.org/av//reference/read_audio.html)
      result directly

2.  If re-encoding or windowing needed:

    - Uses
      [`av::av_audio_transcode()`](https://docs.ropensci.org/av//reference/encoding.html)
      for in-memory transcoding

    - Returns audio data directly (no temporary files)

**Performance:**

- Pure in-memory operation (no temporary files)

- Fast conversion for compatible codecs

- Time windowing reduces memory usage

## References

Ooms J (2024). *av: Working with Audio and Video in R*. rOpenSci. R
package, <https://docs.ropensci.org/av/>.

FFmpeg Developers (2024). *FFmpeg Codecs Documentation*. FFmpeg project.
<https://ffmpeg.org/ffmpeg-codecs.html>.

## See also

[`av_audio_transcode`](https://docs.ropensci.org/av//reference/encoding.html),
[`read_audio_bin`](https://docs.ropensci.org/av//reference/read_audio.html),
[`av_to_asspDataObj`](https://humlab-speech.github.io/superassp/reference/av_to_asspDataObj.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage - convert to WAV PCM
audio <- prep_recode("video.mp4", codec = "pcm_s16le")

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

# Convert to MP3 with specific bit rate
audio_mp3 <- prep_recode("speech.wav",
                         codec = "mp3",
                         bit_rate = 192000)

# Convert to FLAC (lossless compression)
audio_flac <- prep_recode("recording.wav",
                          codec = "flac")

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
