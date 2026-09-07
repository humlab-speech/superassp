# Read an audio file into an AsspDataObj

Unified audio reader. Tries the ASSP C-level reader first (fastest,
supports WAV, AU, NIST, SSFF and other native formats). Falls back to
the `av` package for any format av/FFmpeg can decode (MP3, MP4, FLAC,
OGG, etc.).

## Usage

``` r
read_audio(fname, begin = 0, end = 0, samples = FALSE)
```

## Arguments

- fname:

  Path to the audio or SSFF file.

- begin:

  Start of region to read. Default 0 = file start. In seconds (default)
  or samples (when `samples = TRUE`).

- end:

  End of region to read. Default 0 = file end. In seconds or samples.

- samples:

  Logical. If `TRUE`, `begin` and `end` are interpreted as sample
  indices.

## Value

An `AsspDataObj` containing one `audio` track (`n_samples` x
`n_channels`, INT16) plus the standard sample-rate metadata. Attributes
include `sampleRate`, `startTime`, `startRecord`, `endRecord`,
`trackFormats` and `filePath`. Frame rate equals the audio sample rate;
one record per audio sample.

## Details

For variable-rate encoded formats (e.g. MP3) with `samples = TRUE`,
sample positions are approximated via the nominal sample rate reported
by the container.

## See also

[`read_ssff`](https://humlab-speech.github.io/superassp/reference/read_ssff.md)
for SSFF-only files;
[`read_jstf`](https://humlab-speech.github.io/superassp/reference/read_jstf.md)
for JSTF files.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read entire WAV file
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
au  <- read_audio(wav)

# Read a 1-second window starting at 0.5 s
clip <- read_audio(wav, begin = 0.5, end = 1.5)

# Read by sample index
first_1k <- read_audio(wav, begin = 1, end = 1000, samples = TRUE)

# MP3/MP4 fall back to FFmpeg automatically
# mp3 <- read_audio("recording.mp3")
} # }
```
