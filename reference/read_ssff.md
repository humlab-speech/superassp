# Read an SSFF or audio file into an AsspDataObj

A user-facing wrapper around the internal ASSP C-level reader. Interface
is identical to the legacy `read.AsspDataObj`.

## Usage

``` r
read_ssff(fname, begin = 0, end = 0, samples = FALSE)
```

## Arguments

- fname:

  Path to an SSFF or native ASSP audio file (WAV, AU, NIST, etc.).

- begin:

  Start of region to read (seconds, or samples if `samples=TRUE`).
  Default 0 = file start.

- end:

  End of region to read (seconds, or samples if `samples=TRUE`). Default
  0 = file end.

- samples:

  Logical. If `TRUE`, `begin`/`end` are in samples; otherwise in
  seconds.

## Value

An `AsspDataObj`. For audio files, contains an `audio` track (n_samples
x n_channels). For SSFF tracks, contains one matrix per stored track
(e.g. `F0`, `fm`, `bw`, `rms`) at the analysis frame rate. Standard
attributes include `sampleRate`, `startTime`, `startRecord`,
`endRecord`, `trackFormats` and `filePath`.

## See also

[`read_audio`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
for universal format support including MP3/MP4.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read an audio file from the bundled samples
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
au  <- read_ssff(wav)
names(au)               # "audio"
attr(au, "sampleRate")  # native sample rate

# Read an SSFF parameter track produced earlier by a trk_* function
f0_path <- tempfile(fileext = ".f0")
trk_pitch_rapt(wav, toFile = TRUE, outputDirectory = dirname(f0_path),
               explicitExt = "f0")
f0_obj <- read_ssff(file.path(dirname(f0_path),
                              paste0(tools::file_path_sans_ext(basename(wav)), ".f0")))
names(f0_obj)
} # }
```
