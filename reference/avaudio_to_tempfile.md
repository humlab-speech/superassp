# Convert AVAudio to Temporary WAV File

Write AVAudio object to a temporary WAV file and return the path. Useful
for passing to external DSP functions that require file paths.

## Usage

``` r
avaudio_to_tempfile(audio, verbose = FALSE)
```

## Arguments

- audio:

  AVAudio object

- verbose:

  Logical; show messages (default: FALSE)

## Value

Character; path to temporary WAV file

## Details

The temporary file is created with
[`tempfile()`](https://rdrr.io/r/base/tempfile.html) and will be
automatically deleted when the R session ends, unless explicitly deleted
earlier with [`unlink()`](https://rdrr.io/r/base/unlink.html).

## Examples

``` r
if (FALSE) { # \dontrun{
audio <- read_avaudio("speech.wav")
temp_path <- avaudio_to_tempfile(audio)

# Use temp file with external tools
# ... processing ...

# Clean up (optional - happens automatically at session end)
unlink(temp_path)
} # }
```
