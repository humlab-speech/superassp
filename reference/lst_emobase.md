# Compute the emobase openSMILE feature set

This function applies the emobase openSMILE (Eyben et al. 2010)
configuration to compute 988 acoustic features reasoned to be part of
perception of emotion.

## Usage

``` r
lst_emobase(listOfFiles, beginTime = 0, endTime = 0, explicitExt = "emb", verbose = FALSE, toFile = FALSE, return_jstf = FALSE, outputDirectory = NULL)
```

## Arguments

- listOfFiles:

  The full path to the sound file.

- beginTime:

  The starting time of the section of the sound files that should be
  analysed.

- endTime:

  The end time of the section of the sound files that should be
  analysed.

- explicitExt:

  The file extension of the slice file where the results should be
  stored. Default "emb".

- verbose:

  Logical. Print processing information (default: FALSE).

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- return_jstf:

  Logical. Return JsonTrackObj instead of list? Default FALSE. When both
  toFile and return_jstf are TRUE, writes the file AND returns the
  object.

- outputDirectory:

  Character. Output directory path. Default NULL (use input directory).

## Value

If `return_jstf=FALSE` and `toFile=FALSE` (default), a list of 988
acoustic values. If `toFile=TRUE`, invisibly returns the path(s) to the
written JSTF file(s). If `return_jstf=TRUE`, returns a JsonTrackObj.

## References

Eyben F, Wöllmer M, Schuller B (2010). *Opensmile: the munich versatile
and fast open-source audio feature extractor*, the international
conference. ACM. ISBN 978-1-60558-933-6.
[doi:10.1145/1873951.1874246](https://doi.org/10.1145/1873951.1874246) .
<http://dl.acm.org/citation.cfm?id=1874246>.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using C++ implementation (default, fastest)
emobase <- lst_emobase("audio.wav")

# With time windowing
emobase <- lst_emobase("audio.wav", beginTime = 1.0, endTime = 3.0)

# Write results to JSTF file
lst_emobase("audio.wav", toFile = TRUE)  # Creates audio.emb

# Read back and convert to data.frame
track <- read_track("audio.emb")
df <- as.data.frame(track)
head(df)  # Shows begin_time, end_time, and all 988 emobase features
} # }
```
