# Compute the ComParE 2016 openSMILE feature set

This function applies the "The INTERSPEECH 2016 Computational
Paralinguistics Challenge: Deception, Sincerity & Native Language"
ComParE (Weninger et al. 2013) to a portion of a recording.

## Usage

``` r
lst_ComParE_2016(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  The full path to the sound file.

- beginTime:

  The starting time of the section of the sound files that should be
  analysed.

- endTime:

  The end time of the sound files that should be analysed.

- explicitExt:

  The file extension of the slice file where the results should be
  stored. Default "cmp".

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- return_jstf:

  Logical. Return JsonTrackObj instead of list? Default FALSE. When both
  toFile and return_jstf are TRUE, writes the file AND returns the
  object.

- outputDirectory:

  Character. Output directory path. Default NULL (use input directory).

## Value

If `return_jstf=FALSE` and `toFile=FALSE` (default), a list of 6,373
acoustic values. If `toFile=TRUE`, invisibly returns the path(s) to the
written JSTF file(s). If `return_jstf=TRUE`, returns a JsonTrackObj.
Please consult the (Weninger et al. 2013) for a description of the
features.

## Details

The ComParE feature set consists of of 6 373 static acoustic features
resulting from the computation of various functionals over low-level
descriptor features, and is applied by this function using the openSMILE
(Eyben et al. 2010; Jaimes et al. 2013) acoustic feature extraction
library.

## References

Eyben F, Wöllmer M, Schuller B (2010). *Opensmile: the munich versatile
and fast open-source audio feature extractor*, the international
conference. ACM. ISBN 978-1-60558-933-6.
[doi:10.1145/1873951.1874246](https://doi.org/10.1145/1873951.1874246) .
<http://dl.acm.org/citation.cfm?id=1874246>.  
  
Jaimes A(, Sebe N, Boujemaa N, Gatica-Perez D, Shamma DA, Worring M,
Zimmermann R, Eyben F, Weninger F, Gross F, Schuller B (2013). “Recent
developments in openSMILE, the munich open-source multimedia feature
extractor.” *Proceedings of the 21st ACM international conference on
Multimedia*, 835–838.
[doi:10.1145/2502081.2502224](https://doi.org/10.1145/2502081.2502224)
.  
  
Weninger F, Eyben F, Schuller BW, Mortillaro M, Scherer KR (2013). “On
the Acoustics of Emotion in Audio: What Speech, Music, and Sound have in
Common.” *Frontiers in Psychology*, **4**, 292. ISSN 1664-1078.
[doi:10.3389/fpsyg.2013.00292](https://doi.org/10.3389/fpsyg.2013.00292)
.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using C++ implementation (default, fastest)
compare <- lst_ComParE_2016("audio.wav")

# With time windowing
compare <- lst_ComParE_2016("audio.wav", beginTime = 1.0, endTime = 3.0)

# Write results to JSTF file
lst_ComParE_2016("audio.wav", toFile = TRUE)  # Creates audio.cmp

# Read back and convert to data.frame
track <- read_track("audio.cmp")
df <- as.data.frame(track)
head(df)  # Shows begin_time, end_time, and all 6,373 ComParE features
} # }
```
