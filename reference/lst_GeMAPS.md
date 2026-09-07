# Compute the GeMAPS openSMILE feature set (C++ Implementation)

This function applies the "The Geneva Minimalistic Acoustic Parameter
Set (GeMAPS) for Voice Research and Affective Computing" (the *Geneva
Minimalistic Standard Parameter Set*, GeMAPS v0.1b) (Eyben et al. 2015)
to a portion of a recording using the native OpenSMILE C++ library for
maximum performance.

## Usage

``` r
lst_GeMAPS(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  The full path to the sound file.

- beginTime:

  The starting time of the section of the sound files that should be
  analysed (in seconds).

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- explicitExt:

  The file extension of the slice file where the results should be
  stored. Default "gem".

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

If `return_jstf=FALSE` and `toFile=FALSE` (default), a named list of 62
acoustic values. If `toFile=TRUE`, invisibly returns the path(s) to the
written JSTF file(s). If `return_jstf=TRUE`, returns a JsonTrackObj.

## Details

The GeMAPS feature set consists of 62 static acoustic features (pitch,
jitter, shimmer, formants, HNR, spectral balance, and loudness) computed
via the openSMILE C++ library directly (3-5x faster than the Python
implementation). See the openSMILE GeMAPS paper ((Eyben et al. 2015) )
for the full definition of each descriptor and functional.

## References

Eyben F, Scherer KR, Schuller BW, Sundberg J, Andre E, Busso C,
Devillers LY, Epps J, Laukka P, Narayanan SS, Truong KP (2015). “The
Geneva Minimalistic Acoustic Parameter Set (GeMAPS) for Voice Research
and Affective Computing.” *IEEE Transactions on Affective Computing*,
**7**(2), 190–202. ISSN 1949-3045.
[doi:10.1109/taffc.2015.2457417](https://doi.org/10.1109/taffc.2015.2457417)
.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using C++ implementation (default, fastest)
gemaps <- lst_GeMAPS("audio.wav")

# With time windowing
gemaps <- lst_GeMAPS("audio.wav", beginTime = 1.0, endTime = 3.0)

# Write results to JSTF file
lst_GeMAPS("audio.wav", toFile = TRUE)  # Creates audio.gem

# Read back and convert to data.frame
track <- read_track("audio.gem")
df <- as.data.frame(track)
head(df)  # Shows begin_time, end_time, and all 62 GeMAPS features
} # }
```
