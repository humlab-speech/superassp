# Compute the eGeMAPS openSMILE feature set

This function applies the extended version of the "Geneva Minimalistic
Acoustic Parameter Set (eGeMAPS) for Voice Research and Affective
Computing" (the *Extended Geneva Minimalistic Standard Parameter Set*,
eGeMAPS v02) (Eyben et al. 2015) to a portion of a recording.

## Usage

``` r
lst_eGeMAPS(listOfFiles, beginTime = 0, endTime = 0, explicitExt = "egm", toFile = FALSE, return_jstf = FALSE, outputDirectory = NULL)
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
  stored. Default "egm".

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- return_jstf:

  Logical. Return JsonTrackObj instead of list? Default FALSE. When both
  toFile and return_jstf are TRUE, writes the file AND returns the
  object.

- outputDirectory:

  Character. Output directory path. Default NULL (use input directory).

## Value

If `return_jstf=FALSE` and `toFile=FALSE` (default), a list of 88
acoustic values. If `toFile=TRUE`, invisibly returns the path(s) to the
written JSTF file(s). If `return_jstf=TRUE`, returns a JsonTrackObj. The
extendedacoustic parameter set contains the following compact set of 18
low-level descriptors (LLD), sorted by parameter groups:

Frequency related parameters:

- Pitch, logarithmic f0 on a semitone frequency scale,starting at 27.5
  Hz (semitone 0).

- Jitter, deviations in individual consecutive f0 period lengths.

- Formant 1, 2, and 3 frequency, centre frequency of first, second, and
  third formant

- Formant 1, bandwidth of first formant.Energy/Amplitude related
  parameters:

- Shimmer, difference of the peak amplitudes of consecutive f0 periods.

- Loudness, estimate of perceived signal intensity from an auditory
  spectrum.

- Harmonics-to-noise ratio (HNR), relation of energy in harmonic
  components to energy in noise-like components.

- Formant 2-3 bandwidth

Spectral (balance) parameters:

- Alpha Ratio, ratio of the summed energy from50-1000 Hz and 1-5 kHz

- Hammarberg Index, ratio of the strongest energy peak in the 0-2 kHz
  region to the strongest peak in the 2–5 kHz region.

- Spectral Slope 0-500 Hz and 500-1500 Hz, linear regression slope of
  the logarithmic power spectrum within the two given bands.

- Formant 1, 2, and 3 relative energy, as well as the ratio of the
  energy of the spectral harmonic peak at the first, second, third
  formant's centre frequency to the energy of the spectral peak atF0.

- Harmonic difference H1-H2, ratio of energy of the first f0 harmonic
  (H1) to the energy of the second f0 harmonic (H2).

- Harmonic difference H1-A3, ratio of energy of the first f0harmonic
  (H1) to the energy of the highest harmonic in the third formant range
  (A3).

- MFCC 1-4 Mel-Frequency Cepstral Coefficients 1-4.

- Spectral flux difference of the spectra of two consecutive frames.

which are analysed in terms of mean and coefficient of variation, as
well as 20th, median (50th), and 80th percentile (pitch and loudness),
the arithmetic mean of the Alpha Ratio, the Hammarberg Index, and the
spectral slopes from 0-500 Hz and 500-1500 Hz over all unvoiced
segments, and the equivalent sound level.

Temporal features:

- the rate of loudness peaks, i.e., the number of loudness peaks per
  second,

- the mean length and the standard deviation of continuously voiced
  regions(f0\>0),

- the mean length and the standard deviation of unvoiced regions (f0 ==
  0; approximating pauses),

- the number of continuous voiced regions per second(pseudo syllable
  rate).

Please consult the (Eyben et al. 2015) for a description of the
features.

## Details

The GeMAPS feature set consists of of 88 static acoustic features
resulting from the computation of various functionals over low-level
descriptor features, and is applied by this function using the
openSMILE. (Eyben et al. 2010; Jaimes et al. 2013) acoustic feature
extraction library.

## References

Eyben F, Scherer KR, Schuller BW, Sundberg J, Andre E, Busso C,
Devillers LY, Epps J, Laukka P, Narayanan SS, Truong KP (2015). “The
Geneva Minimalistic Acoustic Parameter Set (GeMAPS) for Voice Research
and Affective Computing.” *IEEE Transactions on Affective Computing*,
**7**(2), 190–202. ISSN 1949-3045.
[doi:10.1109/taffc.2015.2457417](https://doi.org/10.1109/taffc.2015.2457417)
.  
  
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
[doi:10.1145/2502081.2502224](https://doi.org/10.1145/2502081.2502224) .

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute eGeMAPS features
egemaps <- lst_eGeMAPS("audio.wav")

# With time windowing
egemaps <- lst_eGeMAPS("audio.wav", beginTime = 1.0, endTime = 3.0)

# Write results to JSTF file
lst_eGeMAPS("audio.wav", toFile = TRUE)  # Creates audio.egm

# Read back and convert to data.frame
track <- read_track("audio.egm")
df <- as.data.frame(track)
head(df)  # Shows begin_time, end_time, and all 88 eGeMAPS features
} # }
```
