# Vocal Tremor Analysis Using pladdrr

Analyzes vocal tremor from sustained vowel recordings using pladdrr's
Praat bindings. Extracts 18 measures of frequency and amplitude tremor
based on Brückl (2012) autocorrelation algorithm.

## Usage

``` r
lst_voice_tremor(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  analysisTimeStep = 0.015,
  minPitch = 60,
  maxPitch = 350,
  silenceThreshold = 0.03,
  voicingThreshold = 0.3,
  octaveCost = 0.01,
  octaveJumpCost = 0.35,
  voicedUnvoicedCost = 0.14,
  minTremorFreq = 1.5,
  maxTremorFreq = 15,
  tremorMagThresh = 0.01,
  tremorCyclicalThresh = 0.15,
  freqTremorOctaveCost = 0.01,
  ampTremorOctaveCost = 0.01,
  nanAsZero = FALSE,
  toFile = FALSE,
  return_jstf = FALSE,
  explicitExt = "pvt",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector with path(s) to audio file(s)

- beginTime:

  Numeric. Start time in seconds (default 0)

- endTime:

  Numeric. End time in seconds (0 = end of file)

- analysisTimeStep:

  Numeric. Time step for analysis in seconds (default 0.015)

- minPitch:

  Numeric. Minimum pitch for extraction in Hz (default 60)

- maxPitch:

  Numeric. Maximum pitch for extraction in Hz (default 350)

- silenceThreshold:

  Numeric. Threshold for silence detection (default 0.03)

- voicingThreshold:

  Numeric. Threshold for voicing detection (default 0.3)

- octaveCost:

  Numeric. Cost for octave jumps in pitch tracking (default 0.01)

- octaveJumpCost:

  Numeric. Cost for large octave jumps (default 0.35)

- voicedUnvoicedCost:

  Numeric. Cost for voiced/unvoiced transitions (default 0.14)

- minTremorFreq:

  Numeric. Minimum tremor frequency in Hz (default 1.5)

- maxTremorFreq:

  Numeric. Maximum tremor frequency in Hz (default 15)

- tremorMagThresh:

  Numeric. Threshold for contour magnitude (default 0.01)

- tremorCyclicalThresh:

  Numeric. Threshold for cyclicality (default 0.15)

- freqTremorOctaveCost:

  Numeric. Octave cost for frequency tremor (default 0.01)

- ampTremorOctaveCost:

  Numeric. Octave cost for amplitude tremor (default 0.01)

- nanAsZero:

  Logical. Convert undefined measurements to zeros (default FALSE)

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- explicitExt:

  Character. File extension for output. Default "pvt".

- outputDirectory:

  Character. Output directory path. Default NULL (use input directory).

- verbose:

  Logical. Print progress messages (default TRUE)

- return_jstf:

  Logical. Return JsonTrackObj instead of data.frame? Default FALSE.
  When both toFile and return_jstf are TRUE, the file is written AND the
  object returned.

## Value

If `toFile=FALSE` (default), a data.frame (or list of data.frames for
multiple files) with 18 tremor measurements. If `toFile=TRUE`, invisibly
returns the path(s) to the written JSTF file(s).

Each result contains 18 columns:

- FCoM:

  Frequency contour magnitude

- FTrC:

  Frequency tremor cyclicality (0-1)

- FMoN:

  Number of frequency modulation candidates

- FTrF:

  Frequency tremor frequency (Hz)

- FTrI:

  Frequency tremor intensity index (percent)

- FTrP:

  Frequency tremor power index

- FTrCIP:

  Frequency tremor cyclicality-intensity product

- FTrPS:

  Frequency tremor product sum

- FCoHNR:

  Frequency contour HNR (dB)

- ACoM:

  Amplitude contour magnitude

- ATrC:

  Amplitude tremor cyclicality (0-1)

- AMoN:

  Number of amplitude modulation candidates

- ATrF:

  Amplitude tremor frequency (Hz)

- ATrI:

  Amplitude tremor intensity index (percent)

- ATrP:

  Amplitude tremor power index

- ATrCIP:

  Amplitude tremor cyclicality-intensity product

- ATrPS:

  Amplitude tremor product sum

- ACoHNR:

  Amplitude contour HNR (dB)

## Details

This function processes sustained phonations to detect tremor
characteristics in both pitch (frequency) and intensity (amplitude)
contours. It applies Gaussian1 windowing and uses autocorrelation-based
analysis to identify tremor frequency, intensity, and cyclicality.

## References

(Brückl 2012)

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze sustained vowel
result <- lst_voice_tremor("sustained_vowel.wav")
print(result$FTrF)  # Frequency tremor frequency
print(result$FTrI)  # Frequency tremor intensity

# Write to JSTF file
lst_voice_tremor("sustained_vowel.wav", toFile = TRUE)
track <- read_track("sustained_vowel.pvt")
df <- as.data.frame(track)
} # }
```
