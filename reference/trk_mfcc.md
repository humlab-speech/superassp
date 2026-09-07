# Extract Mel-Frequency Cepstral Coefficients (MFCCs) via SPTK

Computes HTK-style MFCCs using the SPTK C++ library. MFCCs are the
standard frame-level feature for speech recognition, speaker
identification, and general audio classification. Covers the full
filterbank-DCT pipeline with optional cepstral liftering.

## Usage

``` r
trk_mfcc(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 10.0 ms.

- windowSize:

  Numeric. Analysis window length in milliseconds. Default 25.0 ms.

- n_mfcc:

  Integer. Number of MFCC coefficients to extract (must be \< n_mels).
  Default 13.

- n_mels:

  Integer. Number of Mel filterbank channels. Default 40.

- fmin:

  Numeric. Lowest filterbank center frequency in Hz. Default 0.0 Hz.

- fmax:

  Numeric. Highest filterbank center frequency in Hz. `NULL` (default)
  uses the Nyquist frequency.

- lifter:

  Integer. Cepstral liftering exponent (HTK default is 22). Set to 0 to
  disable liftering. Default 22.

- floor:

  Numeric. Minimum energy floor for Mel filterbank outputs (prevents
  log(0)). Default 1.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"mfcc"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `mfcc_0` … `mfcc_{n_mfcc-1}`:

  REAL32, cepstral coefficients c0 through c{n_mfcc-1}, n_frames × 1
  each. Dimensionless.

Frame rate: `1000 / windowShift` Hz (default 100 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract 13 MFCCs (default)
trk_mfcc("recording.wav")

# Extract 20 MFCCs with custom parameters
trk_mfcc("speech.mp3", n_mfcc = 20, n_mels = 80)

# Return data without writing file
mfcc_data <- trk_mfcc("audio.wav", toFile = FALSE)

# Process with specific frequency range
trk_mfcc("recording.wav", fmin = 80, fmax = 8000)

# Process video file (extracts audio)
trk_mfcc("interview.mp4")
} # }
```
