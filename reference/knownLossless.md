# Procedure reporting lossless file formats

This procedure reports file extensions that are known to contain
losslessly encoded sound data. These formats preserve the original audio
signal without quality loss, which is essential for accurate DSP
(Digital Signal Processing) analysis.

## Usage

``` r
knownLossless()
```

## Value

Character vector of file extensions (without dots)

## Details

The list includes formats supported by the av package and native
ASSP/wrassp formats:

**Container formats (lossless codecs):**

- wav: WAV / WAVE (Waveform Audio) - most common

- flac: Free Lossless Audio Codec

- aiff: Audio Interchange File Format (Apple)

- wv: WavPack

- ape: Monkey's Audio

- tta: True Audio

- caf: Apple Core Audio Format

- au: Sun/NeXT Audio

- w64: Sony Wave64 (64-bit WAV variant)

**High-resolution audio:**

- dsf: DSD Stream File (Direct Stream Digital)

- dff: DSD Interchange File Format

**Professional/scientific formats:**

- kay: Kay Elemetrics CSL files

- nist: NIST SPHERE

- nsp: NSP (used in speech research)

## Author

Fredrik Nylén

## Examples

``` r
# Get list of lossless formats
knownLossless()
#> Error in knownLossless(): could not find function "knownLossless"

# Check if a file extension is lossless
"flac" %in% knownLossless()  # TRUE
#> Error in knownLossless(): could not find function "knownLossless"
"mp3" %in% knownLossless()   # FALSE
#> Error in knownLossless(): could not find function "knownLossless"
```
