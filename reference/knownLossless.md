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
superassp:::knownLossless()
#>  [1] "wav"  "flac" "aiff" "wv"   "ape"  "tta"  "caf"  "au"   "w64"  "dsf" 
#> [11] "dff"  "kay"  "nist" "nsp" 

# Check if a file extension is lossless
"flac" %in% superassp:::knownLossless()  # TRUE
#> [1] TRUE
"mp3" %in% superassp:::knownLossless()   # FALSE
#> [1] FALSE
```
