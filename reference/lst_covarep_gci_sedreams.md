# SEDREAMS Glottal Closure Instant Detection

Detects glottal closure instants (GCIs) using SEDREAMS algorithm.
Returns event times (GCI instants), not a regular frame grid.

## Usage

``` r
lst_covarep_gci_sedreams(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Vector of file paths (WAV, MP3, MP4, etc.) to analyze

- beginTime:

  Start time in seconds (0 for beginning of file)

- endTime:

  End time in seconds (0 for end of file)

- f0mean:

  Estimated mean F0 in Hz. If NULL, auto-estimated from signal.

- polarity:

  Signal polarity (1 or -1). If NULL, auto-detected.

- verbose:

  Show progress messages (default: TRUE)

## Value

Data frame with columns:

- `file`: Input file path

- `n_gcis`: Number of detected GCIs

- `gci_times`: List column of numeric vectors (GCI times in seconds)

## Details

**SEDREAMS Algorithm** (Ney and Kneser 2002) :

1.  Compute LPC residual (25ms frames, 5ms shift, order ≈ fs/1000 + 2)

2.  Bandpass filter signal around estimated F0 (mean-based signal)

3.  Find maxima/minima pairs in mean-based signal

4.  Locate GCI positions in LP residual peaks within windows

**Typical output**:

- Voiced speech: 100-200 GCIs per second (F0-dependent)

- Unvoiced/silence: 0 GCIs (no glottal closures)

**Use cases**:

- Foundation for GCI-based voice quality (NAQ, QOQ, H1H2 via
  trk_covarep_vq_gci)

- Voice pathology assessment (irregular GCI spacing = vocal pathology)

- Glottal source analysis (GCI-anchored inverse filtering)

- Speech analysis (pitch period estimation, voicing detection)

**Downstream workflow**:

1.  `lst_covarep_gci_sedreams()` — detect GCIs

2.  [`trk_covarep_vq_gci()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_vq_gci.md)
    — compute voice quality per GCI

3.  [`lst_covarep_vq()`](https://humlab-speech.github.io/superassp/reference/lst_covarep_vq.md)
    — summarize to scalars

## References

Ney H, Kneser R (2002). “Speech recognition using continuous-space
embeddings.” *IEEE Signal Processing Magazine*, **19**(1), 33–42. Signal
processing foundations for GCI detection and SEDREAMS algorithm.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single file
gcis <- lst_covarep_gci_sedreams("speech.wav", f0mean = 100)

# Batch process
files <- c("file1.wav", "file2.wav")
results <- lst_covarep_gci_sedreams(files, f0mean = 110)

# View results
results$gci_times[[1]]  # GCI times for first file (in seconds)

# Chain with voice quality analysis
vq <- trk_covarep_vq_gci("speech.wav", gci_times = gcis$gci_times[[1]])
} # }
```
