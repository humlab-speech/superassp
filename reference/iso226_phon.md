# ISO 226:2023 Phon (Loudness Level) Conversions

Functions for converting between sound pressure level (dB) and loudness
level (phon) according to ISO 226:2023 "Acoustics - Normal
equal-loudness-level contours".

## Details

The phon scale represents loudness level - the perceived loudness of a
sound compared to a 1000 Hz reference tone. Unlike other psychoacoustic
scales (Bark, ERB, Mel), phon is a **bivariate measure** requiring both
frequency and sound pressure level.

**Key concept**: A sound at 40 phon sounds as loud as a 1000 Hz tone at
40 dB SPL, but different frequencies require different SPLs to achieve
the same loudness.

**Valid ranges** (per ISO 226:2023):

- Frequency: 20 Hz to 12,500 Hz

- Loudness level: 20 phon to 90 phon (20-4000 Hz), 20-80 phon
  (5000-12500 Hz)

- Below 20 phon: Informative only (near hearing threshold)

- Above 90/80 phon: Limited experimental data

**Implementation notes**:

- Based on ISO 226:2023 Formulas (1) and (2)

- Parameters from ISO 226:2023 Table 1

- Interpolation used for frequencies between standard 1/3-octave values

- Free-field listening conditions (frontal incidence, binaural)

## References

( )
