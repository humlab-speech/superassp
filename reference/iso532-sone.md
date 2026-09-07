# ISO 532 Sone (Loudness) Conversions

Functions for converting between phon (loudness level) and sone
(loudness) according to ISO 532 standards. Sone is a linear unit of
perceived loudness, where 1 sone = 40 phons by definition.

## Details

### Standards

Two ISO 532 methods are supported:

- **ISO 532-1 (Zwicker method)**: Uses piecewise mathematical formulas

- **ISO 532-2 (Moore-Glasberg method)**: Uses lookup table with
  interpolation

### Zwicker Method (ISO 532-1)

The Zwicker method uses different formulas depending on loudness level:

**For sone \< 1 (phon \< 40):** \$\$L_N = 40 \times S^{0.35}\$\$

**For sone \>= 1 (phon \>= 40):** \$\$L_N = 40 + 10 \times \log_2(S)\$\$

Inverse formulas:

**For phon \< 40:** \$\$S = (L_N / 40)^{1/0.35} = (L_N / 40)^{2.857}\$\$

**For phon \>= 40:** \$\$S = 2^{(L_N - 40) / 10}\$\$

### Moore-Glasberg Method (ISO 532-2)

Uses a lookup table with 23 reference points from 0.001 sone (0 phon) to
337.6 sone (120 phon), with log-linear interpolation.

### Key Properties

- **Reference:** 1 sone = 40 phons = 40 dB SPL at 1 kHz

- **Doubling:** Each 10 phon increase ~= doubles loudness in sones

- **Linear perception:** Sones represent linear loudness (2 sones = 2x
  louder)

- **Stevens' power law:** Loudness is proportional to intensity^0.3

## References

( )

( )

(Stevens 1936)
