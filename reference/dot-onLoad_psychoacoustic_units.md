# Initialize Psychoacoustic Units

Registers psychoacoustic units (Bark, ERB, mel, semitone) with the units
package. This function is called automatically when superassp is loaded
if the units package is available.

## Usage

``` r
.onLoad_psychoacoustic_units()
```

## Value

NULL (called for side effects)

## Details

Registers:

- Bark - Critical band rate scale

- ERB - Equivalent Rectangular Bandwidth rate scale

- mel - Mel scale (perceptual pitch)

- semitone - Musical semitone (ST)
