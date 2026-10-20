# Internal WORLD vocoder pipeline

Orchestrate Harvest and CheapTrick on pre-loaded raw PCM signal and
compute normalized linear power from the CheapTrick spectrogram.

## Usage

``` r
.voxit_world_pipeline(wave, fs, frame_period = 5)
```

## Arguments

- wave:

  Numeric vector; raw PCM samples in the range -1 to 1

- fs:

  Integer; original audio sample rate (Hz)

- frame_period:

  Numeric; frame shift in milliseconds (default 5ms)

## Value

List with two elements:

- f0_parameter: list(temporal_positions, f0, vuv)

- spectrum_parameter: list(temporal_positions, fs, linPower)

## Details

D4C aperiodicity is not computed: no voxit feature reads it (only F0 and
linPower are consumed), and it is costly on long recordings. The full
CheapTrick spectrogram is dropped once linPower has been derived from
it.
