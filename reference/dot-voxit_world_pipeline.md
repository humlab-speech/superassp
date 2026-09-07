# Internal WORLD vocoder pipeline

Orchestrate Harvest, CheapTrick, D4C on pre-loaded raw PCM signal.
Computes normalized linear power from CheapTrick spectrogram.

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

List with three elements:

- f0_parameter: list(temporal_positions, f0, vuv)

- spectrum_parameter: list(spectrogram, temporal_positions, fs,
  linPower)

- source_parameter: list(aperiodicity, temporal_positions, fs)
