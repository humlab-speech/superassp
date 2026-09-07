# Internal: Analyze Pharyngeal Voice Quality from Time Range

Core pharyngeal analysis function using pladdrr. Analyzes spectral
measures (H1-H2, H1-A1, etc.) at vowel onset and optionally midpoint.

## Usage

``` r
analyze_pharyngeal_times(
  sound,
  start_time,
  end_time,
  min_pitch_initial = 50,
  max_pitch_initial = 800
)
```

## Arguments

- sound:

  pladdrr Sound object

- start_time:

  Start time of interval (seconds)

- end_time:

  End time of interval (seconds)

- min_pitch_initial:

  Initial minimum pitch (default: 50 Hz)

- max_pitch_initial:

  Initial maximum pitch (default: 800 Hz)

## Value

Named list with 68 pharyngeal measures
