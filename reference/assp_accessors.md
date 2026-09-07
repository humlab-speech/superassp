# Accessor methods for AsspDataObj and JsonTrackObj

Read-only accessors that work on both
[AsspDataObj](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
(SSFF signal data) and
[JsonTrackObj](https://humlab-speech.github.io/superassp/reference/JsonTrackObj.md)
(JSTF summary objects).

|                      |                                  |
|----------------------|----------------------------------|
| Function             | Returns                          |
| `sample_rate(x)`     | Sample rate in Hz                |
| `n_records(x)`       | Number of records / slices       |
| `signal_duration(x)` | Duration in seconds              |
| `start_time(x)`      | Start time of first record/slice |
| `track_names(x)`     | Track (field) names              |
| `file_path(x)`       | Source audio file path           |
| `track_formats(x)`   | Data-type string per track       |

## Usage

``` r
sample_rate(x, ...)

# S3 method for class 'AsspDataObj'
sample_rate(x, ...)

# S3 method for class 'JsonTrackObj'
sample_rate(x, ...)

n_records(x, ...)

# S3 method for class 'AsspDataObj'
n_records(x, ...)

# S3 method for class 'JsonTrackObj'
n_records(x, ...)

signal_duration(x, ...)

# S3 method for class 'AsspDataObj'
signal_duration(x, ...)

# S3 method for class 'JsonTrackObj'
signal_duration(x, ...)

start_time(x, ...)

# S3 method for class 'AsspDataObj'
start_time(x, ...)

# S3 method for class 'JsonTrackObj'
start_time(x, ...)

track_names(x, ...)

# S3 method for class 'AsspDataObj'
track_names(x, ...)

# S3 method for class 'JsonTrackObj'
track_names(x, ...)

file_path(x, ...)

# S3 method for class 'AsspDataObj'
file_path(x, ...)

# S3 method for class 'JsonTrackObj'
file_path(x, ...)

track_formats(x, ...)

# S3 method for class 'AsspDataObj'
track_formats(x, ...)

# S3 method for class 'JsonTrackObj'
track_formats(x, ...)

# S3 method for class 'AsspDataObj'
tracks(x, ...)

# S3 method for class 'AsspDataObj'
dur(x, ...)

# S3 method for class 'AsspDataObj'
numRecs(x, ...)

# S3 method for class 'AsspDataObj'
rate(x, ...)

# S3 method for class 'AsspDataObj'
startTime(x, ...)

# S3 method for class 'JsonTrackObj'
tracks(x, ...)

# S3 method for class 'JsonTrackObj'
dur(x, ...)

# S3 method for class 'JsonTrackObj'
rate(x, ...)

# S3 method for class 'JsonTrackObj'
numRecs(x, ...)

# S3 method for class 'JsonTrackObj'
startTime(x, ...)

dur(x, ...)

numRecs(x, ...)

rate(x, ...)

startTime(x, ...)

tracks(x, ...)
```

## Arguments

- x:

  An
  [AsspDataObj](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
  or
  [JsonTrackObj](https://humlab-speech.github.io/superassp/reference/JsonTrackObj.md).

- ...:

  Currently unused; reserved for future extensions.

## Value

- `sample_rate(x)` — numeric scalar, sample rate in Hz.

- `n_records(x)` — integer scalar, number of records (SSFF) or slices
  (JSTF).

- `signal_duration(x)` — numeric scalar, duration in seconds.

- `start_time(x)` — numeric scalar, start time of the first record/slice
  in seconds.

- `track_names(x)` — character vector of track (field) names.

- `file_path(x)` — character scalar, source audio file path.

- `track_formats(x)` — character vector, one storage-type string (e.g.
  `"INT16"`, `"REAL64"`) per track.

## Deprecated aliases (since v2.8.0)

|                |                      |
|----------------|----------------------|
| Old            | New                  |
| `rate(x)`      | `sample_rate(x)`     |
| `numRecs(x)`   | `n_records(x)`       |
| `dur(x)`       | `signal_duration(x)` |
| `startTime(x)` | `start_time(x)`      |
| `tracks(x)`    | `track_names(x)`     |

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
au <- read_audio(wav)

sample_rate(au)      # 44100
#> [1] 44100
n_records(au)        # sample count
#> [1] 177960
signal_duration(au)  # seconds
#> [1] 4.035374
start_time(au)       # 0
#> [1] 0
track_names(au)      # "audio"
#> [1] "audio"
track_formats(au)    # "INT16"
#> [1] "INT16"
basename(file_path(au))
#> [1] "a1.wav"
```
