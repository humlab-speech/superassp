# Create a JsonTrackObj

Create a JsonTrackObj

## Usage

``` r
create_json_track_obj(
  results,
  function_name,
  file_path,
  sample_rate = NULL,
  audio_duration = NULL,
  beginTime = 0,
  endTime = 0,
  parameters = list()
)
```

## Arguments

- results:

  List of results from lst\_\* function

- function_name:

  Name of the DSP function

- file_path:

  Original audio file path

- sample_rate:

  Audio sample rate in Hz

- audio_duration:

  Total audio duration in seconds

- beginTime:

  Start time in seconds

- endTime:

  End time in seconds

- parameters:

  Function parameters used

## Value

Object of class JsonTrackObj
