# SSFF and JSTF: reading, writing, and round-tripping

`trk_*` functions write **SSFF** (Simple Signal File Format) tracks;
`lst_*` functions write **JSTF** (JSON Track Format) summaries. Both
round-trip through matching `read_*`/`write_*` pairs and both load
straight into `emuR`.

## SSFF: time-series tracks

Write a track to disk, then read it back with
[`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md):

``` r

wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
out_dir <- tempdir()

trk_pitch_rapt(wav, toFile = TRUE, outputDirectory = out_dir, explicitExt = "f0")
#> Applying `trk_pitch_rapt()` to 1 recording
#> Successfully processed 1 of 1 file

f0_path <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(wav)), ".f0"))
f0_obj <- read_ssff(f0_path)
track_names(f0_obj)
#> [1] "f0"
sample_rate(f0_obj)
#> [1] 100
```

[`write_ssff()`](https://humlab-speech.github.io/superassp/reference/write_ssff.md)
is the inverse operation, for writing an `AsspDataObj` you built or
modified in R:

``` r

f0_obj[["f0"]][1:5, ] <- 0  # zero out the first 5 frames
out_path <- tempfile(fileext = ".f0")
write_ssff(f0_obj, out_path)
read_ssff(out_path)[["f0"]][1:5, ]
#> [1] 0 0 0 0 0
```

[`read_track()`](https://humlab-speech.github.io/superassp/reference/read_track.md)/[`write_track()`](https://humlab-speech.github.io/superassp/reference/write_track.md)
are format-agnostic dispatchers: pass any `trk_*` output path and they
detect SSFF vs. JSTF from the file’s own metadata rather than the file
extension.

## JSTF: summary measures

`lst_*` functions produce a `JsonTrackObj` — a self-describing container
with a field schema and one “slice” per analysis window:

``` r

vr <- lst_voice_report(wav, toFile = FALSE, return_jstf = TRUE)
names(vr$field_schema)
#>  [1] "start_time"            "end_time"              "selection_start"      
#>  [4] "selection_end"         "median_pitch"          "mean_pitch"           
#>  [7] "sd_pitch"              "min_pitch"             "max_pitch"            
#> [10] "num_pulses"            "num_periods"           "mean_period"          
#> [13] "sd_period"             "fraction_unvoiced"     "num_voice_breaks"     
#> [16] "degree_voice_breaks"   "jitter_local_percent"  "jitter_local_abs"     
#> [19] "jitter_rap_percent"    "jitter_ppq5_percent"   "jitter_ddp_percent"   
#> [22] "shimmer_local_percent" "shimmer_local_db"      "shimmer_apq3_percent" 
#> [25] "shimmer_apq5_percent"  "shimmer_apq11_percent" "shimmer_dda_percent"  
#> [28] "mean_autocorrelation"  "mean_nhr"              "mean_hnr"
vr$slices[[1]]$values
#> $start_time
#> [1] 0
#> 
#> $end_time
#> [1] 4.035374
#> 
#> $selection_start
#> [1] 0
#> 
#> $selection_end
#> [1] 4.035374
#> 
#> $median_pitch
#> [1] 120.3934
#> 
#> $mean_pitch
#> [1] 120.6306
#> 
#> $sd_pitch
#> [1] 5.429101
#> 
#> $min_pitch
#> [1] 109.6691
#> 
#> $max_pitch
#> [1] 194.3577
#> 
#> $num_pulses
#> [1] 312
#> 
#> $num_periods
#> [1] 310
#> 
#> $mean_period
#> [1] 0.008302197
#> 
#> $sd_period
#> [1] 0.0002832662
#> 
#> $fraction_unvoiced
#> [1] 0.3557772
#> 
#> $num_voice_breaks
#> [1] 1
#> 
#> $degree_voice_breaks
#> [1] 0.007434255
#> 
#> $jitter_local_percent
#> [1] 0.5258646
#> 
#> $jitter_local_abs
#> [1] 4.365831e-05
#> 
#> $jitter_rap_percent
#> [1] 0.2692026
#> 
#> $jitter_ppq5_percent
#> [1] 0.2436183
#> 
#> $jitter_ddp_percent
#> [1] 0.8076078
#> 
#> $shimmer_local_percent
#> [1] 4.226066
#> 
#> $shimmer_local_db
#> [1] 0.3440395
#> 
#> $shimmer_apq3_percent
#> [1] 1.774378
#> 
#> $shimmer_apq5_percent
#> [1] 2.545493
#> 
#> $shimmer_apq11_percent
#> [1] 4.347712
#> 
#> $shimmer_dda_percent
#> [1] 5.323134
#> 
#> $mean_autocorrelation
#> [1] NaN
#> 
#> $mean_nhr
#> [1] 0.005467156
#> 
#> $mean_hnr
#> [1] 22.62239
```

Write it to a `.jstf` file and read it back:

``` r

jstf_path <- tempfile(fileext = ".jstf")
write_jstf(vr, jstf_path)
vr_reloaded <- read_jstf(jstf_path)
identical(vr$field_schema, vr_reloaded$field_schema)
#> [1] TRUE
```

## Building a `JsonTrackObj` from your own results

If you’re wrapping a measure that isn’t already a `lst_*` function,
[`create_json_track_obj()`](https://humlab-speech.github.io/superassp/reference/create_json_track_obj.md)
builds the same self-describing structure that
[`write_jstf()`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)
expects. It is an internal helper (not exported — see the package’s
Export Policy), so package contributors writing a new `lst_*` wrapper
call it via `:::`:

``` r

obj <- superassp:::create_json_track_obj(
  results = list(f0_mean = 150.2, f0_sd = 18.4),
  function_name = "my_custom_summary",
  file_path = wav,
  sample_rate = 44100,
  audio_duration = 1.5
)
write_jstf(obj, tempfile(fileext = ".jstf"))
```

## Loading into emuR

Both
[`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md)-produced
`AsspDataObj`s and
[`read_jstf()`](https://humlab-speech.github.io/superassp/reference/read_jstf.md)-produced
`JsonTrackObj`s carry the sample rate, start time, and track names emuR
expects — write with `toFile = TRUE` into your emuR database’s `_ssff`
directory structure and the tracks load without further conversion.
