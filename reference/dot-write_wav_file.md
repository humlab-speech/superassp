# Write WAV File

Internal function to write audio samples to WAV file format.

## Usage

``` r
.write_wav_file(filename, samples, sample_rate, channels)
```

## Arguments

- filename:

  Character; output file path

- samples:

  Integer vector; audio samples (s32le format)

- sample_rate:

  Integer; sample rate in Hz

- channels:

  Integer; number of channels
