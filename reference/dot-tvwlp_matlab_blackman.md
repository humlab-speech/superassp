# Symmetric Blackman window matching MATLAB blackman(N) av::blackman is a different (DFT-even) variant; SRH/SEDREAMS need the symmetric form: w(n) = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1)).

Symmetric Blackman window matching MATLAB blackman(N) av::blackman is a
different (DFT-even) variant; SRH/SEDREAMS need the symmetric form: w(n)
= 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1)).

## Usage

``` r
.tvwlp_matlab_blackman(N)
```
