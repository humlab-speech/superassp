# Build a geom class from a ggplot2 parent

ggplot2 is a suggested dependency, so ggproto classes cannot be created
at build time; they are constructed per layer instead.

## Usage

``` r
.assp_geom(name, parent)
```

## Arguments

- name:

  Character. Name of the new ggproto class.

- parent:

  Character. Name of the exported ggplot2 parent class.

## Value

A ggproto object.
