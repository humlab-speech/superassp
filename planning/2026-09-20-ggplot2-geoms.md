# ggplot2 geoms for the audio data objects — 2026-09-20

**Status:** implemented (2026-09-20). Verification numbers in §5, deviations in §6.

## 1. The ask

The package could produce plots only through `ggtrack()` — 160 lines in
`R/ggtrack.R` that were never exported (`git grep ggtrack NAMESPACE` was empty),
so no plotting entry point existed for users at all, and the `_pkgdown.yml`
reference index listed `ggtrack`, `get_track_label` and `get_track_label_expr`
as if they were public. The ask: a ggplot2-compliant geom that handles the data
objects the package produces.

## 2. What the objects are

| object | produced by | shape |
|---|---|---|
| `AsspDataObj` | every `trk_*`, `read_audio()`, `read_ssff()`, `read_track()` | named list of equally spaced tracks; one column = a time series (`RMS[dB]`, `audio`), n columns = coefficients (`Fi[Hz]`, `DFT[dB]`, `LPCi`); `sampleRate` = record rate, `startTime`, `origFreq` = source sample rate |
| `JsonTrackObj` | `lst_*(return_jstf = TRUE)`, `read_jstf()` | time-bounded slices with named scalar/vector fields |
| `AVAudio` | internal only | S7 raw sample buffer; not user-reachable |

`as.data.frame.AsspDataObj()` already flattens tracks to a wide table
(`frame_time` + one column per coefficient) and attaches `track_labels`,
`track_descriptions`, `sampleRate`, `startTime`; it now also attaches `origFreq`
(§4.4).

## 3. Constraints found while probing (ggplot2 4.0.3, R 4.6)

These decided the architecture; they are not obvious from the ggplot2 docs:

1. **`fortify()` is the only object hook.** `ggplot(data = obj)` and
   `layer(data = obj)` both call `fortify()`, whose default rejects anything
   whose `dim()` is not length 2. Without a method, `ggplot(obj)` errors with
   *"`data` must be a <data.frame> …"*.
2. **A layer never sees the track columns.** In 4.0 `Layer$compute_aesthetics()`
   evaluates the mapping against the layer data and hands the stat/geom only the
   resulting aesthetic columns (measured: a probe stat's `setup_data()` received
   `x, PANEL, group` for a 1026-column spectrum layer). Any wide → long rewrite
   in `setup_data()` is therefore impossible.
3. **Inherited data can be prepared with a data function.**
   `Layer$layer_data()` accepts `data = function(plot_data) …` (and
   `fortify.function()` explicitly allows it), applied to the plot data at build
   time — that is where an inherited `AsspDataObj` is still whole.
4. **`after_stat()` labels are stripped**: `make_labels()` unwraps
   `after_stat(x)` to `x`, so default axis labels stay clean.
5. **`GeomRaster` is the fast path for spectra**: one rasterGrob per panel
   (1025×404 for a 4 s DFT at 10 ms shift) instead of 400k tiles, and it warns
   loudly if the grid is uneven — the package's spectra are evenly spaced by
   construction.
6. **A units-assigned column breaks ggplot2 scaling.** `scale_type.units()`
   aborts with *"Variable of class 'units' found, but 'units' package is not
   attached"*, so `as.data.frame(obj)` (whose default is `convert_units = TRUE`)
   cannot be plotted unless the caller attaches units. `fortify()` therefore
   passes `convert_units = FALSE`.

## 4. What landed

### 4.1 `R/ggtrack_geoms.R` (new)

`geom_track()` (GeomLine-derived `GeomTrack`) and `geom_spectrogram()`
(GeomRaster-derived `GeomSpectrogram`), both ordinary `ggplot2::layer()`
constructions, so scales, facets, coordinates, legends (`draw_key` inherited)
and themes apply unchanged. Both accept an `AsspDataObj`, the `as.data.frame()`
table, or an already long table, and rewrite to a long table before the layer is
built (§3.2): `frame_time`, `value`, `track` (cleaned column name),
`band` (track template), `bin` (coefficient index, `NA` when not indexed), plus
`freq` for spectrograms. `mapping = NULL` means "draw the tracks"
(`colour = track` and a legend); a supplied mapping refers to that long table,
and `inherit.aes = FALSE` keeps plot-level aesthetics from referring to columns
the rewrite removed. `tracks =` selects by column name (`"F1_Hz"`), template
(`"Fi[Hz]"`) or object track name; `na.zeros = TRUE` turns stored zeros into line
breaks; a track plot of more than 64 columns warns and points at
`geom_spectrogram()`.

### 4.2 Frequency axis (§ "faithful or nothing")

SSFF stores a spectrum from 0 Hz to the Nyquist rate, so a track with `n`
columns puts coefficient `j` at `(j - 1) * origFreq / (2 * (n - 1))`.
Measured with a 1 kHz sine (44.1 kHz, 1 s, written to a temp WAV):

| algorithm | bins | peak bin | Nyquist grid | `origFreq / ncol` |
|---|---|---|---|---|
| `trk_dft_spectrum()` | 1025 | 47 | **990 Hz** | 2022 Hz |
| `trk_css_spectrum()` | 1025 | 47 | **990 Hz** | 2022 Hz |
| `trk_lps_spectrum()` | 1025 | 49 | **1034 Hz** | 2108 Hz |

The `origFreq / ncol` recipe that both spectrum examples documented was twice
the truth; those examples are fixed (§4.5). `freq_hz_per_bin` overrides the
derivation, and `data` without `origFreq` aborts with that instruction.

### 4.3 `fortify()` and `ggtrack()`

`fortify.AsspDataObj` / `fortify.JsonTrackObj` are registered in `.onLoad`
(`R/zzz.R`) with `registerS3method()` into ggplot2's namespace — the pattern the
file already used for `print`/`summary` — because the generic lives in a
suggested package. `ggtrack()` is now exported and additionally labels the axes
of the new workflow: no mapping + one track → that track's label
(`ggtrack(f0) + geom_track()` gives "fo \[Hz\]"), the waveform track →
"Amplitude", `frame_time` → "Time \[s\]".

### 4.4 Export surface and metadata

`geom_*` plus `ggtrack`, `get_track_label`, `get_track_label_expr` are exported
(they were documented, pkgdown-referenced and unreachable before);
`tests/testthat/test-export-policy.R`, `CLAUDE.md` and `_pkgdown.yml` follow.
`as.data.frame.AsspDataObj()` gains the `origFreq` attribute. Version bumped to
3.2.0 with a `NEWS.md` section.

### 4.5 Doc-example bugs fixed

`trk_lps_spectrum()`'s example plotted `res[["CSS[dB]"]]` — missing, so R
silently plotted `NULL` — and both it and `trk_css_spectrum()` used the
doubled frequency axis.

## 5. Verification

`tests/testthat/test-gg-geoms.R` (13 tests, all passing) pins: panel data equals
the object's values and times; the inherited path equals the supplied path
byte-for-byte; track selection by all three name forms; the legend keys; zero →
NA; the spectrogram grid (range `0 … origFreq/2`, spacing `bin_hz`, and every
row's `value` equal to `matrix[frame, bin]`); the raster build and the wide-table
path; the abort messages (no spectral track, two spectral tracks, no
`origFreq`, JSTF as tracks); `fortify()` output and column-wise plotting; and
`ggtrack()` labels in both subscript and plain-text modes.

Rendered output was checked on the drawn grobs, not just the built data:
1/2/3/4/4 selected tracks → 1/2/3/4/4 paths in the panel, waveform → 1 path,
overlay of `geom_spectrogram()` + `geom_track(data = audio)` → 1 rasterGrob
(1025×404) + 1 path. Five figures are rendered to `/tmp/ggproof/*.png`
(`/tmp/render.R`); the terminal-browser pane could not be opened here
(`herdr pane split` rejects `--right-click`), so the figures were not shown in a
pane.

## 6. Deviations and notes

* `geom_spectrogram()` requires the `origFreq` convention above; multi-column
  tracks of anything else (formants, LPC coefficients) are rejected for the
  raster and belong to `geom_track()`. There is deliberately no track-type
  registry — the frequency axis is derived from `origFreq` and the column count.
* `JsonTrackObj` gets `fortify()` (and `CLAUDE.md` shows the `geom_segment()`
  recipe) but no geom: slices are bounded intervals, not equally spaced tracks.
* The long table is row-per-frame-per-coefficient, so a long recording's
  spectrogram is a large data frame (12 M rows for 60 s at 5 ms shift);
  decimation and windowing of the object are the caller's job.
* `inherit.aes = FALSE` on both layers is a deliberate divergence from the
  ggplot2 default, forced by the rewrite (§3.2) — same choice `geom_sf()` makes.
