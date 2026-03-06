# healpix-plotting

> **Fast, practical static plots of [HEALPix](https://healpix.sourceforge.io/) data with [matplotlib](https://matplotlib.org/) and [cartopy](https://scitools.org.uk/cartopy/) — with ellipsoidal (WGS84) support for geoscience workflows.**

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.13%2B-blue)](https://www.python.org/)

---

## Overview

`healpix-plotting` prioritises **getting a usable figure quickly** over perfectly accurate cell-geometry rendering. It rasterises HEALPix data via nearest-neighbour resampling onto a regular lon/lat grid and renders the result with Cartopy's `imshow`.

Unlike astronomy-focused HEALPix tools, this library is built with **Earth observation and geoscience in mind**: the underlying coordinate operations are provided by [healpix-geo](https://healpix-geo.readthedocs.io/en/latest/), which supports geodetically correct reference ellipsoids such as WGS84.

The library is well suited for:

- Exploratory analysis and quality control of EO / climate data
- Debugging and sanity checks
- Quick snapshots for reports and discussions

It is **not** intended as a true boundary-polygon renderer.

---

## Dependencies

| Package                                                      | Role                                                  |
| ------------------------------------------------------------ | ----------------------------------------------------- |
| [healpix-geo](https://healpix-geo.readthedocs.io/en/latest/) | Core HEALPix ↔ lon/lat conversions, ellipsoid support |
| [cartopy](https://scitools.org.uk/cartopy/)                  | Map projections and rendering                         |
| [matplotlib](https://matplotlib.org/)                        | Figure/axes backend                                   |
| [numpy](https://numpy.org/)                                  | Array operations                                      |
| [numpy-groupies](https://github.com/ml31415/numpy-groupies)  | Aggregation during duplicate-cell-id deduplication    |
| [affine](https://github.com/rasterio/affine)                 | Affine transform support for `AffineSamplingGrid`     |
| [scipy](https://scipy.org/)                                  | Reserved for future bilinear interpolation            |

---

## Installation

```bash
pip install healpix-plotting
```

---

## How it works

1. **Build a target sampling grid** — a regular lon/lat grid inferred from the data, defined by a bounding box, or specified via an affine transform.
2. **Resample** — each sampling point is mapped to a HEALPix cell id and filled by nearest-neighbour lookup (with optional aggregation for duplicate ids).
3. **Render** — the resulting raster is drawn on a Cartopy axis with `imshow(..., transform=PlateCarree())`.

Because the library rasterises via nearest-neighbour resampling, it does **not** attempt exact polygon boundary filling.

## Ellipsoidal support

Standard HEALPix was originally defined on the unit sphere. For geoscience applications —
such as Earth observation, numerical weather prediction, or climate modelling — coordinates
are commonly expressed in geodetic lon/lat on a reference ellipsoid (most often WGS84).

`healpix-plotting` delegates coordinate conversions to [healpix-geo](https://healpix-geo.readthedocs.io/en/latest/),
which implements ellipsoidal HEALPix conversions. Passing `ellipsoid="WGS84"` to `HealpixGrid`
propagates this choice through all internal operations.

For a full explanation of how ellipsoidal HEALPix works, see the
[healpix-geo ellipsoids tutorial](https://healpix-geo.readthedocs.io/en/latest/tutorials/).

---

## Things to keep in mind

HEALPix is a spherical tessellation. Depending on your use case you may need to consider:

- **True cell boundaries** — cells are not lon/lat rectangles; boundaries are curved on the sphere.
- **Projection effects** — any map projection changes apparent cell geometry.
- **Resolution and aliasing** — the target sampling-grid resolution determines sharpness and artifact level.
- **Dateline and polar behaviour** — extent wrapping and polar distortion can introduce discontinuities.
- **Sphere vs. ellipsoid** — for precise geoscience work, prefer `ellipsoid="WGS84"` over the default `"sphere"`.

---

## Limitations

- Rendering is raster-based; exact HEALPix cell boundaries are **not** drawn.
- Output quality depends on the sampling-grid resolution (speed vs. aliasing trade-off).
- `"bilinear"` interpolation is currently **not implemented** (raises `NotImplementedError`).
- The helper pattern `healpix_grid.operations.healpix_to_lonlat(cell_ids, **healpix_grid.as_keyword_params())` does **not** work for `indexing_scheme="zuniq"` — see the caveat in [HealpixGrid](#healpixgrid).

---

## Alternatives

| Library                                                      | Best for                                                                            |
| ------------------------------------------------------------ | ----------------------------------------------------------------------------------- |
| [**healpy**](https://healpy.readthedocs.io/)                 | Conventional HEALPix map visualisation; `mollview` and friends. Sphere only.        |
| [**earthkit-plots**](https://earthkit-plots.readthedocs.io/) | Publication-quality figures in the ECMWF / earthkit stack.                          |
| [**xdggs**](https://xdggs.readthedocs.io/)                   | Analysis/selection workflows and polygon-based rendering with true cell boundaries. |

### Choosing a tool

```
Need a quick geoscience plot (EO, NWP, climate)?  → healpix-plotting  ✓  (WGS84 support)
Working in the ECMWF/earthkit ecosystem?           → earthkit-plots    ✓
Need exact cell boundaries / polygon operations?   → xdggs             ✓
Standard healpy full-sky (astronomy) workflow?     → healpy            ✓
```
---

## Implementation notes

### Nearest-neighbour resampling

For each sampling-grid point, the corresponding HEALPix cell id is computed. Ids not present
in the source data are masked; the remaining values are placed into the raster using
`searchsorted`-based indexing.

### Aggregation

Before resampling, duplicate `cell_ids` are collapsed with `numpy-groupies` using the
function specified by `agg` (default: `"mean"`). The sorted unique ids are then used
as the lookup table for the raster fill.

### Cartopy rendering

Plotting uses `imshow` with `transform=ccrs.PlateCarree()` and `interpolation="nearest"`.
For non-global subsets, the extent is set **before** plotting to obtain a smoother result
with Cartopy.

