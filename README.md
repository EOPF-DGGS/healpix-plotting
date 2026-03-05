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

| Package | Role |
|---|---|
| [healpix-geo](https://healpix-geo.readthedocs.io/en/latest/) | Core HEALPix ↔ lon/lat conversions, ellipsoid support |
| [cartopy](https://scitools.org.uk/cartopy/) | Map projections and rendering |
| [matplotlib](https://matplotlib.org/) | Figure/axes backend |
| [numpy](https://numpy.org/) | Array operations |
| [numpy-groupies](https://github.com/ml31415/numpy-groupies) | Aggregation during duplicate-cell-id deduplication |
| [affine](https://github.com/rasterio/affine) | Affine transform support for `AffineSamplingGrid` |
| [scipy](https://scipy.org/) | Reserved for future bilinear interpolation |

---

## Installation
```bash
pip install healpix-plotting
```

---

## How it works

1. **Build a target sampling grid** — a regular lon/lat grid inferred from the data, defined by a bounding box, or specified via an affine transform (see [Sampling grid variants](#sampling-grid-variants)).
2. **Resample** — each sampling point is mapped to a HEALPix cell id and filled by nearest-neighbour lookup (with optional aggregation for duplicate ids).
3. **Render** — the resulting raster is drawn on a Cartopy axis with `imshow(..., transform=PlateCarree())`.

Because the library rasterises via nearest-neighbour resampling, it does **not** attempt exact polygon boundary filling.

---

### Default example (sphere)
```python
import numpy as np
import healpix_plotting as hpplt

# Define the HEALPix grid on the unit sphere
healpix_grid = hpplt.HealpixGrid(level=6, indexing_scheme="nested", ellipsoid="sphere")

# All cell ids at this level (global)
cell_ids = np.arange(12 * 4 ** healpix_grid.level, dtype="uint64")

# Synthetic data as a function of lon/lat
lon, lat = healpix_grid.operations.healpix_to_lonlat(
    cell_ids, **healpix_grid.as_keyword_params()
)
data = np.cos(8 * np.deg2rad(lon)) * np.sin(4 * np.deg2rad(lat))

# Plot — returns a matplotlib AxesImage mappable
mappable = hpplt.plot(
    cell_ids,
    data,
    healpix_grid=healpix_grid,
    sampling_grid={"shape": 1024},
    projection="Mollweide",
    colorbar=True,
)
```

### WGS84 example (recommended for geoscience)

For Earth observation and geoscience data, use `ellipsoid="WGS84"` to perform coordinate
conversions on the geodetically correct reference ellipsoid rather than a sphere.
This matters most at higher resolutions and higher latitudes.
See the [healpix-geo ellipsoids tutorial](https://healpix-geo.readthedocs.io/en/latest/tutorials/)
for background on why this matters.
```python
import numpy as np
import healpix_plotting as hpplt

# Use WGS84 — the standard geodetic ellipsoid for Earth observation
healpix_grid = hpplt.HealpixGrid(level=6, indexing_scheme="nested", ellipsoid="WGS84")

cell_ids = np.arange(12 * 4 ** healpix_grid.level, dtype="uint64")

# lon/lat are now geodetic coordinates on WGS84
lon, lat = healpix_grid.operations.healpix_to_lonlat(
    cell_ids, **healpix_grid.as_keyword_params()
)
data = np.cos(8 * np.deg2rad(lon)) * np.sin(4 * np.deg2rad(lat))

mappable = hpplt.plot(
    cell_ids,
    data,
    healpix_grid=healpix_grid,
    sampling_grid={"shape": 1024},
    projection="Mollweide",
    colorbar=True,
)
```

> **Note:** The examples above use `indexing_scheme="nested"`. The schemes `"ring"` and `"zuniq"` are also valid — see [HealpixGrid](#healpixgrid) below.

---

## `HealpixGrid`

```python
hpplt.HealpixGrid(level, indexing_scheme, ellipsoid="sphere")
```

| Parameter | Type | Description |
|---|---|---|
| `level` | `int` | HEALPix depth. Must be in **[0, 29]**. |
| `indexing_scheme` | `str` | Cell numbering scheme: `"nested"`, `"ring"`, or `"zuniq"`. |
| `ellipsoid` | `str` or dict | Reference ellipsoid: `"sphere"` (default), `"WGS84"`, or a custom dict with `radius` (sphere) or `semimajor_axis` + `inverse_flattening` (ellipsoid). |

`healpix_grid.as_keyword_params()` returns `{"depth": level, "ellipsoid": ellipsoid}` for direct unpacking into `healpix-geo` calls.

> **`"zuniq"` caveat:** `zuniq.healpix_to_lonlat` does not accept a `depth` argument, so the pattern
> `healpix_grid.operations.healpix_to_lonlat(cell_ids, **healpix_grid.as_keyword_params())`
> will fail when `indexing_scheme="zuniq"`. Pass only `ellipsoid` explicitly in that case:
> ```python
> lon, lat = healpix_grid.operations.healpix_to_lonlat(cell_ids, ellipsoid=healpix_grid.ellipsoid)
> ```

---

## Key parameters of `plot()`

`plot()` returns a `matplotlib.image.AxesImage` (the mappable). Use `.axes` to access the underlying `Axes` object.

| Parameter | Description |
|---|---|
| `cell_ids` | `numpy.ndarray` of cell ids describing spatial positions. |
| `data` | 1-D array for scalar data (color-mapped), or 2-D array of shape `(N, 3)` / `(N, 4)` for RGB / RGBA. |
| `healpix_grid` | A `HealpixGrid` instance (or equivalent dict). |
| `sampling_grid` | Target raster resolution and extent. Pass a dict such as `{"shape": 1024}`; missing `center` / `resolution` are inferred from the data. See [Sampling grid variants](#sampling-grid-variants). |
| `projection` | A Cartopy CRS name (e.g. `"Mollweide"`) or an actual CRS object. Unknown names raise a `ValueError`. |
| `agg` | Aggregation function applied when `cell_ids` contains duplicates before resampling. Accepted values: `"mean"` (default), `"median"`, `"std"`, `"var"`, `"min"`, `"max"`, `"first"`, `"last"`. |
| `interpolation` | Resampling method. `"nearest"` (default and only implemented). `"bilinear"` is accepted by the API but raises `NotImplementedError`. |
| `background_value` | Fill value for grid points with no matching cell id. Default: `numpy.nan`. |
| `ax` | An existing Cartopy `Axes` to draw on. If omitted, a new figure is created using `projection`. |
| `title` | Optional string title for the axes. |
| `cmap` | Colormap (name or `Colormap` object). Default: `"viridis"`. |
| `vmin`, `vmax` | Scalar data range for colour normalisation. |
| `norm` | A `matplotlib.colors.Normalize` instance for finer colour control. |
| `colorbar` | `True` to add a colorbar, or a dict of kwargs forwarded to `figure.colorbar()`. Default: `False`. |
| `axis_labels` | `None` (default, uses `"Longitude"` / `"Latitude"`), a dict with `"x"` / `"y"` keys, or `"none"` to suppress labels entirely. |

---

## Sampling grid variants

Three ways to define the target raster:

### 1. Parametrised (dict or `ParametrizedSamplingGrid`)

Pass a dict to `sampling_grid`. Recognised keys:

| Key | Default | Description |
|---|---|---|
| `shape` | `1024` | Output array size. An `int` produces a square grid; a 2-tuple sets `(width, height)`. |
| `resolution` | inferred | Step size in degrees. A `float` expands to equal x/y steps. |
| `center` | inferred | `(lon, lat)` centre of the grid in degrees. |

```python
sampling_grid={"shape": (2048, 1024), "center": (0.0, 0.0)}
```

### 2. Bounding box (`ParametrizedSamplingGrid.from_bbox`)

```python
from healpix_plotting import SamplingGrid
from healpix_plotting.sampling_grid import ParametrizedSamplingGrid

grid = ParametrizedSamplingGrid.from_bbox(
    bbox=(-30.0, 35.0, 45.0, 75.0),  # (xmin, ymin, xmax, ymax) in degrees
    shape=512,
)
hpplt.plot(..., sampling_grid=grid)
```

### 3. Affine transform (`AffineSamplingGrid`)

For raster-aligned outputs (e.g. matching an existing GeoTIFF grid):

```python
from affine import Affine
from healpix_plotting.sampling_grid import AffineSamplingGrid

transform = Affine(0.1, 0, -180, 0, -0.1, 90)  # 0.1° resolution, global
grid = AffineSamplingGrid.from_transform(transform, shape=(3600, 1800))
hpplt.plot(..., sampling_grid=grid)
```

---

## RGB / RGBA data

`plot()` accepts multi-band data directly. If `data` has shape `(N, 3)` (RGB) or `(N, 4)` (RGBA), the values are composited per-pixel rather than colour-mapped. Colourmap parameters (`cmap`, `vmin`, `vmax`, `norm`) are ignored in this mode.

```python
rgb = np.stack([red, green, blue], axis=-1)  # shape (N, 3), values in [0, 1]
hpplt.plot(cell_ids, rgb, healpix_grid=healpix_grid, sampling_grid={"shape": 1024})
```

---

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

## Recipe: overlay true HEALPix boundaries

For sanity checks or presentations you can combine fast raster plots from `healpix-plotting`
with boundary polygons from [`xdggs`](https://xdggs.readthedocs.io/):
```python
import xdggs
import cartopy.crs as ccrs

info = xdggs.HealpixInfo(
    level=healpix_grid.level,
    indexing_scheme=healpix_grid.indexing_scheme,
)
polys = info.cell_boundaries(cell_ids_subset)  # shapely polygons

# plot() returns an AxesImage; use .axes to get the Axes object
ax = hpplt.plot(...).axes
ax.add_geometries(
    polys,
    crs=ccrs.PlateCarree(),
    facecolor="none",
    edgecolor="k",
    linewidth=0.2,
)
```

---
## Alternatives

| Library | Best for |
|---|---|
| [**healpy**](https://healpy.readthedocs.io/) | Conventional HEALPix map visualisation; `mollview` and friends. Sphere only. |
| [**earthkit-plots**](https://earthkit-plots.readthedocs.io/) | Publication-quality figures in the ECMWF / earthkit stack. |
| [**xdggs**](https://xdggs.readthedocs.io/) | Analysis/selection workflows and polygon-based rendering with true cell boundaries. |

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

---

## License

Apache License 2.0 — see [LICENSE](LICENSE) for details.
