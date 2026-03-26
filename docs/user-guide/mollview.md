# `healpix_mollview` — HEALPix Mollweide Visualisation

> A `healpy`-free replacement for `healpy.mollview` and `healpy.gnomview`,
> built on top of **healpix-geo** for all coordinate-system geometry.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation & Dependencies](#2-installation--dependencies)
3. [Quick Start](#3-quick-start)
4. [Public API](#4-public-api)
   - 4.1 [`mollview`](#41-mollview)
   - 4.2 [`mollgnomview`](#42-mollgnomview)
5. [Figure Layout — `hold` and `sub`](#5-figure-layout--hold-and-sub)
6. [RING vs NESTED pixel ordering](#6-ring-vs-nested-pixel-ordering)
7. [Ellipsoid support](#7-ellipsoid-support)
8. [Rendering backends](#8-rendering-backends)
   - 8.1 [Fast path (default)](#81-fast-path-default)
   - 8.2 [Cartopy path (`coastlines=True`)](#82-cartopy-path-coastlinestrue)
9. [Algorithm — projection math](#9-algorithm--projection-math)
   - 9.1 [Mollweide projection](#91-mollweide-projection)
   - 9.2 [Gnomonic projection](#92-gnomonic-projection)
10. [Internal helpers](#10-internal-helpers)
11. [Comparison with `healpy`](#11-comparison-with-healpy)
12. [Common recipes](#12-common-recipes)
13. [Known limitations](#13-known-limitations)

---

## 1. Overview

`mollview.py` provides two visualisation functions for HEALPix sky/sphere maps:

| Function       | Projection                           | healpy equivalent |
| -------------- | ------------------------------------ | ----------------- |
| `mollview`     | Mollweide (equal-area, full sky)     | `healpy.mollview` |
| `mollgnomview` | Gnomonic (tangent-plane, local zoom) | `healpy.gnomview` |

Key differences from healpy:

- **No healpy dependency.** All HEALPix geometry is handled by
  [`healpix-geo`](https://github.com/EOPF-DGGS/healpix-geo).
- **Depth is inferred automatically** from the map size — you never pass it
  explicitly.
- **RING order by default**, like healpy. Pass `nest=True` for NESTED maps.
- **Non-spherical ellipsoids** (e.g. WGS84) are supported via healpix-geo.
- **No return value** — the function draws into the current matplotlib state,
  exactly like `healpy.mollview`.
- **`hold` and `sub` parameters** control where the plot appears, mirroring
  healpy's interface.
- **Two rendering backends**, selected automatically:
  - _Fast path_ (default): pure matplotlib, no cartopy required, ~70 ms.
  - _Cartopy path_: activated by `coastlines=True`, ~500 ms.

---

## 2. Installation & Dependencies

**Always required:**

```bash
pip install healpix-geo numpy matplotlib
```

**Optional — only needed when `coastlines=True`:**

```bash
pip install cartopy
```

`cartopy` is imported **lazily**: it is never loaded unless you explicitly pass
`coastlines=True`. If cartopy is not installed and `coastlines=False` (the
default), the module works normally. If cartopy is absent and `coastlines=True`,
a clear `ImportError` is raised with installation instructions.

Python >= 3.10 is required (uses `X | Y` union type hints).

---

## 3. Quick Start

```python
import numpy as np
import matplotlib.pyplot as plt
from mollview import mollview, mollgnomview

# --- Synthetic RING map (depth 5, nside 32) ---
depth = 5
npix  = 12 * 4**depth          # 12 288 pixels
m     = np.random.default_rng(0).standard_normal(npix)

# Fast full-sky view (pure matplotlib, no cartopy)
mollview(m, title="My map", cmap="RdBu_r", unit="K")
plt.show()

# With coastlines (requires cartopy)
mollview(m, title="With coastlines", coastlines=True,
         coastline_kwargs={"linewidth": 0.8, "edgecolor": "cyan"})
plt.show()

# Local zoom centred on (lon=45deg, lat=30deg)
mollgnomview(m, lon_center=45.0, lat_center=30.0, fov_deg=20.0,
             title="Zoom 20deg", cmap="plasma")
plt.show()
```

---

## 4. Public API

### 4.1 `mollview`

```python
mollview(
    hpx_map,
    nest             = False,
    title            = "",
    cmap             = "viridis",
    vmin             = None,
    vmax             = None,
    rot              = 0.0,
    ellipsoid        = "sphere",
    graticule        = True,
    graticule_step   = 30.0,
    unit             = "",
    bgcolor          = "black",
    n_lon            = 1800,
    n_lat            = 900,
    norm             = None,
    bad_color        = "gray",
    flip             = "geo",
    figsize          = (14, 7),
    colorbar         = True,
    hold             = False,
    sub              = None,
    coastlines       = False,
    coastline_kwargs = None,
)
```

Renders a HEALPix map in the **Mollweide equal-area projection** (full sky).

The rendering backend is chosen **automatically**:

| `coastlines`      | Backend         | cartopy required | Typical time |
| ----------------- | --------------- | ---------------- | ------------ |
| `False` (default) | Pure matplotlib | No               | ~70 ms       |
| `True`            | Cartopy         | Yes              | ~500 ms      |

#### Parameters

| Parameter          | Type                                | Default     | Description                                                                                                                                       |
| ------------------ | ----------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hpx_map`          | `np.ndarray`, shape `(12·4^depth,)` | —           | Input HEALPix map. RING order by default; use `nest=True` for NESTED.                                                                             |
| `nest`             | `bool`                              | `False`     | Pixel ordering. `False` = RING (healpy default), `True` = NESTED.                                                                                 |
| `title`            | `str`                               | `""`        | Title displayed above the map.                                                                                                                    |
| `cmap`             | `str` or `Colormap`                 | `"viridis"` | Matplotlib colormap. `"RdBu_r"` is recommended for CMB/temperature maps.                                                                          |
| `vmin`             | `float` or `None`                   | `None`      | Lower bound of the colour scale. Defaults to the 2nd percentile of finite values.                                                                 |
| `vmax`             | `float` or `None`                   | `None`      | Upper bound of the colour scale. Defaults to the 98th percentile of finite values.                                                                |
| `rot`              | `float`                             | `0.0`       | Central longitude of the map in degrees. `rot=180` centres on the anti-meridian.                                                                  |
| `ellipsoid`        | `str`                               | `"sphere"`  | Reference ellipsoid for healpix-geo. `"sphere"` gives results identical to healpy. Other values: `"WGS84"`, `"GRS80"`.                            |
| `graticule`        | `bool`                              | `True`      | Draw meridians and parallels.                                                                                                                     |
| `graticule_step`   | `float`                             | `30.0`      | Spacing of graticule lines in degrees.                                                                                                            |
| `unit`             | `str`                               | `""`        | Unit string shown below the colorbar.                                                                                                             |
| `bgcolor`          | `str`                               | `"black"`   | Background colour outside the Mollweide ellipse.                                                                                                  |
| `n_lon`            | `int`                               | `1800`      | Number of sample columns in the internal raster grid. Increase for high-depth maps (depth >= 8).                                                  |
| `n_lat`            | `int`                               | `900`       | Number of sample rows in the internal raster grid.                                                                                                |
| `norm`             | `Normalize` or `None`               | `None`      | Custom matplotlib normalisation (e.g. `LogNorm()`). Overrides `vmin`/`vmax`.                                                                      |
| `bad_color`        | `str`                               | `"gray"`    | Colour for `NaN` values.                                                                                                                          |
| `flip`             | `str`                               | `"geo"`     | East/west convention. `"geo"`: east to the right (default). `"astro"`: east to the left (astronomical convention).                                |
| `figsize`          | `(float, float)`                    | `(14, 7)`   | Figure size in inches. Only used when a new figure is created.                                                                                    |
| `colorbar`         | `bool`                              | `True`      | Show a horizontal colorbar below the map.                                                                                                         |
| `hold`             | `bool`                              | `False`     | If `True`, draw into the current axes. Ignored when `sub` is provided.                                                                            |
| `sub`              | `(int, int, int)` or `None`         | `None`      | `(nrows, ncols, index)` subplot position. Overrides `hold`.                                                                                       |
| `coastlines`       | `bool`                              | `False`     | Overlay Natural Earth coastlines. Activates the cartopy backend — cartopy must be installed.                                                      |
| `coastline_kwargs` | `dict` or `None`                    | `None`      | Extra kwargs forwarded to `ax.add_feature(COASTLINE, ...)`. Only used when `coastlines=True`. Example: `{"linewidth": 0.8, "edgecolor": "cyan"}`. |

#### Returns

`None`. Access the current figure with `plt.gcf()`.

#### Raises

| Exception     | Condition                                       |
| ------------- | ----------------------------------------------- |
| `ValueError`  | `hpx_map.size` is not of the form `12·4^depth`. |
| `ValueError`  | `flip` is not `"astro"` or `"geo"`.             |
| `ImportError` | `coastlines=True` but cartopy is not installed. |

---

### 4.2 `mollgnomview`

```python
mollgnomview(
    hpx_map,
    lon_center,
    lat_center,
    nest             = False,
    fov_deg          = 10.0,
    title            = "",
    cmap             = "viridis",
    vmin             = None,
    vmax             = None,
    ellipsoid        = "sphere",
    unit             = "",
    n_lon            = 800,
    n_lat            = 800,
    figsize          = (7, 7),
    colorbar         = True,
    hold             = False,
    sub              = None,
    coastlines       = False,
    coastline_kwargs = None,
)
```

Renders a local zoom in the **gnomonic (tangent-plane) projection**.
The depth is inferred automatically from the map size.

The same backend selection logic as `mollview` applies: `coastlines=False`
uses the fast matplotlib path; `coastlines=True` activates cartopy.

#### Parameters specific to `mollgnomview`

| Parameter    | Type    | Default | Description                                   |
| ------------ | ------- | ------- | --------------------------------------------- |
| `lon_center` | `float` | —       | Longitude of the view centre in degrees.      |
| `lat_center` | `float` | —       | Latitude of the view centre in degrees.       |
| `fov_deg`    | `float` | `10.0`  | Total field of view (square side) in degrees. |

All other parameters (`nest`, `cmap`, `vmin`, `vmax`, `ellipsoid`, `unit`,
`n_lon`, `n_lat`, `figsize`, `colorbar`, `hold`, `sub`, `coastlines`,
`coastline_kwargs`) behave identically to those of `mollview`.

---

## 5. Figure Layout — `hold` and `sub`

The three modes below are mutually exclusive, with priority order:
**`sub` > `hold=True` > `hold=False`**.

### New figure (default)

```python
mollview(m, title="My map")
plt.show()
```

A fresh `plt.figure()` is created on each call. Same behaviour as
`healpy.mollview`.

### `hold=True` — reuse the current axes

```python
plt.figure(figsize=(14, 7))
mollview(m, hold=True, title="Overlay")
plt.show()
```

The map is drawn into whichever axes is currently active.

### `sub=(nrows, ncols, idx)` — subplot grid

Place multiple maps in a single figure using standard matplotlib subplot
indexing (1-based).

```python
fig = plt.figure(figsize=(18, 5), facecolor="black")
mollview(m1, sub=(1, 2, 1), title="Map A", cmap="plasma")
mollview(m2, sub=(1, 2, 2), title="Map B", cmap="RdBu_r")
plt.tight_layout()
plt.savefig("comparison.png", dpi=150, bbox_inches="tight", facecolor="black")
plt.close()
```

More complex layouts:

```python
plt.figure(figsize=(18, 10), facecolor="black")
mollview(m1, sub=(2, 3, 1), title="(1,1)")
mollview(m2, sub=(2, 3, 2), title="(1,2)")
mollview(m3, sub=(2, 3, 3), title="(1,3)")
mollview(m4, sub=(2, 3, 4), title="(2,1)")
mollview(m5, sub=(2, 3, 5), title="(2,2)")
mollview(m6, sub=(2, 3, 6), title="(2,3)")
plt.tight_layout()
plt.show()
```

---

## 6. RING vs NESTED pixel ordering

HEALPix maps can be stored in two pixel orderings:

| Ordering   | Description                                          | Default in            |
| ---------- | ---------------------------------------------------- | --------------------- |
| **RING**   | Pixels ordered in iso-latitude rings, west to east   | `healpy`, this module |
| **NESTED** | Pixels ordered along a space-filling (Z-order) curve | `healpix-geo`         |

This module always uses RING order by default (`nest=False`), matching
`healpy.mollview`.

```python
# RING map (healpy default)
mollview(m_ring)

# NESTED map
mollview(m_nested, nest=True)
```

---

## 7. Ellipsoid support

The `ellipsoid` parameter is forwarded directly to
`healpix_geo.nested.lonlat_to_healpix`. On a sphere (`ellipsoid="sphere"`,
the default), results are numerically identical to healpy. On a non-spherical
ellipsoid, the authalic latitude is used for the lon/lat → cell-ID conversion.

```python
# Geographic data referenced to WGS84
mollview(m_ring, ellipsoid="WGS84", title="WGS84 Mollweide")
```

> **Note:** changing the ellipsoid only affects which HEALPix cell each
> image pixel maps to. The visual shape of the Mollweide projection is always
> the same mathematical ellipse regardless of the ellipsoid choice.

---

## 8. Rendering backends

The backend is selected **automatically** based on the `coastlines` parameter.

### 8.1 Fast path (default)

Activated when `coastlines=False` (the default). **cartopy is not imported.**

The grid is sampled uniformly in **Mollweide (x, y) space** using the
analytic inverse formula. The resulting `data_img` array is therefore already
in display coordinates, and `ax.imshow()` requires no `transform` argument and
no render-time reprojection.

```
Uniform grid in Mollweide (x, y) space   [n_lon × n_lat]
         |
         v  _mollweide_inverse()          — analytic, vectorised
  (lon, lat) for each sample point        [degrees]
         |
         v  lonlat_to_healpix()           — healpix-geo
  cell_id
         |
         v  hpx_map[cell_id]
  data_img  [n_lat × n_lon]              — already in display space
         |
         v  ax.imshow()                  — no transform, no reprojection
  Final figure
```

The graticule is drawn by pre-projecting meridians and parallels with
`_mollweide_forward()` (Newton iteration) and calling `ax.plot()` directly.
This replaces `ax.gridlines()`, which alone costs ~400 ms in cartopy.

Typical time (nside=64, 1800×900): **~70 ms**.

### 8.2 Cartopy path (`coastlines=True`)

Activated when `coastlines=True`. **cartopy is imported lazily** at this point.
If cartopy is not installed, a clear `ImportError` is raised.

The grid must be sampled uniformly in **PlateCarree (lon, lat) space** so that
it matches the `extent=[-180, 180, -90, 90]` declared to
`imshow(transform=PlateCarree())`. Sampling in Mollweide space and passing a
geographic extent would misalign the raster with the coastline overlay (the
root cause of the artefact visible in earlier versions).

```
Uniform grid in PlateCarree (lon, lat)   [n_lon × n_lat]
         |
         v  lonlat_to_healpix()           — healpix-geo
  cell_id
         |
         v  hpx_map[cell_id]
  data_img  [n_lat × n_lon]              — in PlateCarree space
         |
         v  GeoAxes.imshow(transform=PlateCarree())
            cartopy reprojects to Mollweide at render time
         |
         v  ax.add_feature(COASTLINE)     — Natural Earth overlay
         |
         v  pre-projected graticule       — _mollweide_forward(), no ax.gridlines()
  Final figure
```

Typical time (nside=64, 1800×900): **~500 ms** (vs ~1.1 s for a naive
all-cartopy implementation that also calls `ax.gridlines()`).

The remaining overhead is unavoidable when coastlines are needed: `GeoAxes`
initialisation (~200–400 ms) and cartopy's render-time reprojection of the
image (~200 ms) are both required by `add_feature()`.

---

## 9. Algorithm — projection math

### 9.1 Mollweide projection

**Inverse formula** (Mollweide (x, y) → lon/lat, used for the fast-path grid):

The standard Mollweide ellipse spans `x ∈ [-2√2, +2√2]`, `y ∈ [-√2, +√2]`.
Points outside satisfy `x²/8 + y²/2 > 1` and are masked as NaN.

The auxiliary angle `θ = arcsin(y / √2)`, then:

```
sin(lat) = (2θ + sin(2θ)) / π
lon      = lon_0 + π·x / (2√2·cos(θ))
```

**Forward formula** (lon/lat → Mollweide (x, y), used for graticule lines):

Requires Newton iteration to solve `2θ + sin(2θ) = π·sin(lat)`:

```
x = (2√2 / π) · (lon − lon_0) · cos(θ)
y = √2 · sin(θ)
```

Convergence is reached in fewer than 10 iterations to `tol=1e-9`.

### 9.2 Gnomonic projection

The gnomonic (tangent-plane) projection is centred on `(lon_center, lat_center)`.
Sampling coordinates `(x, y)` are in degrees of arc from the tangent point.
The inverse formula (tangent-plane → lon/lat):

```
c  = arctan(rho),    rho = sqrt(x^2 + y^2)

lat = arcsin( cos(c)·sin(lat_0) + y·sin(c)·cos(lat_0)/rho )
lon = lon_0 + arctan2( x·sin(c),  rho·cos(lat_0)·cos(c) − y·sin(lat_0)·sin(c) )
```

For the fast path the gnomonic grid is also sampled directly in tangent-plane
space, so `ax.imshow()` requires no `transform`. The cartopy path samples in
PlateCarree space and uses a `Gnomonic` GeoAxes, for the same reason as
`mollview` (coastline alignment).

---

## 10. Internal helpers

| Function                     | Signature                                                                | Description                                                                               |
| ---------------------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------- |
| `_depth_from_npix`           | `(npix) → int`                                                           | Infers HEALPix depth. Raises `ValueError` for invalid sizes.                              |
| `_mollweide_inverse`         | `(x, y, central_lon_rad) → (lon_deg, lat_deg, valid)`                    | Analytic inverse Mollweide. Returns a boolean mask for points inside the ellipse.         |
| `_mollweide_forward`         | `(lon_deg, lat_deg, central_lon_rad) → (x, y)`                           | Forward Mollweide via Newton iteration. Used for graticule drawing.                       |
| `_gnomonic_inverse`          | `(x, y, lon0_deg, lat0_deg) → (lon_deg, lat_deg)`                        | Inverse gnomonic projection. `x`, `y` in degrees of arc.                                  |
| `_rasterise`                 | `(hpx_map, depth, lon_grid, lat_grid, valid, nest, ellipsoid) → ndarray` | HEALPix lookup for valid grid points. Returns `(H, W)` float64 with NaN outside the mask. |
| `_build_cmap`                | `(cmap, bad_color) → Colormap`                                           | Builds the colormap object and sets the bad-value colour.                                 |
| `_build_norm`                | `(data_img, vmin, vmax, norm) → Normalize`                               | Builds the normalisation object (2nd/98th percentile defaults).                           |
| `_add_colorbar`              | `(fig, ax, norm, cmap_obj, unit, text_color)`                            | Adds a horizontal colorbar.                                                               |
| `_make_axes_plain`           | `(hold, sub, figsize, bgcolor) → (fig, ax)`                              | Creates a plain `Axes` (fast path).                                                       |
| `_make_axes_geo`             | `(crs, hold, sub, figsize, bgcolor) → (fig, ax)`                         | Creates a cartopy `GeoAxes` (cartopy path).                                               |
| `_draw_graticule_moll`       | `(ax, step_deg, central_lon_rad, ...)`                                   | Pre-projected graticule on a plain Axes (fast path).                                      |
| `_draw_graticule_on_geoaxes` | `(ax, step_deg, central_lon_rad, ...)`                                   | Pre-projected graticule on a GeoAxes (cartopy path, avoids `ax.gridlines()`).             |
| `_finalize_moll_axes`        | `(ax)`                                                                   | Adds the oval boundary, clips the image to the ellipse, and sets axis limits.             |
| `_require_cartopy`           | `() → (ccrs, cfeature)`                                                  | Lazy cartopy import with a clear error message if not installed.                          |

---

## 11. Comparison with `healpy`

| Feature              | `healpy.mollview`             | `mollview` (this module)                             |
| -------------------- | ----------------------------- | ---------------------------------------------------- |
| Pixel ordering       | RING by default               | RING by default (`nest=False`)                       |
| Depth/nside          | inferred from map length      | inferred from map length                             |
| Rotation             | `rot=(lon, lat, psi)` 3-tuple | `rot=lon_deg` scalar                                 |
| Ellipsoid            | sphere only                   | `"sphere"`, `"WGS84"`, `"GRS80"`, …                  |
| Image resolution     | fixed internal grid           | configurable: `n_lon`, `n_lat`                       |
| East/west default    | `"astro"` (east left)         | `"geo"` (east right)                                 |
| Return value         | `None`                        | `None`                                               |
| `hold` / `sub`       | supported                     | supported                                            |
| Coastlines           | not supported                 | `coastlines=True` (requires cartopy)                 |
| healpy dependency    | required                      | **not required**                                     |
| cartopy dependency   | not required                  | optional (only for `coastlines=True`)                |
| NaN handling         | `bad_color`                   | `bad_color`                                          |
| Colour normalisation | `min`/`max`                   | percentiles 2/98 by default; custom `norm` supported |

### Migration from healpy

```python
# healpy
import healpy as hp
hp.mollview(m, title="My map", nest=False, cmap="RdBu_r", unit="K",
            min=-3, max=3)

# This module — drop-in equivalent
from mollview import mollview
mollview(m, title="My map", nest=False, cmap="RdBu_r", unit="K",
         vmin=-3, vmax=3)
```

Parameter renames: `min` → `vmin`, `max` → `vmax`, `width_px`/`height_px` → `n_lon`/`n_lat`.

The `rot` parameter accepts only a longitude scalar (not a 3-tuple).

The `flip` default is `"geo"` here (east to the right), whereas healpy defaults
to `"astro"` (east to the left). Add `flip="astro"` to reproduce the healpy
convention.

---

## 12. Common recipes

### Save to file (no display)

```python
import matplotlib
matplotlib.use("Agg")

from mollview import mollview
import matplotlib.pyplot as plt

mollview(m, title="My map", cmap="RdBu_r")
plt.savefig("map.png", dpi=150, bbox_inches="tight", facecolor="black")
plt.close()
```

### Astronomical convention (east to the left)

```python
mollview(m, flip="astro", title="Astro convention")
```

### With coastlines (geographic data)

```python
# Requires: pip install cartopy
mollview(m, flip="geo", ellipsoid="WGS84",
         coastlines=True,
         coastline_kwargs={"linewidth": 0.8, "edgecolor": "white"},
         title="Geographic — WGS84")
```

### Custom colour normalisation (log scale)

```python
import matplotlib.colors as mcolors
from mollview import mollview

mollview(m_positive, norm=mcolors.LogNorm(vmin=1e-3, vmax=1.0),
         title="Log scale", cmap="inferno")
```

### Symmetric diverging scale

```python
import numpy as np
from mollview import mollview

absmax = np.nanpercentile(np.abs(m), 98)
mollview(m, vmin=-absmax, vmax=absmax, cmap="RdBu_r",
         title="Symmetric +/-{:.2f}".format(absmax))
```

### Compare two maps side by side

```python
import matplotlib.pyplot as plt
from mollview import mollview

fig = plt.figure(figsize=(18, 5), facecolor="black")
mollview(m1, sub=(1, 2, 1), title="Map A", cmap="plasma", vmin=-3, vmax=3)
mollview(m2, sub=(1, 2, 2), title="Map B", cmap="plasma", vmin=-3, vmax=3)
plt.tight_layout()
plt.savefig("comparison.png", dpi=150, bbox_inches="tight", facecolor="black")
plt.close()
```

### Full-sky + local zoom in one figure

```python
import matplotlib.pyplot as plt
from mollview import mollview, mollgnomview

fig = plt.figure(figsize=(18, 8), facecolor="black")
mollview(m,     sub=(1, 2, 1), title="Full sky",   cmap="RdBu_r")
mollgnomview(m, lon_center=45.0, lat_center=30.0,
             fov_deg=20.0, sub=(1, 2, 2), title="Zoom 20deg", cmap="RdBu_r")
plt.tight_layout()
plt.show()
```

### Increase raster resolution for a high-depth map

```python
# depth=8 -> nside=256 -> 786 432 pixels
# Default 1800x900 may be too coarse; use 3600x1800
mollview(m_high_res, n_lon=3600, n_lat=1800,
         title="High-res map (depth=8)")
```

---

## 13. Known limitations

- **`rot` is a scalar longitude only.** Unlike `healpy.mollview` which accepts
  a 3-tuple `(lon, lat, psi)` for full rotation, only longitude rotation is
  supported. Latitude rotation and roll are not implemented.

- **No partial-sky maps.** The input must be a full-sky map of exactly
  `12·4^depth` pixels. Partial-sky maps must be zero-padded or filled with
  `NaN` to full size before plotting.

- **Performance.** The HEALPix lookup is vectorised over all `n_lon × n_lat`
  sample points with a single `lonlat_to_healpix` call. For
  `n_lon × n_lat >= 4 × 10^6` this may take a few seconds on CPU. Reduce
  the resolution or pre-cache the index grid if you need to render many maps
  at the same resolution.

- **Coastlines require cartopy.** When `coastlines=True` the cartopy backend
  adds ~430 ms of fixed overhead (GeoAxes init + render-time reprojection)
  that is unavoidable because `add_feature()` requires a `GeoAxes`. This
  overhead is independent of the map resolution.

- **`flip="astro"` on the cartopy path** reverses the column order of the
  raster before passing it to cartopy. Cartopy is not aware of this flip, so
  the coastlines are drawn in the correct (geographic) orientation while the
  data is mirrored. Use `flip="geo"` (the default) with `coastlines=True`
  for correct coastline alignment.
