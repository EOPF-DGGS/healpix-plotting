# `healpix_mollview` — HEALPix Mollweide Visualisation

> A `healpy`-free replacement for `healpy.mollview` and `healpy.gnomview`,
> built on top of **healpix-geo** for all coordinate-system geometry.

## Overview

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
  - _Cartopy path_: activated by `coastlines=True`, ~500 ms.\*

## Figure Layout — `hold` and `sub`

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

## RING vs NESTED pixel ordering

HEALPix maps can be stored in two pixel orderings:

| Ordering   | Description                                          | Default in            |
| ---------- | ---------------------------------------------------- | --------------------- |
| **RING**   | Pixels ordered in iso-latitude rings, west to east   | `healpy`, this module |
| **NESTED** | Pixels ordered along a space-filling (Z-order) curve | `healpix-geo`         |

This module always uses RING order by default (`nest=False`), matching `healpy.mollview`.

```python
# RING map (healpy default)
mollview(m_ring)

# NESTED map
mollview(m_nested, nest=True)
```

## Ellipsoid support

The `ellipsoid` parameter is forwarded directly to `healpix_geo.nested.lonlat_to_healpix`. On a sphere (`ellipsoid="sphere"`, the default), results are numerically identical to healpy. On a non-spherical ellipsoid, the authalic latitude is used for the lon/lat → cell-ID conversion.

```python
# Geographic data referenced to WGS84
mollview(m_ring, ellipsoid="WGS84", title="WGS84 Mollweide")
```

> **Note:** changing the ellipsoid only affects which HEALPix cell each
> image pixel maps to. The visual shape of the Mollweide projection is always
> the same mathematical ellipse regardless of the ellipsoid choice.

## Rendering backends

The backend is selected **automatically** based on the `coastlines` parameter.

### Fast path (default)

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

### Cartopy path (`coastlines=True`)

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

## Comparison with `healpy`

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

hp.mollview(m, title="My map", nest=False, cmap="RdBu_r", unit="K", min=-3, max=3)

# This module — drop-in equivalent
from mollview import mollview

mollview(m, title="My map", nest=False, cmap="RdBu_r", unit="K", vmin=-3, vmax=3)
```

Parameter renames: `min` → `vmin`, `max` → `vmax`, `width_px`/`height_px` → `n_lon`/`n_lat`.

The `rot` parameter accepts only a longitude scalar (not a 3-tuple).

The `flip` default is `"geo"` here (east to the right), whereas healpy defaults
to `"astro"` (east to the left). Add `flip="astro"` to reproduce the healpy
convention.

## Known limitations

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
