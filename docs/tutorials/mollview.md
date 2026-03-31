---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  name: python3
  display_name: Python 3
---

# Mollweide visualisation tutorial

This tutorial walks through the main use cases for `mollview` and `mollgnomview`, from a minimal first plot to advanced layouts and recipes.

## Setup

**Always required:**

```{code-cell} python
pip install healpix-geo numpy matplotlib
```

**Optional — only needed when `coastlines=True`:**

```{code-cell} python
pip install cartopy
```

`cartopy` is imported **lazily**: it is never loaded unless you explicitly pass `coastlines=True`. If cartopy is not installed and `coastlines=False` (the default), the module works normally. If cartopy is absent and `coastlines=True`, a clear `ImportError` is raised with installation instructions.

Python >= 3.10 is required (uses `X | Y` union type hints).

```{code-cell} python
import numpy as np
import matplotlib.pyplot as plt
from mollview import mollview, mollgnomview
```

## Quick start

```{code-cell} python
# --- Synthetic RING map (depth 5, nside 32) ---
depth = 5
npix = 12 * 4**depth  # 12 288 pixels
m = np.random.default_rng(0).standard_normal(npix)

# Fast full-sky view (pure matplotlib, no cartopy)
mollview(m, title="My map", cmap="RdBu_r", unit="K")
plt.show()
```

```{code-cell} python
# With coastlines (requires cartopy)
mollview(
    m,
    title="With coastlines",
    coastlines=True,
    coastline_kwargs={"linewidth": 0.8, "edgecolor": "cyan"},
)
plt.show()
```

```{code-cell} python
# Local zoom centred on (lon=45deg, lat=30deg)
mollgnomview(
    m, lon_center=45.0, lat_center=30.0, fov_deg=20.0, title="Zoom 20deg", cmap="plasma"
)
plt.show()
```

---

## Recipes

### Save to file (no display)

```{code-cell} python
import matplotlib
matplotlib.use("Agg")

mollview(m, title="My map", cmap="RdBu_r")
plt.savefig("map.png", dpi=150, bbox_inches="tight", facecolor="black")
plt.close()
```

### Astronomical convention (east to the left)

```{code-cell} python
mollview(m, flip="astro", title="Astro convention")
```

### With coastlines (geographic data)

```{code-cell} python
# Requires: pip install cartopy
mollview(
    m,
    flip="geo",
    ellipsoid="WGS84",
    coastlines=True,
    coastline_kwargs={"linewidth": 0.8, "edgecolor": "white"},
    title="Geographic — WGS84",
)
```

### Custom colour normalisation (log scale)

```{code-cell} python
import matplotlib.colors as mcolors

mollview(
    m_positive,
    norm=mcolors.LogNorm(vmin=1e-3, vmax=1.0),
    title="Log scale",
    cmap="inferno",
)
```

### Symmetric diverging scale

```{code-cell} python
absmax = np.nanpercentile(np.abs(m), 98)
mollview(
    m,
    vmin=-absmax,
    vmax=absmax,
    cmap="RdBu_r",
    title="Symmetric +/-{:.2f}".format(absmax),
)
```

### Compare two maps side by side

```{code-cell} python
fig = plt.figure(figsize=(18, 5), facecolor="black")
mollview(m1, sub=(1, 2, 1), title="Map A", cmap="plasma", vmin=-3, vmax=3)
mollview(m2, sub=(1, 2, 2), title="Map B", cmap="plasma", vmin=-3, vmax=3)
plt.tight_layout()
plt.savefig("comparison.png", dpi=150, bbox_inches="tight", facecolor="black")
plt.close()
```

### Full-sky + local zoom in one figure

```{code-cell} python
fig = plt.figure(figsize=(18, 8), facecolor="black")
mollview(m, sub=(1, 2, 1), title="Full sky", cmap="RdBu_r")
mollgnomview(
    m,
    lon_center=45.0,
    lat_center=30.0,
    fov_deg=20.0,
    sub=(1, 2, 2),
    title="Zoom 20deg",
    cmap="RdBu_r",
)
plt.tight_layout()
plt.show()
```

### Multi-panel layout (6 maps)

```{code-cell} python
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

### Increase raster resolution for a high-depth map

```{code-cell} python
# depth=8 -> nside=256 -> 786 432 pixels
# Default 1800x900 may be too coarse; use 3600x1800
mollview(m_high_res, n_lon=3600, n_lat=1800, title="High-res map (depth=8)")
```
