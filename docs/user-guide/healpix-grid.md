# The HEALPix grid

healpix-plot uses [healpix-geo](https://healpix-geo.readthedocs.io) internally for all
HEALPix coordinate conversions. The healpix-geo documentation covers everything you need to
know about the HEALPix grid: resolution levels, indexing schemes (nested, ring, zuniq),
ellipsoids, and coordinate conversions.

| Parameter         | Type          | Description                                                                                                                                           |
| ----------------- | ------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `level`           | `int`         | HEALPix depth. Must be in **[0, 29]**.                                                                                                                |
| `indexing_scheme` | `str`         | Cell numbering scheme: `"nested"`, `"ring"`, or `"zuniq"`.                                                                                            |
| `ellipsoid`       | `str` or dict | Reference ellipsoid: `"sphere"` (default), `"WGS84"`, or a custom dict with `radius` (sphere) or `semimajor_axis` + `inverse_flattening` (ellipsoid). |

## The HealpixGrid object

In healpix-plot, you describe your data with a `HealpixGrid` object that wraps
the healpix-geo parameters:

```python
import healpix_plot

grid = healpix_plot.HealpixGrid(
    level=4,  # resolution level, 0-29
    indexing_scheme="nested",  # "nested", "ring", or "zuniq"
    ellipsoid="sphere",  # default; "WGS84" for geodesic accuracy
)
```

It exposes a `.operations` property that delegates directly to healpix-geo:

```python
import numpy as np

cell_ids = np.arange(12 * 4**grid.level, dtype="uint64")
lon, lat = grid.operations.healpix_to_lonlat(cell_ids, **grid.as_keyword_params())
```

For everything else — what levels mean, how indexing schemes differ, ellipsoid support —
refer to the [healpix-geo documentation](https://healpix-geo.readthedocs.io).
