# Sampling grids

The sampling grid defines the regular pixel grid onto which HEALPix data is resampled before being rendered by Matplotlib. It controls the spatial extent, pixel resolution, and output image shape.

Three ways to define the target raster:

### 1. Parametrised (dict or `ParametrizedSamplingGrid`)

Pass a dict to `sampling_grid`. Recognised keys:

| Key          | Default  | Description                                                                           |
| ------------ | -------- | ------------------------------------------------------------------------------------- |
| `shape`      | `1024`   | Output array size. An `int` produces a square grid; a 2-tuple sets `(width, height)`. |
| `resolution` | inferred | Step size in degrees. A `float` expands to equal x/y steps.                           |
| `center`     | inferred | `(lon, lat)` centre of the grid in degrees.                                           |

```python
sampling_grid = {"shape": (2048, 1024), "center": (0.0, 0.0)}
```

The spatial extent and pixel resolution are inferred automatically from the bounding box of your `cell_ids`.

### 2. Bounding box (`ParametrizedSamplingGrid.from_bbox`)

Use `SamplingGrid.from_bbox()` to pin the output to a fixed region regardless of the data:

```python
sampling_grid = healpix_plotting.SamplingGrid.from_bbox(
    bbox=(5, 36.0, 13.0, 45.0),  # (lon_min, lat_min, lon_max, lat_max) in degrees
    shape=512,
)
```

This is the right choice when comparing multiple datasets or animating over time.

### 3. Affine transform (`AffineSamplingGrid`)

Use `AffineSamplingGrid` when the output pixels must align with a reference raster (e.g. a GeoTIFF):

```python
from healpix_plotting.sampling_grid import AffineSamplingGrid
from affine import Affine

transform = Affine(0.01, 0, 5, 0, -0.01, 45)  # 0.01 deg/pixel, top-left at (-10, 60)
sampling_grid = AffineSamplingGrid.from_transform(transform, shape=(4000, 2500))
```

:::{seealso}
See {doc}`../tutorials/quickstart` for more information
:::
