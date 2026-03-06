# The plot() function

{func}`healpix_plotting.plot` is the main entry point of the library. It resamples HEALPix data onto a regular pixel grid and renders it with Matplotlib and Cartopy.

`plot()` returns a `matplotlib.image.AxesImage` (the mappable). Use `.axes` to access the underlying `Axes` object.

| Parameter          | Description                                                                                                                                                                                   |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cell_ids`         | `numpy.ndarray` of cell ids describing spatial positions.                                                                                                                                     |
| `data`             | 1-D array for scalar data (color-mapped), or 2-D array of shape `(N, 3)` / `(N, 4)` for RGB / RGBA.                                                                                           |
| `healpix_grid`     | A `HealpixGrid` instance (or equivalent dict).                                                                                                                                                |
| `sampling_grid`    | Target raster resolution and extent. Pass a dict such as `{"shape": 1024}`; missing `center` / `resolution` are inferred from the data.                                                       |
| `projection`       | A Cartopy CRS name (e.g. `"Mollweide"`) or an actual CRS object. Unknown names raise a `ValueError`.                                                                                          |
| `agg`              | Aggregation function applied when `cell_ids` contains duplicates before resampling. Accepted values: `"mean"` (default), `"median"`, `"std"`, `"var"`, `"min"`, `"max"`, `"first"`, `"last"`. |
| `interpolation`    | Resampling method. `"nearest"` (default and only implemented). `"bilinear"` is accepted by the API but raises `NotImplementedError`.                                                          |
| `background_value` | Fill value for grid points with no matching cell id. Default: `numpy.nan`.                                                                                                                    |
| `ax`               | An existing Cartopy `Axes` to draw on. If omitted, a new figure is created using `projection`.                                                                                                |
| `title`            | Optional string title for the axes.                                                                                                                                                           |
| `cmap`             | Colormap (name or `Colormap` object). Default: `"viridis"`.                                                                                                                                   |
| `vmin`, `vmax`     | Scalar data range for colour normalisation.                                                                                                                                                   |
| `norm`             | A `matplotlib.colors.Normalize` instance for finer colour control.                                                                                                                            |
| `colorbar`         | `True` to add a colorbar, or a dict of kwargs forwarded to `figure.colorbar()`. Default: `False`.                                                                                             |
| `axis_labels`      | `None` (default, uses `"Longitude"` / `"Latitude"`), a dict with `"x"` / `"y"` keys, or `"none"` to suppress labels entirely.                                                                 |

## Minimal call

```python
healpix_plotting.plot(
    cell_ids,  # 1-D array of HEALPix cell IDs (uint64)
    data,  # 1-D scalar array, or (N, 3)/(N, 4) for RGB/RGBA
    healpix_grid=grid,  # HealpixGrid object (or dict)
    sampling_grid={"shape": 1024},
)
```

## Interpolation

```python
healpix_plotting.plot(..., interpolation="nearest")  # default
healpix_plotting.plot(
    ..., interpolation="bilinear"
)  # smoother, better for continuous fields
```

## Plot on an existing axis

When you pass `ax`, the `projection` parameter is ignored:

```python
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.Robinson()})
healpix_plotting.plot(..., ax=ax)
plt.show()
```

## Title and axis labels

```python
healpix_plotting.plot(
    ...,
    title="Surface temperature",
    axis_labels={"x": "Longitude", "y": "Latitude"},
)
healpix_plotting.plot(..., axis_labels="none")  # hide labels
```

## RGB / RGBA data

`plot()` accepts multi-band data directly. If `data` has shape `(N, 3)` (RGB) or `(N, 4)` (RGBA), the values are composited per-pixel rather than colour-mapped. Colourmap parameters (`cmap`, `vmin`, `vmax`, `norm`) are ignored in this mode.

```python
rgb = np.stack([r, g, b], axis=1)  # shape (N, 3)
healpix_plotting.plot(
    cell_ids, rgb, healpix_grid=healpix_grid, sampling_grid={"shape": 1024}
)
```
