from __future__ import annotations

from typing import TYPE_CHECKING

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from healpix_plotting.healpix import HealpixGrid
from healpix_plotting.resampling import resample
from healpix_plotting.sampling_grid import ParametrizedSamplingGrid

if TYPE_CHECKING:
    from typing import Any, Literal

    import cartopy.crs as ccrs
    from matplotlib.axis import Axis
    from matplotlib.cm import ColorMap
    from matplotlib.norm import Norm

    from healpix_plotting.sampling_grid import SamplingGrid, SamplingGridParameters


def plot(
    cell_ids: np.ndarray,
    data: np.ndarray,
    *,
    healpix_grid: HealpixGrid,
    sampling_grid: SamplingGridParameters | SamplingGrid,
    projection: str | ccrs.CRS = "Mollweide",
    agg: str = "mean",
    interpolation: str = "nearest",
    background_value: float = np.nan,
    rgb_clip: tuple[float, float] = (0.0, 1.0),
    ax: Axis | None = None,
    title: str | None = None,
    colorbar: bool | dict[str, Any] = False,
    cmap: str | ColorMap = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Norm | None = None,
    axis_labels: dict[str, str] | Literal["none"] | None = None,
) -> Axis:
    """resample and plot healpix data

    Parameters
    ----------
    cell_ids : numpy.ndarray
        The cell ids describing the spatial position of the data.
    data : numpy.ndarray
        The data to plot. If 1D, will be color-coded using the standard
        matplotlib mechanisms. If 2D, the last axis must have a size of 3 (for
        RGB) or 4 (for RGBA).
    healpix_grid : HealpixGrid or dict of str to any
        The healpix grid parameters necessary to interpret ``cell_ids``.
    sampling_grid : SamplingGrid or dict of str to any
        The target grid.
    projection : str or cartopy.crs.CRS
        The projection used to construct a new axis. Ignored if ``ax`` is given.
    agg : str, default: "mean"
        Aggregation to deduplicate the data.
    interpolation : str, default: "nearest"
        The algorithm used to interpolate from healpix to the target grid. Available values:

        - ``"nearest"``: nearest-neighbour resampling
        - ``"bilinear"``: bilinear resampling

    background_value : float, default: numpy.nan
        The background value for missing values.
    ax : matplotlib.axis.Axis, optional
        The axis to plot on. If not passed, a new figure with a single axis is
        created using ``projection`` and ``figure_params``.
    vmin : float, optional
        Minimum value to color-code.
    vmax : float, optional
        Maximum value to color-code.
    norm : matplotlib.norm.Norm, optional
        Normalization class for more control.
    cmap : str or matplotlib.colors.Colormap, default: "viridis"
        The colormap to use for plotting.
    axis_labels : dict of str to str or "none", optional
        Axis labels. Possible values:

        - if ``None`` or not passed, ``"Longitude"`` and ``"Latitude"`` are used.
        - dict: the keys ``"x"`` and ``"y"`` are used
        - ``"none"``: no axis labels

    Returns
    -------
    mappable : matplotlib.image.AxisImage
        The mappable of the image to allow further processing.

    Examples
    --------
    >>> import healpix_plotting
    >>> import numpy as np

    Define the source grid:

    >>> healpix_params = healpix_plotting.HealpixParameters(
    ...     level=4,
    ...     indexing_scheme="nested",
    ... )
    >>> cell_ids = np.arange(12 * 4 ** healpix_params["level"], dtype="uint64")

    Create the data:

    >>> lon, lat = healpix_params.operations.healpix_to_lonlat(
    ...     cell_ids,
    ...     **healpix_params.as_keyword_params(),
    ... )
    >>> data = np.cos(8 * np.deg2rad(lon)) * np.sin(4 * np.deg2rad(lat))

    Plot the data

    >>> healpix_plotting.plot(
    ...     cell_ids,
    ...     data,
    ...     sampling_grid={"shape": 1024},
    ...     healpix_grid=healpix_params,
    ... )  # doctest: +ELLIPSIS
    <matplotlib.image.AxesImage at 0x...>
    """
    if isinstance(sampling_grid, dict):
        sampling_grid = ParametrizedSamplingGrid.from_dict(sampling_grid)
    if isinstance(healpix_grid, dict):
        healpix_grid = HealpixGrid(**healpix_grid)

    target_grid, image = resample(
        cell_ids,
        data,
        sampling_grid=sampling_grid,
        healpix_grid=healpix_grid,
        interpolation=interpolation,
        agg=agg,
        background_value=background_value,
    )
    if isinstance(projection, str):
        _projection = getattr(ccrs, projection, None)
        if _projection is None:
            raise ValueError(f"unknown projection: {projection}")
        projection = _projection()

    if ax is None:
        fig, ax = plt.subplots(
            figsize=(12, 10),
            subplot_kw={"projection": projection},
            layout="constrained",
        )

    mappable = ax.imshow(
        image,
        extent=target_grid.extent,
        origin="lower",
        interpolation="nearest",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
        norm=norm,
        cmap=cmap,
        transform=ccrs.PlateCarree(),
    )
    ax.set_extent(target_grid.extent, crs=ccrs.PlateCarree())
    if title is not None:
        ax.set_title(title)

    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    if colorbar:
        colorbar_kwargs = colorbar if isinstance(colorbar, dict) else {}
        ax.figure.colorbar(mappable, **colorbar_kwargs)

    if axis_labels != "none":
        if axis_labels is None:
            axis_labels = {"x": "Longitude", "y": "Latitude"}

        for axis in ["x", "y"]:
            getattr(ax, f"set_{axis}label")(axis_labels[axis])

    return mappable
