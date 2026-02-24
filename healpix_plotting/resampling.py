from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy_groupies

from healpix_plotting.sampling_grid import ParametrizedSamplingGrid

if TYPE_CHECKING:
    from typing import Literal

    from healpix_plotting.healpix import HealpixGrid
    from healpix_plotting.sampling_grid import (
        ConcreteSamplingGrid,
        SamplingGrid,
        SamplingGridParameters,
    )


def is_rgb(x):
    return x.ndim == 3 and x.shape[-1] in (3, 4)


def nearest_neighbour_resampling(
    data: np.ndarray,
    target_grid: ConcreteSamplingGrid,
    source_cell_ids: np.ndarray,
    params: HealpixGrid,
    background_value: float,
) -> (np.ndarray, np.ndarray):
    target_cell_ids = params.operations.lonlat_to_healpix(
        np.reshape(target_grid.x, -1),
        np.reshape(target_grid.y, -1),
        **params.as_keyword_params(),
    )

    raw_shape = (target_cell_ids.size,)
    shape = target_grid.shape
    if is_rgb(data):
        raw_shape += data.shape[-1:]
        shape += data.shape[-1:]

    # indices that correspond to the target cells
    #
    # searchsorted(a, b) returns insert indices that insert a into b, not search
    # values from a in b. We thus need to mask all target cell ids not in the source.
    valid = np.isin(target_cell_ids, source_cell_ids)
    indices = np.searchsorted(source_cell_ids, target_cell_ids, side="left")
    valid_indices = indices[valid]

    # actual interpolation
    image = np.full(raw_shape, fill_value=background_value)
    image[valid, ...] = data[valid_indices, ...]

    return np.reshape(image, shape)


def bilinear_resampling(
    data: np.ndarray,
    target_grid: ConcreteSamplingGrid,
    source_cell_ids: np.ndarray,
    params: HealpixGrid,
) -> (np.ndarray, np.ndarray):
    raise NotImplementedError


def resample(
    cell_ids: np.ndarray,
    data: np.ndarray,
    *,
    sampling_grid: SamplingGrid | SamplingGridParameters,
    healpix_grid: HealpixGrid,
    interpolation: Literal["nearest", "bilinear"],
    agg: Literal["mean", "median", "std", "var", "min", "max", "first", "last"],
    background_value: float = np.nan,
) -> np.ndarray:
    # parameter validation
    interpolation_methods = {
        "nearest": nearest_neighbour_resampling,
        "bilinear": bilinear_resampling,
    }
    interpolator = interpolation_methods.get(interpolation)
    if interpolator is None:
        raise ValueError(f"unknown interpolation method: {interpolation}")
    # concrete sampling grid
    if isinstance(sampling_grid, dict):
        sampling_grid = ParametrizedSamplingGrid.from_dict(sampling_grid)
    target_grid = sampling_grid.resolve(cell_ids, healpix_grid)

    # deduplication
    source_cell_ids, inverse_indices = np.unique(
        cell_ids, sorted=True, return_inverse=True
    )
    deduplicated = numpy_groupies.aggregate(inverse_indices, data, func=agg)

    # interpolation
    return target_grid, interpolator(
        deduplicated,
        target_grid,
        source_cell_ids,
        healpix_grid,
        background_value=background_value,
    )
