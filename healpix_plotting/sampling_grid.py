from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, TypedDict

import numpy as np

if TYPE_CHECKING:
    from typing import Self

    import affine

    from healpix_plotting.healpix import HealpixGrid


class SamplingGridParameters(TypedDict):
    """Sampling parameters as a dict

    Parameters
    ----------
    shape : int or tuple of int, default: 1024
        The shape of the array. If a int, the shape is a square of equal size.
    resolution : float or tuple of float, optional
        The resolution or step size of the sampling grid. If a float, expands to
        a 2-tuple of equal values. If missing, derived from the spatial extent
        of the data and the given shape.
    center : tuple of float, optional
        The center of the sampling grid. If missing, this is inferred from the data.
    """

    shape: int | tuple[int, int] = 1024
    resolution: float | tuple[float, float] | None = None
    center: tuple[float, float] | None = None


class SamplingGrid:
    def resolve(
        self, cell_ids: np.ndarray, parameters: HealpixGrid
    ) -> ConcreteSamplingGrid:  # pragma: no cover
        raise NotImplementedError


def _infer_parameters(
    grid: ParametrizedSamplingGrid, cell_ids: np.ndarray, params: HealpixGrid
) -> (tuple[int, int], tuple[float, float], tuple[float, float]):
    # TODO: actually resolve the parameters
    center = grid.center
    shape = grid.shape
    resolution = grid.resolution
    if resolution is None or center is None:
        lon, lat = params.operations.healpix_to_lonlat(
            cell_ids, **params.as_keyword_params()
        )
        if center is None:
            center = (np.mean(lon).item(), np.mean(lat).item())

        if resolution is None:
            size_x, size_y = shape
            min_x, max_x = np.min(lon).item(), np.max(lon).item()
            min_y, max_y = np.min(lat).item(), np.max(lat).item()

            dx = (max_x - min_x) / (size_x - 1)
            dy = (max_y - min_y) / (size_y - 1)

            resolution = (dx, dy)

    return shape, resolution, center


@dataclass
class ParametrizedSamplingGrid:
    shape: tuple[int, int]
    resolution: tuple[float, float] | None
    center: tuple[float, float] | None

    @classmethod
    def from_parameters(
        cls,
        shape: int | tuple[int, int],
        resolution: float | tuple[float, float] | None = None,
        center: tuple[float, float] | None = None,
    ) -> Self:
        if isinstance(shape, int):
            shape = (shape, shape)
        if isinstance(resolution, float):
            resolution = (resolution, resolution)

        return cls(shape=shape, resolution=resolution, center=center)

    @classmethod
    def from_dict(cls, mapping: SamplingGridParameters) -> Self:
        return cls.from_parameters(**mapping)

    @classmethod
    def from_bbox(
        cls,
        bbox: tuple[float, float, float, float],
        shape: int | tuple[int, int],
    ) -> Self:
        if isinstance(shape, int):
            shape = (shape, shape)

        xmin, ymin, xmax, ymax = bbox

        center = (float(np.mean([xmin, xmax])), float(np.mean([ymin, ymax])))

        resolution = (
            (xmax - xmin) / (shape[0] - 1),
            (ymax - ymin) / (shape[1] - 1),
        )

        return cls(shape=shape, center=center, resolution=resolution)

    def resolve(
        self, cell_ids: np.ndarray, parameters: HealpixGrid
    ) -> ConcreteSamplingGrid:
        shape, resolution, center = _infer_parameters(self, cell_ids, parameters)

        size_x, size_y = shape
        half_x = resolution[0] * size_x / 2
        half_y = resolution[1] * size_y / 2
        center_x, center_y = center

        xs = (
            np.linspace(-half_x, +half_x, size_x, endpoint=False)
            + (half_x / size_x)
            + center_x
        )
        ys = (
            np.linspace(-half_y, +half_y, size_y, endpoint=False)
            + (half_y / size_y)
            + center_y
        )

        x, y = np.meshgrid(xs, ys)

        return ConcreteSamplingGrid(x, y)


@dataclass
class AffineSamplingGrid(SamplingGrid):
    transform: affine.Affine
    shape: tuple[int, int]

    @classmethod
    def from_transform(
        cls,
        transform: affine.Affine,
        shape: int | tuple[int, int],
    ) -> Self:
        if isinstance(shape, int):
            shape = (shape, shape)

        return cls(transform, shape)

    def resolve(
        self, cell_ids: np.ndarray, parameters: HealpixGrid
    ) -> ConcreteSamplingGrid:
        pixel_x, pixel_y = np.mgrid[: self.shape[0], : self.shape[1]]

        x, y = self.transform * (pixel_x, pixel_y)

        return ConcreteSamplingGrid(x, y)


@dataclass
class ConcreteSamplingGrid:
    x: np.ndarray
    y: np.ndarray

    @property
    def shape(self):
        return self.x.shape
