import numpy as np
import pytest
from affine import Affine

from healpix_plotting import sampling_grid as sg


class TestAffineSamplingGrid:
    @pytest.mark.parametrize(
        ["transform", "shape", "expected_shape"],
        (
            (Affine.translation(1, 1), 3, (3, 3)),
            (Affine.scale(2, 0.5), (10, 15), (10, 15)),
        ),
    )
    def test_from_transform(self, transform, shape, expected_shape):
        actual = sg.AffineSamplingGrid.from_transform(transform, shape)
        assert actual.transform == transform
        assert actual.shape == expected_shape

    @pytest.mark.parametrize(
        ["transform", "shape", "expected"],
        (
            (Affine.translation(1, 1), (3, 3), np.mgrid[:3, :3].astype("float64") + 1),
            (
                Affine.scale(2, 1),
                (4, 5),
                np.mgrid[:4, :5].astype("float64") * np.array([2, 1])[:, None, None],
            ),
        ),
    )
    def test_resolve(self, transform, shape, expected):
        grid = sg.AffineSamplingGrid(transform, shape)

        # the parameters are ignored
        actual = grid.resolve(0, None)
        expected_x, expected_y = expected
        np.testing.assert_equal(actual.x, expected_x)
        np.testing.assert_equal(actual.y, expected_y)
