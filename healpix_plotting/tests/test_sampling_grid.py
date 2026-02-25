import numpy as np
import pytest
from affine import Affine

from healpix_plotting import sampling_grid as sg
from healpix_plotting.healpix import HealpixGrid


@pytest.mark.parametrize(
    ["shape", "resolution", "center", "expected_resolution", "expected_center"],
    (
        pytest.param(
            (15, 7),
            (0.5, 1.2),
            (9.5, 10.0),
            (0.5, 1.2),
            (9.5, 10.0),
        ),
        pytest.param(
            (7, 5),
            None,
            (10.0, 10.0),
            (1.875, 3.590377890729141),
            (10.0, 10.0),
        ),
        pytest.param(
            (7, 5),
            (0.5, 0.5),
            None,
            (0.5, 0.5),
            (11.25, 0.0),
        ),
    ),
)
def test_infer_parameters(
    shape, resolution, center, expected_resolution, expected_center
):
    healpix_params = HealpixGrid(level=4, indexing_scheme="nested", ellipsoid="sphere")
    cell_ids = np.array(
        [
            1079,
            1099,
            1120,
            1121,
            1122,
            1123,
            1124,
            1126,
            1127,
            1128,
            1129,
            1131,
            1132,
            1133,
            1134,
            1135,
            1144,
            1220,
        ],
        dtype="uint64",
    )
    grid = sg.ParametrizedSamplingGrid(
        shape=shape, resolution=resolution, center=center
    )
    actual_shape, actual_resolution, actual_center = sg._infer_parameters(
        grid, cell_ids, healpix_params
    )

    assert actual_shape == shape
    assert actual_resolution == pytest.approx(expected_resolution)
    assert actual_center == pytest.approx(expected_center)


class TestParametrizedSamplingGrid:
    @pytest.mark.parametrize(
        ["shape", "expected_shape"], ((3, (3, 3)), ((2, 5), (2, 5)))
    )
    @pytest.mark.parametrize(
        ["resolution", "expected_resolution"],
        ((1.1, (1.1, 1.1)), ((0.5, 1.2), (0.5, 1.2)), (None, None)),
    )
    @pytest.mark.parametrize(
        ["center", "expected_center"],
        (((10.1, 13.6), (10.1, 13.6)), ((0.5, 1.2), (0.5, 1.2)), (None, None)),
    )
    def test_from_parameters(
        self,
        shape,
        expected_shape,
        resolution,
        expected_resolution,
        center,
        expected_center,
    ):
        actual = sg.ParametrizedSamplingGrid.from_parameters(
            shape=shape, resolution=resolution, center=center
        )

        assert actual.shape == expected_shape
        assert actual.resolution == expected_resolution
        assert actual.center == expected_center

    def test_from_dict(self):
        shape = (3, 4)
        resolution = (0.5, 1.0)
        center = (45.3, -4.3)
        actual = sg.ParametrizedSamplingGrid.from_dict(
            {"shape": shape, "resolution": resolution, "center": center}
        )

        assert actual.shape == shape
        assert actual.resolution == resolution
        assert actual.center == center

    @pytest.mark.parametrize(
        ["bbox", "shape", "expected_resolution", "expected_center"],
        (
            ((0, 0, 5, 5), 10, (5 / 9, 5 / 9), (2.5, 2.5)),
            ((0, -5, 10, 15), (11, 21), (1.0, 1.0), (5, 5)),
        ),
    )
    def test_from_bbox(self, bbox, shape, expected_resolution, expected_center):
        actual = sg.ParametrizedSamplingGrid.from_bbox(bbox, shape)
        assert actual.shape == (shape if isinstance(shape, tuple) else (shape, shape))
        assert actual.resolution == expected_resolution
        assert actual.center == expected_center

    @pytest.mark.parametrize(
        ["shape", "resolution", "center", "expected"],
        (
            pytest.param(
                (5, 5),
                (1, 1),
                (0, 0),
                np.stack(
                    np.meshgrid(
                        0.0 + 1.0 * np.arange(-2, 3),
                        0.0 + 1.0 * np.arange(-2, 3),
                    ),
                    axis=0,
                ),
            ),
            pytest.param(
                (15, 7),
                (0.5, 1.2),
                (9.5, 10.0),
                np.stack(
                    np.meshgrid(
                        10.0 + 1.2 * np.arange(-3, 4),
                        9.5 + 0.5 * np.arange(-7, 8),
                    ),
                    axis=0,
                ),
            ),
        ),
    )
    def test_resolve(self, shape, resolution, center, expected):
        params = HealpixGrid(level=4, indexing_scheme="nested", ellipsoid="WGS84")
        cell_ids = np.array(
            [
                1079,
                1099,
                1120,
                1121,
                1122,
                1123,
                1124,
                1126,
                1127,
                1128,
                1129,
                1131,
                1132,
                1133,
                1134,
                1135,
                1144,
                1220,
            ],
            dtype="uint64",
        )
        grid = sg.ParametrizedSamplingGrid(
            shape=shape, resolution=resolution, center=center
        )
        actual = grid.resolve(cell_ids, params)
        expected_y, expected_x = expected
        np.testing.assert_allclose(actual.x, expected_x)
        np.testing.assert_allclose(actual.y, expected_y)


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
