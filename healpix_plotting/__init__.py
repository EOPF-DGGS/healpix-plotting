from importlib.metadata import version

try:
    __version__ = version("healpix_plotting")
except Exception:
    __version__ = "9999"

from healpix_plotting.ellipsoid import EllipsoidLike
from healpix_plotting.healpix import HealpixGrid
from healpix_plotting.plotting import plot
from healpix_plotting.resample import resample
from healpix_plotting.sampling_grid import SamplingGrid

__all__ = ["HealpixGrid", "plot", "resample", "SamplingGrid", "EllipsoidLike"]
