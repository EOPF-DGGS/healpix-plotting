API reference
=============

Main function
-------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   healpix_plotting.plot

Low-level function
------------------

`resample` is called internally by `plot`. Use it directly only if you need the raw image array without rendering (e.g. to post-process it or display it with a different tool).

.. autosummary::
   :toctree: generated
   :nosignatures:

   healpix_plotting.resample

Classes
-------

.. autosummary::
   :toctree: generated
   :nosignatures:

   healpix_plotting.HealpixGrid
   healpix_plotting.sampling_grid.ParametrizedSamplingGrid
   healpix_plotting.sampling_grid.AffineSamplingGrid
   healpix_plotting.sampling_grid.ConcreteSamplingGrid
