API reference
=============

Main function
-------------

.. autosummary::
   :toctree: generated

   healpix_plot.plot

Low-level function
------------------

`resample` is called internally by `plot`. Use it directly only if you need the raw image array without rendering (e.g. to post-process it or display it with a different tool).

.. autosummary::
   :toctree: generated

   healpix_plot.resample

Classes
-------

.. autosummary::
   :toctree: generated

   healpix_plot.HealpixGrid
   healpix_plot.sampling_grid.ParametrizedSamplingGrid
   healpix_plot.sampling_grid.AffineSamplingGrid
   healpix_plot.sampling_grid.ConcreteSamplingGrid

mollview / mollgnomview
-----------------------

.. autosummary::
   :toctree: generated

   healpix_plot.mollview
   healpix_plot.mollgnomview

mollview parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 25 15 40

   * - Parameter
     - Type
     - Default
     - Description
   * - ``hpx_map``
     - ``np.ndarray``, shape ``(12·4^depth,)``
     - —
     - Input HEALPix map. RING order by default; use ``nest=True`` for NESTED.
   * - ``nest``
     - ``bool``
     - ``False``
     - Pixel ordering. ``False`` = RING (healpy default), ``True`` = NESTED.
   * - ``title``
     - ``str``
     - ``""``
     - Title displayed above the map.
   * - ``cmap``
     - ``str`` or ``Colormap``
     - ``"viridis"``
     - Matplotlib colormap. ``"RdBu_r"`` is recommended for CMB/temperature maps.
   * - ``vmin``
     - ``float`` or ``None``
     - ``None``
     - Lower bound of the colour scale. Defaults to the 2nd percentile of finite values.
   * - ``vmax``
     - ``float`` or ``None``
     - ``None``
     - Upper bound of the colour scale. Defaults to the 98th percentile of finite values.
   * - ``rot``
     - ``float``
     - ``0.0``
     - Central longitude of the map in degrees. ``rot=180`` centres on the anti-meridian.
   * - ``ellipsoid``
     - ``str``
     - ``"sphere"``
     - Reference ellipsoid for healpix-geo. ``"sphere"`` gives results identical to healpy. Other values: ``"WGS84"``, ``"GRS80"``.
   * - ``graticule``
     - ``bool``
     - ``True``
     - Draw meridians and parallels.
   * - ``graticule_step``
     - ``float``
     - ``30.0``
     - Spacing of graticule lines in degrees.
   * - ``unit``
     - ``str``
     - ``""``
     - Unit string shown below the colorbar.
   * - ``bgcolor``
     - ``str``
     - ``"black"``
     - Background colour outside the Mollweide ellipse.
   * - ``n_lon``
     - ``int``
     - ``1800``
     - Number of sample columns in the internal raster grid. Increase for high-depth maps (depth >= 8).
   * - ``n_lat``
     - ``int``
     - ``900``
     - Number of sample rows in the internal raster grid.
   * - ``norm``
     - ``Normalize`` or ``None``
     - ``None``
     - Custom matplotlib normalisation (e.g. ``LogNorm()``). Overrides ``vmin``/``vmax``.
   * - ``bad_color``
     - ``str``
     - ``"gray"``
     - Colour for ``NaN`` values.
   * - ``flip``
     - ``str``
     - ``"geo"``
     - East/west convention. ``"geo"``: east to the right (default). ``"astro"``: east to the left (astronomical convention).
   * - ``figsize``
     - ``(float, float)``
     - ``(14, 7)``
     - Figure size in inches. Only used when a new figure is created.
   * - ``colorbar``
     - ``bool``
     - ``True``
     - Show a horizontal colorbar below the map.
   * - ``hold``
     - ``bool``
     - ``False``
     - If ``True``, draw into the current axes. Ignored when ``sub`` is provided.
   * - ``sub``
     - ``(int, int, int)`` or ``None``
     - ``None``
     - ``(nrows, ncols, index)`` subplot position. Overrides ``hold``.
   * - ``coastlines``
     - ``bool``
     - ``False``
     - Overlay Natural Earth coastlines. Activates the cartopy backend — cartopy must be installed.
   * - ``coastline_kwargs``
     - ``dict`` or ``None``
     - ``None``
     - Extra kwargs forwarded to ``ax.add_feature(COASTLINE, ...)``. Only used when ``coastlines=True``.

**Returns:** ``None``. Access the current figure with ``plt.gcf()``.

**Raises:**

- ``ValueError`` — ``hpx_map.size`` is not of the form ``12·4^depth``.
- ``ValueError`` — ``flip`` is not ``"astro"`` or ``"geo"``.
- ``ImportError`` — ``coastlines=True`` but cartopy is not installed.

mollgnomview additional parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mollgnomview`` shares all parameters of ``mollview`` and adds:

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lon_center``
     - ``float``
     - —
     - Longitude of the view centre in degrees.
   * - ``lat_center``
     - ``float``
     - —
     - Latitude of the view centre in degrees.
   * - ``fov_deg``
     - ``float``
     - ``10.0``
     - Total field of view (square side) in degrees.

Internal helpers
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Function
     - Signature
     - Description
   * - ``_depth_from_npix``
     - ``(npix) → int``
     - Infers HEALPix depth. Raises ``ValueError`` for invalid sizes.
   * - ``_mollweide_inverse``
     - ``(x, y, central_lon_rad) → (lon_deg, lat_deg, valid)``
     - Analytic inverse Mollweide. Returns a boolean mask for points inside the ellipse.
   * - ``_mollweide_forward``
     - ``(lon_deg, lat_deg, central_lon_rad) → (x, y)``
     - Forward Mollweide via Newton iteration. Used for graticule drawing.
   * - ``_gnomonic_inverse``
     - ``(x, y, lon0_deg, lat0_deg) → (lon_deg, lat_deg)``
     - Inverse gnomonic projection. ``x``, ``y`` in degrees of arc.
   * - ``_rasterise``
     - ``(hpx_map, depth, lon_grid, lat_grid, valid, nest, ellipsoid) → ndarray``
     - HEALPix lookup for valid grid points. Returns ``(H, W)`` float64 with NaN outside the mask.
   * - ``_build_cmap``
     - ``(cmap, bad_color) → Colormap``
     - Builds the colormap object and sets the bad-value colour.
   * - ``_build_norm``
     - ``(data_img, vmin, vmax, norm) → Normalize``
     - Builds the normalisation object (2nd/98th percentile defaults).
   * - ``_add_colorbar``
     - ``(fig, ax, norm, cmap_obj, unit, text_color)``
     - Adds a horizontal colorbar.
   * - ``_make_axes_plain``
     - ``(hold, sub, figsize, bgcolor) → (fig, ax)``
     - Creates a plain ``Axes`` (fast path).
   * - ``_make_axes_geo``
     - ``(crs, hold, sub, figsize, bgcolor) → (fig, ax)``
     - Creates a cartopy ``GeoAxes`` (cartopy path).
   * - ``_draw_graticule_moll``
     - ``(ax, step_deg, central_lon_rad, ...)``
     - Pre-projected graticule on a plain Axes (fast path).
   * - ``_draw_graticule_on_geoaxes``
     - ``(ax, step_deg, central_lon_rad, ...)``
     - Pre-projected graticule on a GeoAxes (cartopy path, avoids ``ax.gridlines()``).
   * - ``_finalize_moll_axes``
     - ``(ax)``
     - Adds the oval boundary, clips the image to the ellipse, and sets axis limits.
   * - ``_require_cartopy``
     - ``() → (ccrs, cfeature)``
     - Lazy cartopy import with a clear error message if not installed.

Projection math
~~~~~~~~~~~~~~~

**Mollweide — Inverse formula** (Mollweide (x, y) → lon/lat, used for the fast-path grid):

The standard Mollweide ellipse spans ``x ∈ [-2√2, +2√2]``, ``y ∈ [-√2, +√2]``.
Points outside satisfy ``x²/8 + y²/2 > 1`` and are masked as NaN.

The auxiliary angle ``θ = arcsin(y / √2)``, then::

    sin(lat) = (2θ + sin(2θ)) / π
    lon      = lon_0 + π·x / (2√2·cos(θ))

**Mollweide — Forward formula** (lon/lat → Mollweide (x, y), used for graticule lines):

Requires Newton iteration to solve ``2θ + sin(2θ) = π·sin(lat)``::

    x = (2√2 / π) · (lon − lon_0) · cos(θ)
    y = √2 · sin(θ)

Convergence is reached in fewer than 10 iterations to ``tol=1e-9``.

**Gnomonic — Inverse formula** (tangent-plane → lon/lat):

The gnomonic projection is centred on ``(lon_center, lat_center)``.
Sampling coordinates ``(x, y)`` are in degrees of arc from the tangent point::

    c  = arctan(rho),    rho = sqrt(x^2 + y^2)

    lat = arcsin( cos(c)·sin(lat_0) + y·sin(c)·cos(lat_0)/rho )
    lon = lon_0 + arctan2( x·sin(c),  rho·cos(lat_0)·cos(c) − y·sin(lat_0)·sin(c) )
