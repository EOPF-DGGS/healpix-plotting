"""
mollview.py
===========
HEALPix Mollweide and gnomonic visualisation backed by **cartopy**.

Replaces ~130 lines of hand-rolled projection math (_mollweide_inverse,
_mollweide_forward, _make_mollweide_grid, _draw_graticule,
_draw_graticule_labels) with cartopy's built-in CRS, which gives us:

  - Correct Mollweide / gnomonic math maintained by cartopy.
  - Automatic map boundary and clipping.
  - Gridlines and tick labels via ax.gridlines().
  - Optional coastlines / features for geoscience use-cases.

The HEALPix side (pixel-lookup) is unchanged: we still build a regular
lon/lat grid, look up the HEALPix cell for every grid point with
healpix_geo (vectorised), and hand the resulting 2-D array to
pcolormesh with transform=PlateCarree().

Dependencies
------------
    pip install healpix-geo cartopy numpy matplotlib

Public API
----------
  mollview(hpx_map, *, nest, title, cmap, vmin, vmax, rot, ellipsoid,
           graticule, graticule_step, unit, bgcolor, n_lon, n_lat,
           norm, bad_color, flip, figsize, colorbar, hold, sub,
           coastlines, coastline_kwargs)

  mollgnomview(hpx_map, lon_center, lat_center, *, fov_deg, ...)
"""

from __future__ import annotations

import math
from typing import Optional, Union

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import healpix_geo


# ---------------------------------------------------------------------------
# Depth inference  (unchanged)
# ---------------------------------------------------------------------------

def _depth_from_npix(npix: int) -> int:
    """Infer HEALPix depth from map size (npix = 12 · 4**depth)."""
    if npix <= 0 or npix % 12 != 0:
        raise ValueError(
            f"Invalid map size: {npix}. "
            "npix must be a multiple of 12 of the form 12 * 4**depth."
        )
    nside_sq   = npix // 12
    log2_nside = math.log2(nside_sq) / 2.0
    depth      = int(round(log2_nside))
    if 12 * (4**depth) != npix:
        raise ValueError(
            f"Invalid map size: {npix}. "
            f"Closest valid size is depth={depth} → npix={12 * 4**depth}."
        )
    return depth


# ---------------------------------------------------------------------------
# Shared rasterisation kernel
# ---------------------------------------------------------------------------

def _rasterise(
    hpx_map:   np.ndarray,
    depth:     int,
    lon_grid:  np.ndarray,   # (H, W)  degrees, PlateCarree
    lat_grid:  np.ndarray,   # (H, W)  degrees
    nest:      bool,
    ellipsoid: str,
) -> np.ndarray:
    """
    Look up the HEALPix value for every (lon, lat) grid point.

    Returns a (H, W) float64 array; NaN where the map value is NaN.
    """
    H, W = lon_grid.shape

    if nest:
        ipix = healpix_geo.nested.lonlat_to_healpix(
            lon_grid.ravel(), lat_grid.ravel(), depth, ellipsoid=ellipsoid,
        )
    else:
        ipix = healpix_geo.ring.lonlat_to_healpix(
            lon_grid.ravel(), lat_grid.ravel(), depth, ellipsoid=ellipsoid,
        )

    return hpx_map[ipix].reshape(H, W)


# ---------------------------------------------------------------------------
# Axes creation helper
# ---------------------------------------------------------------------------

def _make_axes(
    crs,
    hold:    bool,
    sub:     Optional[tuple[int, int, int]],
    figsize: tuple[float, float],
    bgcolor: str,
):
    """
    Return (fig, ax) for a cartopy GeoAxes.

    Priority: sub > hold=True > hold=False (new figure).

    Notes
    -----
    When hold=True the *current* axes is reused.  The caller is responsible
    for having created it with the matching CRS — cartopy axes cannot change
    their projection after creation.  If plt.gca() is not a GeoAxes an
    informative error will be raised by cartopy at draw time.
    """
    if sub is not None:
        nrows, ncols, idx = sub
        fig = plt.gcf()
        ax  = fig.add_subplot(nrows, ncols, idx, projection=crs,
                              facecolor=bgcolor)
    elif hold:
        ax  = plt.gca()
        fig = ax.get_figure()
    else:
        fig = plt.figure(figsize=figsize, facecolor=bgcolor)
        ax  = fig.add_subplot(1, 1, 1, projection=crs, facecolor=bgcolor)

    fig.patch.set_facecolor(bgcolor)
    return fig, ax


# ---------------------------------------------------------------------------
# mollview
# ---------------------------------------------------------------------------

def mollview(
    hpx_map:          np.ndarray,
    *,
    nest:             bool                      = False,
    title:            str                       = "",
    cmap:             Union[str, mcolors.Colormap] = "viridis",
    vmin:             Optional[float]           = None,
    vmax:             Optional[float]           = None,
    rot:              float                     = 0.0,
    ellipsoid:        str                       = "sphere",
    graticule:        bool                      = True,
    graticule_step:   float                     = 30.0,
    unit:             str                       = "",
    bgcolor:          str                       = "black",
    # Resolution of the raster grid (replaces width_px / height_px).
    # n_lon × n_lat points are sampled in PlateCarree space and reprojected
    # by cartopy — independent of the final figure DPI.
    n_lon:            int                       = 1800,
    n_lat:            int                       = 900,
    norm:             Optional[mcolors.Normalize] = None,
    bad_color:        str                       = "gray",
    flip:             str                       = "geo",
    figsize:          tuple[float, float]       = (14, 7),
    colorbar:         bool                      = True,
    hold:             bool                      = False,
    sub:              Optional[tuple[int,int,int]] = None,
    # Extras made possible by cartopy
    coastlines:       bool                      = False,
    coastline_kwargs: Optional[dict]            = None,
) -> None:
    """
    Display a HEALPix map in the Mollweide equal-area projection.

    Uses cartopy's ``Mollweide`` CRS for all projection math.  No
    custom Mollweide formulae are needed in this module.

    Parameters
    ----------
    hpx_map : np.ndarray, shape (12 · 4**depth,)
        Input HEALPix map.  RING order by default; pass ``nest=True``
        for NESTED maps.
    nest : bool, default False
        Pixel ordering of the input map.
    title : str
        Figure title.
    cmap : str or Colormap, default "viridis"
        Matplotlib colormap.
    vmin, vmax : float or None
        Colour-scale limits.  Default: 2nd / 98th percentile.
    rot : float, default 0.0
        Central longitude in degrees.  Forwarded to
        ``cartopy.crs.Mollweide(central_longitude=rot)``.
    ellipsoid : str, default "sphere"
        Reference ellipsoid for healpix_geo.  Use ``"WGS84"`` for
        geographic data.
    graticule : bool, default True
        Draw meridians and parallels via ``ax.gridlines()``.
    graticule_step : float, default 30.0
        Graticule spacing in degrees.
    unit : str
        Colorbar label.
    bgcolor : str, default "black"
        Figure and axes background colour.
    n_lon, n_lat : int, default 1800, 900
        Resolution of the internal lon/lat sampling grid.  Higher values
        produce sharper maps at the cost of more memory and CPU.
    norm : Normalize or None
        Custom matplotlib normalisation (e.g. ``LogNorm()``).
    bad_color : str, default "gray"
        Colour for NaN pixels.
    flip : str, default "geo"
        ``"astro"`` — east to the left (astronomical convention).
        ``"geo"``   — east to the right (geographic convention, cartopy default).
    figsize : (float, float), default (14, 7)
        Figure size in inches.  Used only when a new figure is created.
    colorbar : bool, default True
        Show a horizontal colorbar.
    hold : bool, default False
        If True, draw into the current axes (must already be a cartopy
        GeoAxes with Mollweide projection).
    sub : (nrows, ncols, idx) or None
        Subplot position within the current figure.
    coastlines : bool, default False
        Overlay Natural Earth coastlines (useful for geographic maps).
    coastline_kwargs : dict or None
        Extra kwargs forwarded to ``ax.add_feature(COASTLINE, ...)``.

    Returns
    -------
    None  — same behaviour as healpy.mollview.
    """
    if flip not in ("astro", "geo"):
        raise ValueError(f"flip must be 'astro' or 'geo', got {flip!r}.")

    hpx_map = np.asarray(hpx_map, dtype=np.float64)
    depth   = _depth_from_npix(hpx_map.size)

    # ------------------------------------------------------------------
    # 1. Regular lon/lat sampling grid in PlateCarree space
    #    (replaces _make_mollweide_grid + _mollweide_inverse)
    #
    #    For the "astro" flip we reverse the longitude direction so that
    #    east appears on the left after cartopy reprojects the data.
    # ------------------------------------------------------------------
    lons_1d = np.linspace(-180.0, 180.0, n_lon, endpoint=False)
    lats_1d = np.linspace(-90.0, 90.0, n_lat)

    # Always keep lons monotonically increasing so cartopy doesn't misread
    # the extent.  The astro flip (east-left) is applied to the data array
    # AFTER the HEALPix lookup, not to the coordinate axis.
    lon_grid, lat_grid = np.meshgrid(lons_1d, lats_1d)   # (n_lat, n_lon)

    # ------------------------------------------------------------------
    # 2. HEALPix lookup — single vectorised call per ordering
    # ------------------------------------------------------------------
    data_img = _rasterise(hpx_map, depth, lon_grid, lat_grid, nest, ellipsoid)

    # Astro flip: mirror columns so east is to the left.
    # We also need to reverse lons_1d so the imshow extent is consistent.
    if flip == "astro":
        data_img = data_img[:, ::-1]
        lons_1d  = -lons_1d[::-1].copy()

    # ------------------------------------------------------------------
    # 3. Colour normalisation
    # ------------------------------------------------------------------
    valid = data_img[np.isfinite(data_img)]
    if norm is None:
        if vmin is None:
            vmin = float(np.percentile(valid, 2))  if valid.size else 0.0
        if vmax is None:
            vmax = float(np.percentile(valid, 98)) if valid.size else 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    cmap_obj = plt.get_cmap(cmap).copy() if isinstance(cmap, str) \
               else (cmap.copy() if hasattr(cmap, "copy") else cmap)
    cmap_obj.set_bad(bad_color, alpha=1.0)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)

    # ------------------------------------------------------------------
    # 4. Axes  — cartopy Mollweide CRS handles all projection math
    # ------------------------------------------------------------------
    crs      = ccrs.Mollweide(central_longitude=rot)
    plate_cr = ccrs.PlateCarree()

    fig, ax = _make_axes(crs, hold, sub, figsize, bgcolor)
    ax.set_global()

    # ------------------------------------------------------------------
    # 5. Plot — imshow is orders of magnitude faster than pcolormesh here.
    #
    #    pcolormesh with transform=PlateCarree() forces cartopy to reproject
    #    every single quad cell (~1.6M for the default grid), making it
    #    extremely slow.
    #
    #    imshow with transform=PlateCarree() and a global extent lets cartopy
    #    treat the data as a single image warp — one reprojection operation
    #    for the whole array, regardless of n_lon × n_lat.
    #
    #    The extent must be [lon_min, lon_max, lat_min, lat_max] in the
    #    PlateCarree frame.  We flip lat so that row 0 = south (origin="lower").
    # ------------------------------------------------------------------
    ax.imshow(
        data_img,
        origin="lower",
        extent=[lons_1d[0], lons_1d[-1], lats_1d[0], lats_1d[-1]],
        transform=plate_cr,
        cmap=cmap_obj,
        norm=norm,
        interpolation="nearest",
        aspect="auto",
    )

    # ------------------------------------------------------------------
    # 6. Optional extras — all handled by cartopy, zero custom math
    # ------------------------------------------------------------------
    if coastlines:
        ckw = dict(linewidth=0.6, edgecolor="white") if coastline_kwargs is None \
              else coastline_kwargs
        ax.add_feature(cfeature.COASTLINE, **ckw)

    if graticule:
        gl = ax.gridlines(
            crs=plate_cr,
            draw_labels=False,
            linewidth=0.5,
            color="white",
            linestyle="--",
            alpha=0.7,
            xlocs=np.arange(-180, 181, graticule_step),
            ylocs=np.arange(-90,   91, graticule_step),
        )
        # Equator slightly thicker
        ax.gridlines(
            crs=plate_cr,
            draw_labels=False,
            linewidth=0.8,
            color="white",
            linestyle="-",
            alpha=0.7,
            xlocs=[],
            ylocs=[0],
        )

    if title:
        ax.set_title(title, color="white", fontsize=11, pad=6)

    # ------------------------------------------------------------------
    # 7. Colorbar
    # ------------------------------------------------------------------
    if colorbar:
        cbar = fig.colorbar(
            sm, ax=ax,
            orientation="horizontal",
            fraction=0.025, pad=0.04, aspect=40, shrink=0.6,
        )
        cbar.ax.tick_params(colors="white", labelsize=8)
        cbar.outline.set_edgecolor("white")
        if unit:
            cbar.set_label(unit, color="white", fontsize=9)

    # No return — same behaviour as healpy.mollview


# ---------------------------------------------------------------------------
# mollgnomview  — gnomonic zoom
# ---------------------------------------------------------------------------

def mollgnomview(
    hpx_map:    np.ndarray,
    lon_center: float,
    lat_center: float,
    *,
    nest:       bool                         = False,
    fov_deg:    float                        = 10.0,
    title:      str                          = "",
    cmap:       Union[str, mcolors.Colormap] = "viridis",
    vmin:       Optional[float]              = None,
    vmax:       Optional[float]              = None,
    ellipsoid:  str                          = "sphere",
    unit:       str                          = "",
    n_lon:      int                          = 800,
    n_lat:      int                          = 800,
    figsize:    tuple[float, float]          = (7, 7),
    colorbar:   bool                         = True,
    hold:       bool                         = False,
    sub:        Optional[tuple[int,int,int]] = None,
    coastlines: bool                         = False,
    coastline_kwargs: Optional[dict]         = None,
) -> None:
    """
    Local zoom in the gnomonic (tangent-plane) projection.

    Uses ``cartopy.crs.Gnomonic`` — replaces the hand-rolled deprojection
    that was in the original ``mollgnomview``.

    Parameters
    ----------
    hpx_map : np.ndarray
        HEALPix map (RING by default; NESTED if ``nest=True``).
    lon_center, lat_center : float
        Centre of the view in degrees.
    fov_deg : float, default 10.0
        Total field of view (square side) in degrees.
    n_lon, n_lat : int, default 800
        Resolution of the internal sampling grid.
    All other parameters : see ``mollview``.
    """
    hpx_map = np.asarray(hpx_map, dtype=np.float64)
    depth   = _depth_from_npix(hpx_map.size)

    # Build a lon/lat grid that covers the requested FOV.
    # We use PlateCarree for the grid and let cartopy reproject to Gnomonic.
    half    = fov_deg / 2.0
    lons_1d = np.linspace(lon_center - half, lon_center + half, n_lon)
    lats_1d = np.linspace(lat_center - half, lat_center + half, n_lat)
    lon_grid, lat_grid = np.meshgrid(lons_1d, lats_1d)

    data_img = _rasterise(hpx_map, depth, lon_grid, lat_grid, nest, ellipsoid)

    valid = data_img[np.isfinite(data_img)]
    if vmin is None:
        vmin = float(np.percentile(valid, 2))  if valid.size else 0.0
    if vmax is None:
        vmax = float(np.percentile(valid, 98)) if valid.size else 1.0
    norm_obj = mcolors.Normalize(vmin=vmin, vmax=vmax)

    cmap_obj = plt.get_cmap(cmap).copy() if isinstance(cmap, str) \
               else (cmap.copy() if hasattr(cmap, "copy") else cmap)
    cmap_obj.set_bad("gray", alpha=1.0)

    crs      = ccrs.Gnomonic(central_latitude=lat_center,
                              central_longitude=lon_center)
    plate_cr = ccrs.PlateCarree()

    fig, ax  = _make_axes(crs, hold, sub, figsize, "black")
    ax.set_extent(
        [lon_center - half, lon_center + half,
         lat_center - half, lat_center + half],
        crs=plate_cr,
    )

    im = ax.imshow(
        data_img,
        origin="lower",
        extent=[lons_1d[0], lons_1d[-1], lats_1d[0], lats_1d[-1]],
        transform=plate_cr,
        cmap=cmap_obj,
        norm=norm_obj,
        interpolation="nearest",
        aspect="auto",
    )

    if coastlines:
        ckw = dict(linewidth=0.6, edgecolor="white") if coastline_kwargs is None \
              else coastline_kwargs
        ax.add_feature(cfeature.COASTLINE, **ckw)

    ax.gridlines(crs=plate_cr, draw_labels=True,
                 linewidth=0.5, color="gray", linestyle="--", alpha=0.6)

    if title:
        ax.set_title(title, color="white", fontsize=11)

    if colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(colors="white")
        cbar.outline.set_edgecolor("white")
        if unit:
            cbar.set_label(unit, color="white")
