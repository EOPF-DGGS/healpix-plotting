"""
mollview.py  (hybrid: fast by default, cartopy on demand)
==========================================================
HEALPix Mollweide and gnomonic visualisation.

Backend selection
-----------------
The rendering backend is chosen **automatically** based on whether
coastlines are requested:

  coastlines=False  (default)
      Pure-matplotlib fast path.  cartopy is NOT imported at all.
      All projection math is done in numpy (Mollweide / gnomonic inverse
      formulas).  The 2-D image is built in the output coordinate system,
      so ax.imshow() needs no transform and no render-time reprojection.
      Typical time (nside=64, 1800x900): ~70 ms

  coastlines=True
      Cartopy path.  A GeoAxes is created, HEALPix data is plotted via
      imshow(transform=PlateCarree()), and Natural Earth coastlines are
      overlaid with ax.add_feature().
      The graticule still uses the fast pre-projected approach to avoid
      calling ax.gridlines() (~400 ms penalty).
      Typical time (nside=64, 1800x900): ~500 ms
      (vs ~1.1 s for the original all-cartopy code)

Why cartopy is slow
-------------------
Three sources of latency, all unrelated to HEALPix:

  1. GeoAxes init              ~200-400 ms  (CRS pipeline, clip path, bbox)
  2. ax.gridlines()            ~300-600 ms  (clips each line to the ellipse)
  3. imshow(transform=...)     ~200-400 ms  (full-image resample at draw time)

The fast path eliminates all three.
The cartopy path eliminates (2) by reusing the pre-projected graticule.

Public API
----------
  mollview(hpx_map, *, nest, title, cmap, vmin, vmax, rot, ellipsoid,
           graticule, graticule_step, unit, bgcolor, n_lon, n_lat,
           norm, bad_color, flip, figsize, colorbar, hold, sub,
           coastlines, coastline_kwargs)

  mollgnomview(hpx_map, lon_center, lat_center, *, fov_deg, ...)

Dependencies
------------
  Always required   : healpix-geo, numpy, matplotlib
  Only when needed  : cartopy  (only when coastlines=True)
"""

from __future__ import annotations

import math

import healpix_geo
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

# ---------------------------------------------------------------------------
# Depth inference
# ---------------------------------------------------------------------------


def _depth_from_npix(npix: int) -> int:
    """Infer HEALPix depth from map size (npix = 12 * 4**depth)."""
    if npix <= 0 or npix % 12 != 0:
        raise ValueError(f"Invalid map size: {npix}.")
    nside_sq = npix // 12
    log2_nside = math.log2(nside_sq) / 2.0
    depth = int(round(log2_nside))
    if 12 * (4**depth) != npix:
        raise ValueError(f"Invalid map size: {npix}. Closest: depth={depth}.")
    return depth


# ---------------------------------------------------------------------------
# Mollweide projection math
# ---------------------------------------------------------------------------

_SQRT2 = math.sqrt(2.0)
_2SQRT2 = 2.0 * _SQRT2
_INV_PI = 1.0 / math.pi


def _mollweide_inverse(x, y, central_lon_rad=0.0):
    """
    Inverse Mollweide: (x, y) -> (lon_deg, lat_deg, valid_mask).
    Ellipse spans x in [-2*sqrt2, +2*sqrt2], y in [-sqrt2, +sqrt2].
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    valid = (x * x / 8.0 + y * y / 2.0) <= 1.0
    sin_theta = np.clip(np.where(valid, y / _SQRT2, 0.0), -1.0, 1.0)
    theta = np.arcsin(sin_theta)
    lat_rad = np.arcsin(
        np.clip((2.0 * theta + np.sin(2.0 * theta)) * _INV_PI, -1.0, 1.0)
    )
    cos_theta = np.cos(theta)
    safe_cos = np.where(np.abs(cos_theta) < 1e-12, 1.0, cos_theta)
    lon_rad = central_lon_rad + (math.pi * x) / (_2SQRT2 * safe_cos)
    lon_deg = (np.degrees(lon_rad) + 180.0) % 360.0 - 180.0
    lat_deg = np.degrees(lat_rad)
    return lon_deg, lat_deg, valid


def _mollweide_forward(lon_deg, lat_deg, central_lon_rad=0.0, max_iter=10):
    """
    Forward Mollweide: (lon_deg, lat_deg) -> (x, y).
    Newton iteration for  2*theta + sin(2*theta) = pi*sin(lat).
    Used only for graticule drawing, not in the hot path.
    """
    lon_rad = np.radians(lon_deg)
    lat_rad = np.radians(lat_deg)
    rhs = math.pi * np.sin(lat_rad)
    theta = lat_rad.copy()
    for _ in range(max_iter):
        denom = 2.0 * (1.0 + np.cos(2.0 * theta))
        denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
        delta = -(2.0 * theta + np.sin(2.0 * theta) - rhs) / denom
        theta += delta
        if np.max(np.abs(delta)) < 1e-9:
            break
    theta = np.clip(theta, -math.pi / 2, math.pi / 2)
    dlon = (lon_rad - central_lon_rad + math.pi) % (2 * math.pi) - math.pi
    return _2SQRT2 * _INV_PI * dlon * np.cos(theta), _SQRT2 * np.sin(theta)


# ---------------------------------------------------------------------------
# Gnomonic projection math
# ---------------------------------------------------------------------------


def _gnomonic_inverse(x, y, lon0_deg, lat0_deg):
    """Inverse gnomonic: tangent-plane (x, y) in degrees -> (lon_deg, lat_deg)."""
    x = np.radians(np.asarray(x, dtype=np.float64))
    y = np.radians(np.asarray(y, dtype=np.float64))
    lat0 = math.radians(lat0_deg)
    lon0 = math.radians(lon0_deg)
    rho = np.sqrt(x * x + y * y)
    c = np.arctan(rho)
    cos_c = np.cos(c)
    sin_c = np.sin(c)
    cos_lat0 = math.cos(lat0)
    sin_lat0 = math.sin(lat0)
    safe_rho = np.where(rho < 1e-12, 1.0, rho)
    lat = np.arcsin(cos_c * sin_lat0 + y * sin_c * cos_lat0 / safe_rho)
    lon = lon0 + np.arctan2(
        x * sin_c, safe_rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c
    )
    return np.degrees(lon), np.degrees(lat)


# ---------------------------------------------------------------------------
# Shared: HEALPix rasterisation
# ---------------------------------------------------------------------------


def _rasterise(hpx_map, depth, lon_grid, lat_grid, valid, nest, ellipsoid):
    """HEALPix lookup. Returns (H, W) float64; NaN outside valid mask."""
    H, W = lon_grid.shape
    result = np.full(H * W, np.nan)
    if valid.any():
        flat = valid.ravel()
        fn = (healpix_geo.nested if nest else healpix_geo.ring).lonlat_to_healpix
        ipix = fn(
            lon_grid.ravel()[flat], lat_grid.ravel()[flat], depth, ellipsoid=ellipsoid
        )
        result[flat] = hpx_map[ipix]
    return result.reshape(H, W)


# ---------------------------------------------------------------------------
# Shared: colour helpers
# ---------------------------------------------------------------------------


def _build_cmap(cmap, bad_color):
    obj = (
        plt.get_cmap(cmap).copy()
        if isinstance(cmap, str)
        else (cmap.copy() if hasattr(cmap, "copy") else cmap)
    )
    obj.set_bad(bad_color, alpha=1.0)
    return obj


def _build_norm(data_img, vmin, vmax, norm):
    if norm is not None:
        return norm
    finite = data_img[np.isfinite(data_img)]
    lo = float(np.percentile(finite, 2)) if finite.size else 0.0
    hi = float(np.percentile(finite, 98)) if finite.size else 1.0
    return mcolors.Normalize(
        vmin=vmin if vmin is not None else lo, vmax=vmax if vmax is not None else hi
    )


def _add_colorbar(fig, ax, norm, cmap_obj, unit, text_color="white"):
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(
        sm,
        ax=ax,
        orientation="horizontal",
        fraction=0.025,
        pad=0.04,
        aspect=40,
        shrink=0.6,
    )
    cbar.ax.tick_params(colors=text_color, labelsize=8)
    cbar.outline.set_edgecolor(text_color)
    if unit:
        cbar.set_label(unit, color=text_color, fontsize=9)


# ---------------------------------------------------------------------------
# Fast-path helpers  (pure matplotlib, no cartopy)
# ---------------------------------------------------------------------------


def _make_axes_plain(hold, sub, figsize, bgcolor):
    if sub is not None:
        nrows, ncols, idx = sub
        fig = plt.gcf()
        ax = fig.add_subplot(nrows, ncols, idx, facecolor=bgcolor)
    elif hold:
        ax = plt.gca()
        fig = ax.get_figure()
    else:
        fig = plt.figure(figsize=figsize, facecolor=bgcolor)
        ax = fig.add_subplot(1, 1, 1, facecolor=bgcolor)
    fig.patch.set_facecolor(bgcolor)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def _draw_graticule_moll(
    ax, step_deg, central_lon_rad, lw=0.5, color="white", ls="--", alpha=0.7, n_pts=360
):
    """Pre-projected graticule on a plain Axes -- no CRS pipeline."""
    step = float(step_deg)
    lons_r = np.linspace(-180.0, 180.0, n_pts)
    lats_r = np.linspace(-90.0, 90.0, n_pts)
    for lat in np.arange(-90.0, 90.0 + 1e-9, step):
        xp, yp = _mollweide_forward(lons_r, np.full_like(lons_r, lat), central_lon_rad)
        ax.plot(
            xp,
            yp,
            color=color,
            lw=lw * 1.6 if lat == 0.0 else lw,
            ls="-" if lat == 0.0 else ls,
            alpha=alpha,
            zorder=2,
        )
    for lon in np.arange(-180.0, 180.0 + 1e-9, step):
        xp, yp = _mollweide_forward(np.full_like(lats_r, lon), lats_r, central_lon_rad)
        ax.plot(xp, yp, color=color, lw=lw, ls=ls, alpha=alpha, zorder=2)


def _finalize_moll_axes(ax):
    """Add oval boundary, clip image to ellipse, set axis limits."""
    ell_border = Ellipse(
        xy=(0, 0),
        width=2 * _2SQRT2,
        height=2 * _SQRT2,
        facecolor="none",
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
        transform=ax.transData,
    )
    ax.add_patch(ell_border)
    clip_ell = Ellipse(
        xy=(0, 0), width=2 * _2SQRT2, height=2 * _SQRT2, transform=ax.transData
    )
    for img in ax.get_images():
        img.set_clip_path(clip_ell)
    ax.set_xlim(-_2SQRT2 * 1.02, _2SQRT2 * 1.02)
    ax.set_ylim(-_SQRT2 * 1.08, _SQRT2 * 1.08)


# ---------------------------------------------------------------------------
# Cartopy-path helpers  (imported lazily)
# ---------------------------------------------------------------------------


def _require_cartopy():
    """Import cartopy lazily; raise a clear error if not installed."""
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        return ccrs, cfeature
    except ImportError as exc:
        raise ImportError(
            "cartopy is required when coastlines=True.\n"
            "Install it with:  pip install cartopy\n"
            f"Original error: {exc}"
        ) from exc


def _make_axes_geo(crs, hold, sub, figsize, bgcolor):
    if sub is not None:
        nrows, ncols, idx = sub
        fig = plt.gcf()
        ax = fig.add_subplot(nrows, ncols, idx, projection=crs, facecolor=bgcolor)
    elif hold:
        ax = plt.gca()
        fig = ax.get_figure()
    else:
        fig = plt.figure(figsize=figsize, facecolor=bgcolor)
        ax = fig.add_subplot(1, 1, 1, projection=crs, facecolor=bgcolor)
    fig.patch.set_facecolor(bgcolor)
    return fig, ax


def _draw_graticule_on_geoaxes(
    ax, step_deg, central_lon_rad, lw=0.5, color="white", ls="--", alpha=0.7, n_pts=360
):
    """
    Pre-projected graticule on a GeoAxes, WITHOUT calling ax.gridlines().
    ax.gridlines() clips each line through the CRS pipeline (~400 ms).
    Instead we project lines ourselves and call ax.plot() with
    transform=ax.projection so cartopy places them correctly.
    """
    native = ax.projection
    step = float(step_deg)
    lons_r = np.linspace(-180.0, 180.0, n_pts)
    lats_r = np.linspace(-89.99, 89.99, n_pts)
    for lat in np.arange(-90.0, 90.0 + 1e-9, step):
        lat_c = np.clip(lat, -89.99, 89.99)
        xp, yp = _mollweide_forward(
            lons_r, np.full_like(lons_r, lat_c), central_lon_rad
        )
        ax.plot(
            xp,
            yp,
            color=color,
            lw=lw * 1.6 if lat == 0.0 else lw,
            ls="-" if lat == 0.0 else ls,
            alpha=alpha,
            zorder=2,
            transform=native,
        )
    for lon in np.arange(-180.0, 180.0 + 1e-9, step):
        xp, yp = _mollweide_forward(np.full_like(lats_r, lon), lats_r, central_lon_rad)
        ax.plot(
            xp, yp, color=color, lw=lw, ls=ls, alpha=alpha, zorder=2, transform=native
        )


# ===========================================================================
# mollview
# ===========================================================================


def mollview(
    hpx_map,
    *,
    nest=False,
    title="",
    cmap="viridis",
    vmin=None,
    vmax=None,
    rot=0.0,
    ellipsoid="sphere",
    graticule=True,
    graticule_step=30.0,
    unit="",
    bgcolor="black",
    n_lon=1800,
    n_lat=900,
    norm=None,
    bad_color="gray",
    flip="geo",
    figsize=(14, 7),
    colorbar=True,
    hold=False,
    sub=None,
    coastlines=False,
    coastline_kwargs=None,
):
    """
    Display a HEALPix map in the Mollweide projection.

    Backend is selected automatically:

    * coastlines=False (default) -> fast pure-matplotlib path (~70 ms).
      cartopy is not imported.

    * coastlines=True -> cartopy path (~500 ms).
      cartopy is imported lazily and must be installed.
      The graticule still uses the fast pre-projected approach.

    Parameters
    ----------
    hpx_map : np.ndarray, shape (12 * 4**depth,)
        Input HEALPix map. RING by default; nest=True for NESTED.
    nest : bool, default False
    title : str
    cmap : str or Colormap, default "viridis"
    vmin, vmax : float or None
        Colour limits. Default: 2nd / 98th percentile.
    rot : float, default 0.0
        Central longitude in degrees.
    ellipsoid : str, default "sphere"
    graticule : bool, default True
    graticule_step : float, default 30.0
    unit : str -- colorbar label
    bgcolor : str, default "black"
    n_lon, n_lat : int, default 1800, 900
        Sampling grid resolution in Mollweide space.
    norm : Normalize or None
    bad_color : str, default "gray"
    flip : {"geo", "astro"}, default "geo"
        "astro" puts east to the left.
    figsize : (float, float), default (14, 7)
    colorbar : bool, default True
    hold : bool, default False
    sub : (nrows, ncols, idx) or None
    coastlines : bool, default False
        Overlay Natural Earth coastlines.
        Triggers the cartopy backend -- cartopy must be installed.
    coastline_kwargs : dict or None
        Kwargs forwarded to ax.add_feature(COASTLINE, ...).
        Only used when coastlines=True.
    """
    if flip not in ("astro", "geo"):
        raise ValueError(f"flip must be 'astro' or 'geo', got {flip!r}.")

    hpx_map = np.asarray(hpx_map, dtype=np.float64)
    depth = _depth_from_npix(hpx_map.size)
    central_lon_rad = math.radians(rot)

    cmap_obj = _build_cmap(cmap, bad_color)

    # -------------------------------------------------------------------
    if not coastlines:
        # ===============================================================
        # FAST PATH -- pure matplotlib, cartopy not imported
        # ===============================================================
        # Grid built in Mollweide (x, y) space via the inverse formula.
        # The resulting data_img is already in display coordinates, so
        # ax.imshow() needs no transform argument.
        xs = np.linspace(-_2SQRT2, _2SQRT2, n_lon)
        ys = np.linspace(-_SQRT2, _SQRT2, n_lat)
        if flip == "astro":
            xs = xs[::-1]
        xg, yg = np.meshgrid(xs, ys)

        lon_grid, lat_grid, valid = _mollweide_inverse(xg, yg, central_lon_rad)
        data_img = _rasterise(
            hpx_map, depth, lon_grid, lat_grid, valid, nest, ellipsoid
        )
        norm_obj = _build_norm(data_img, vmin, vmax, norm)

        fig, ax = _make_axes_plain(hold, sub, figsize, bgcolor)

        ax.imshow(
            data_img,
            origin="lower",
            extent=[-_2SQRT2, _2SQRT2, -_SQRT2, _SQRT2],
            cmap=cmap_obj,
            norm=norm_obj,
            interpolation="nearest",
            aspect="auto",
            zorder=1,
        )

        if graticule:
            _draw_graticule_moll(ax, graticule_step, central_lon_rad)

        _finalize_moll_axes(ax)

    else:
        # ===============================================================
        # CARTOPY PATH -- only when coastlines=True
        # cartopy is imported here, lazily.
        # ===============================================================
        # Grid built in PlateCarree (lon/lat) space so that the extent
        # [-180, 180, -90, 90] passed to imshow(transform=PlateCarree())
        # matches the actual sampling layout.
        # (Using a Mollweide-space grid here would cause misalignment
        # between the raster data and the coastline overlay.)
        lons_1d = np.linspace(-180.0, 180.0, n_lon)
        lats_1d = np.linspace(-90.0, 90.0, n_lat)
        if flip == "astro":
            lons_1d = lons_1d[::-1]
        lon_grid, lat_grid = np.meshgrid(lons_1d, lats_1d)
        valid = np.ones((n_lat, n_lon), dtype=bool)
        data_img = _rasterise(
            hpx_map, depth, lon_grid, lat_grid, valid, nest, ellipsoid
        )
        norm_obj = _build_norm(data_img, vmin, vmax, norm)

        ccrs, cfeature = _require_cartopy()

        crs = ccrs.Mollweide(central_longitude=rot)
        plate_cr = ccrs.PlateCarree()

        fig, ax = _make_axes_geo(crs, hold, sub, figsize, bgcolor)
        ax.set_global()

        # imshow with transform=PlateCarree(): cartopy reprojects at
        # render time from the uniform lon/lat grid to Mollweide.
        ax.imshow(
            data_img,
            origin="lower",
            extent=[-180, 180, -90, 90],
            transform=plate_cr,
            cmap=cmap_obj,
            norm=norm_obj,
            interpolation="nearest",
            aspect="auto",
            zorder=1,
        )

        # Coastlines
        ckw = {"linewidth": 0.6, "edgecolor": "white"}
        if coastline_kwargs:
            ckw.update(coastline_kwargs)
        ax.add_feature(cfeature.COASTLINE, **ckw)

        # Graticule: pre-projected (avoids ax.gridlines() ~400 ms cost)
        if graticule:
            _draw_graticule_on_geoaxes(ax, graticule_step, central_lon_rad)

    # -------------------------------------------------------------------
    # Common: title + colorbar
    # -------------------------------------------------------------------
    if title:
        ax.set_title(title, color="white", fontsize=11, pad=6)
    if colorbar:
        _add_colorbar(fig, ax, norm_obj, cmap_obj, unit)


# ===========================================================================
# mollgnomview
# ===========================================================================


def mollgnomview(
    hpx_map,
    lon_center,
    lat_center,
    *,
    nest=False,
    fov_deg=10.0,
    title="",
    cmap="viridis",
    vmin=None,
    vmax=None,
    ellipsoid="sphere",
    unit="",
    n_lon=800,
    n_lat=800,
    figsize=(7, 7),
    colorbar=True,
    hold=False,
    sub=None,
    coastlines=False,
    coastline_kwargs=None,
):
    """
    Local zoom in the gnomonic (tangent-plane) projection.

    Backend selection -- same logic as mollview:

    * coastlines=False -> fast pure-matplotlib path.
    * coastlines=True  -> cartopy path (requires cartopy).

    Parameters
    ----------
    hpx_map : np.ndarray
    lon_center, lat_center : float  -- centre in degrees
    fov_deg : float, default 10.0  -- square side-length in degrees
    n_lon, n_lat : int, default 800
    All other parameters: see mollview.
    """
    hpx_map = np.asarray(hpx_map, dtype=np.float64)
    depth = _depth_from_npix(hpx_map.size)
    half = fov_deg / 2.0

    # -------------------------------------------------------------------
    # Sampling grid + HEALPix lookup (shared by both backends)
    # -------------------------------------------------------------------
    xs = np.linspace(-half, half, n_lon)
    ys = np.linspace(-half, half, n_lat)
    xg, yg = np.meshgrid(xs, ys)

    lon_grid, lat_grid = _gnomonic_inverse(xg, yg, lon_center, lat_center)
    valid = np.ones((n_lat, n_lon), dtype=bool)
    data_img = _rasterise(hpx_map, depth, lon_grid, lat_grid, valid, nest, ellipsoid)

    norm_obj = _build_norm(data_img, vmin, vmax, None)
    cmap_obj = _build_cmap(cmap, "gray")

    # -------------------------------------------------------------------
    if not coastlines:
        # ===============================================================
        # FAST PATH
        # ===============================================================
        fig, ax = _make_axes_plain(hold, sub, figsize, "black")

        ax.imshow(
            data_img,
            origin="lower",
            extent=[-half, half, -half, half],
            cmap=cmap_obj,
            norm=norm_obj,
            interpolation="nearest",
            aspect="auto",
            zorder=1,
        )

        # Light tangent-plane graticule (straight lines -- no projection needed)
        step = max(1.0, round(fov_deg / 4.0))
        for v in np.arange(-half, half + 1e-9, step):
            ax.axhline(v, color="gray", lw=0.5, ls="--", alpha=0.6, zorder=2)
            ax.axvline(v, color="gray", lw=0.5, ls="--", alpha=0.6, zorder=2)
            lo_v, _ = _gnomonic_inverse(
                np.array([v]), np.array([0.0]), lon_center, lat_center
            )
            _, la_v = _gnomonic_inverse(
                np.array([0.0]), np.array([v]), lon_center, lat_center
            )
            ax.text(
                v,
                -half * 1.02,
                f"{lo_v[0]:.1f}deg",
                color="gray",
                fontsize=6,
                ha="center",
                va="top",
            )
            ax.text(
                -half * 1.02,
                v,
                f"{la_v[0]:.1f}deg",
                color="gray",
                fontsize=6,
                ha="right",
                va="center",
            )

        ax.set_xlim(-half * 1.08, half * 1.08)
        ax.set_ylim(-half * 1.08, half * 1.08)

    else:
        # ===============================================================
        # CARTOPY PATH
        # ===============================================================
        ccrs, cfeature = _require_cartopy()

        crs = ccrs.Gnomonic(central_latitude=lat_center, central_longitude=lon_center)
        plate_cr = ccrs.PlateCarree()

        fig, ax = _make_axes_geo(crs, hold, sub, figsize, "black")
        ax.set_extent(
            [
                lon_center - half,
                lon_center + half,
                lat_center - half,
                lat_center + half,
            ],
            crs=plate_cr,
        )

        ax.imshow(
            data_img,
            origin="lower",
            extent=[
                lon_center - half,
                lon_center + half,
                lat_center - half,
                lat_center + half,
            ],
            transform=plate_cr,
            cmap=cmap_obj,
            norm=norm_obj,
            interpolation="nearest",
            aspect="auto",
            zorder=1,
        )

        ckw = {"linewidth": 0.6, "edgecolor": "white"}
        if coastline_kwargs:
            ckw.update(coastline_kwargs)
        ax.add_feature(cfeature.COASTLINE, **ckw)

        # For gnomonic zoom, gridlines is acceptable (small FOV = few lines)
        ax.gridlines(
            crs=plate_cr,
            draw_labels=True,
            linewidth=0.5,
            color="gray",
            linestyle="--",
            alpha=0.6,
        )

    # -------------------------------------------------------------------
    # Common: title + colorbar
    # -------------------------------------------------------------------
    if title:
        ax.set_title(title, color="white", fontsize=11)
    if colorbar:
        _add_colorbar(fig, ax, norm_obj, cmap_obj, unit)
