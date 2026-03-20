"""
mollview — projection de Mollweide d'une carte HEALPix
via healpix-geo pour les conversions de coordonnées.

Équivalent à healpy.mollview, sans dépendance à healpy.
Supporte les ordres RING (défaut) et NESTED, ainsi que les ellipsoïdes
de référence (ex. WGS84) grâce à healpix-geo.

Dépendances :
    pip install healpix-geo numpy matplotlib
"""

from __future__ import annotations

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Ellipse
from matplotlib.axes import Axes
from typing import Optional

import healpix_geo 


# ---------------------------------------------------------------------------
# Utilitaire : déduction de depth depuis la taille de la carte
# ---------------------------------------------------------------------------

def _depth_from_npix(npix: int) -> int:
    """
    Déduit la profondeur HEALPix depuis le nombre de pixels.

    npix = 12 * 4**depth  →  depth = log2(sqrt(npix / 12))

    Raises
    ------
    ValueError si npix n'est pas une taille HEALPix valide.
    """
    if npix <= 0 or npix % 12 != 0:
        raise ValueError(
            f"Taille de carte invalide : {npix}. "
            "npix doit être un multiple de 12 de la forme 12 * 4**depth."
        )
    nside_sq = npix // 12          # = 4**depth = nside²
    log2_nside = math.log2(nside_sq) / 2.0
    depth = int(round(log2_nside))
    if 12 * (4**depth) != npix:
        raise ValueError(
            f"Taille de carte invalide : {npix}. "
            f"La valeur la plus proche serait depth={depth} "
            f"→ npix={12 * 4**depth}."
        )
    return depth


# ---------------------------------------------------------------------------
# Projection de Mollweide — formules directe / inverse
# ---------------------------------------------------------------------------

def _mollweide_inverse(
    x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """(x,y) normalisé [-2,2]×[-1,1]  →  (lon_deg, lat_deg)."""
    sin_theta = np.clip(y, -1.0, 1.0)
    theta = np.arcsin(sin_theta)
    two_theta = 2.0 * theta
    sin_lat = np.clip((two_theta + np.sin(two_theta)) / np.pi, -1.0, 1.0)
    lat_deg = np.degrees(np.arcsin(sin_lat))
    cos_theta = np.cos(theta)
    with np.errstate(invalid="ignore", divide="ignore"):
        lon_rad = np.where(
            cos_theta > 1e-12,
            (x / 2.0) * np.pi / cos_theta,
            0.0,
        )
    return np.degrees(lon_rad), lat_deg


def _mollweide_forward(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """(lon_deg, lat_deg)  →  (x, y) normalisé [-2,2]×[-1,1]."""
    lat_rad = np.radians(lat_deg)
    lon_rad = np.radians(lon_deg)
    target = np.pi * np.sin(lat_rad)
    theta = lat_rad.copy()
    for _ in range(max_iter):
        f = 2.0 * theta + np.sin(2.0 * theta) - target
        df = 2.0 + 2.0 * np.cos(2.0 * theta)
        delta = f / np.where(np.abs(df) > 1e-15, df, 1e-15)
        theta -= delta
        if np.max(np.abs(delta)) < tol:
            break
    return (2.0 / np.pi) * lon_rad * np.cos(theta), np.sin(theta)


# ---------------------------------------------------------------------------
# Grille image → lon/lat
# ---------------------------------------------------------------------------

def _make_mollweide_grid(
    width: int, height: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Renvoie (lon_grid, lat_grid, inside) pour une image width×height.
    lon_grid en degrés [-180,180], lat_grid en degrés [-90,90].
    inside : masque booléen — True si dans l'ellipse de Mollweide.
    """
    xs = np.linspace(-2.0, 2.0, width)
    ys = np.linspace(1.0, -1.0, height)   # nord en haut
    x_grid, y_grid = np.meshgrid(xs, ys)

    inside = (x_grid / 2.0) ** 2 + y_grid ** 2 <= 1.0

    lon_flat = np.zeros(width * height)
    lat_flat = np.zeros(width * height)
    idx = inside.ravel()
    lon_flat[idx], lat_flat[idx] = _mollweide_inverse(
        x_grid.ravel()[idx], y_grid.ravel()[idx]
    )
    return (
        lon_flat.reshape(height, width),
        lat_flat.reshape(height, width),
        inside,
    )


# ---------------------------------------------------------------------------
# Graticule
# ---------------------------------------------------------------------------

def _draw_graticule(
    ax: Axes,
    graticule_step_deg: float = 30.0,
    rot: float = 0.0,
    line_kwargs: Optional[dict] = None,
) -> None:
    """Trace méridiens et parallèles en projection Mollweide sur *ax*."""
    if line_kwargs is None:
        line_kwargs = dict(color="white", linewidth=0.5, linestyle="--", alpha=0.7)

    step = graticule_step_deg
    n_pts = 300

    # Méridiens
    lat_curve = np.linspace(-90.0, 90.0, n_pts)
    for lon_deg in np.arange(-180.0, 180.0 + 1e-6, step):
        lon_s = (lon_deg - rot + 180.0) % 360.0 - 180.0
        x, y = _mollweide_forward(np.full(n_pts, lon_s), lat_curve)
        ax.plot(x, y, **line_kwargs)

    # Parallèles
    lon_curve = np.linspace(-180.0, 180.0, n_pts)
    for lat_deg in np.arange(-90.0 + step, 90.0, step):
        x, y = _mollweide_forward(lon_curve, np.full(n_pts, lat_deg))
        ax.plot(x, y, **line_kwargs)

    # Équateur plus épais
    eq_kw = {**line_kwargs, "linewidth": 0.8, "linestyle": "-"}
    x_eq, y_eq = _mollweide_forward(lon_curve, np.zeros(n_pts))
    ax.plot(x_eq, y_eq, **eq_kw)


def _draw_graticule_labels(
    ax: Axes, step: float, rot: float, color: str = "white"
) -> None:
    """Annotations lon/lat du graticule."""
    # Longitudes (bas de l'ellipse)
    for lon_label in np.arange(-180.0 + step, 180.0, step):
        lon_s = (lon_label - rot + 180.0) % 360.0 - 180.0
        x_lbl, _ = _mollweide_forward(
            np.array([lon_s]), np.array([0.0])
        )
        ax.text(
            x_lbl[0], -1.05,
            f"{int(lon_label):+d}°" if lon_label != 0 else "0°",
            ha="center", va="top", fontsize=7, color=color, clip_on=False,
        )
    # Latitudes (bord gauche)
    for lat_label in np.arange(-90.0 + step, 90.0, step):
        if abs(lat_label) < 1:
            continue
        x_l, y_l = _mollweide_forward(
            np.array([-180.0]), np.array([lat_label])
        )
        ax.text(
            x_l[0] - 0.05, y_l[0],
            f"{int(lat_label):+d}°",
            ha="right", va="center", fontsize=7, color=color, clip_on=False,
        )


# ---------------------------------------------------------------------------
# Fonction principale : mollview
# ---------------------------------------------------------------------------

def mollview(
    hpx_map: np.ndarray,
    *,
    nest: bool = False,
    title: str = "",
    cmap: str | mcolors.Colormap = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    rot: float = 0.0,
    ellipsoid: str = "sphere",
    graticule: bool = True,
    graticule_step: float = 30.0,
    unit: str = "",
    bgcolor: str = "black",
    width_px: int = 1600,
    height_px: int = 800,
    norm: Optional[mcolors.Normalize] = None,
    bad_color: str = "gray",
    flip: str = "astro",
    figsize: tuple[float, float] = (14, 7),
    colorbar: bool = True,
    hold: bool = False,
    sub: Optional[tuple[int, int, int]] = None,
) -> None:
    """
    Affiche une carte HEALPix en projection de Mollweide.

    Équivalent à ``healpy.mollview``, sans dépendance à healpy.
    La profondeur HEALPix est déduite automatiquement de la taille de la carte.

    Parameters
    ----------
    hpx_map : np.ndarray, shape (12 * 4**depth,)
        Carte HEALPix. Ordre RING par défaut (comme healpy), ou NESTED si
        ``nest=True``.
    nest : bool, optional
        False (défaut) → carte en ordre RING.
        True            → carte en ordre NESTED.
    title : str
        Titre affiché au-dessus de la carte.
    cmap : str ou Colormap
        Palette de couleurs matplotlib. Défaut ``"viridis"``.
        Recommandé pour CMB : ``"RdBu_r"``.
    vmin, vmax : float, optional
        Limites de la colormap. Si None, utilise les percentiles 2/98.
    rot : float
        Longitude centrale de la carte en degrés (défaut 0°).
        Exemple : ``rot=180`` centre la carte sur le méridien 180°.
    ellipsoid : str
        Ellipsoïde de référence. Défaut ``"sphere"`` (identique à healpy).
        Exemples : ``"WGS84"``, ``"GRS80"``.
    graticule : bool
        Trace le graticule (méridiens + parallèles). Défaut True.
    graticule_step : float
        Pas du graticule en degrés. Défaut 30°.
    unit : str
        Unité affichée sous la colorbar.
    bgcolor : str
        Couleur de fond hors ellipse. Défaut ``"black"``.
    width_px, height_px : int
        Résolution de l'image rasterisée. Défaut 1600 × 800.
        Augmenter pour des cartes à haute résolution (depth ≥ 8).
    norm : matplotlib.colors.Normalize, optional
        Normalisation personnalisée (ex. ``LogNorm``). Prend le dessus sur
        vmin/vmax si fournie.
    bad_color : str
        Couleur pour les valeurs NaN. Défaut ``"gray"``.
    flip : str
        ``"astro"`` (est à gauche, convention astronomique, défaut) ou
        ``"geo"`` (est à droite, convention géographique).
    figsize : tuple (float, float)
        Taille de la figure en pouces. Utilisé uniquement si une nouvelle
        figure est créée (``hold=False`` et ``sub=None``).
    colorbar : bool
        Affiche la colorbar horizontale. Défaut True.
    hold : bool
        False (défaut) : crée une nouvelle figure, comme healpy.
        True           : dessine sur la figure/axes courant (``plt.gca()``).
    sub : tuple (nrows, ncols, idx), optional
        Positionne la carte dans un subplot de la figure courante.
        Exemple : ``sub=(2, 3, 4)`` → 2 lignes × 3 colonnes, panneau 4.
        Si fourni, ``hold`` est ignoré.

    Returns
    -------
    None
        Comme healpy.mollview, la fonction ne retourne rien.
        Utilisez ``plt.savefig()``, ``plt.show()`` ou ``plt.gcf()`` pour
        accéder à la figure courante.

    Examples
    --------
    >>> import numpy as np, matplotlib.pyplot as plt
    >>> from healpix_mollview import mollview
    >>>
    >>> depth  = 5
    >>> m_ring = np.random.default_rng(0).standard_normal(12 * 4**depth)
    >>>
    >>> # Carte simple — nouvelle figure automatique
    >>> mollview(m_ring, title="Carte RING", cmap="RdBu_r")
    >>> plt.show()
    >>>
    >>> # Deux cartes côte à côte dans la même figure
    >>> plt.figure(figsize=(18, 5))
    >>> mollview(m_ring, title="RING sphere",  sub=(1, 2, 1))
    >>> mollview(m_ring, title="RING rot=180", sub=(1, 2, 2), rot=180)
    >>> plt.show()
    >>>
    >>> # Carte NESTED, WGS84, centrée sur 180°
    >>> mollview(m_nested, nest=True, rot=180, ellipsoid="WGS84",
    ...          title="WGS84 NESTED", cmap="plasma")
    >>> plt.savefig("out.png", dpi=150, bbox_inches="tight", facecolor="black")
    """
    if flip not in ("astro", "geo"):
        raise ValueError(f"flip doit être 'astro' ou 'geo', pas '{flip}'.")

    hpx_map = np.asarray(hpx_map, dtype=np.float64)
    depth = _depth_from_npix(hpx_map.size)

    # ------------------------------------------------------------------
    # 1. Grille image → (lon, lat) en degrés
    # ------------------------------------------------------------------
    lon_grid, lat_grid, inside = _make_mollweide_grid(width_px, height_px)

    # ------------------------------------------------------------------
    # 2. Rotation : décalage de longitude
    # ------------------------------------------------------------------
    lon_rotated = (lon_grid + rot + 180.0) % 360.0 - 180.0

    # ------------------------------------------------------------------
    # 3. lon/lat → indice HEALPix NESTED (healpix-geo)
    #    puis conversion NESTED → RING si la carte est en ordre RING
    # ------------------------------------------------------------------
    idx_flat = inside.ravel()
    if nest:
        ipix_lookup = healpix_geo.nested.lonlat_to_healpix(
            lon_rotated.ravel()[idx_flat],
            lat_grid.ravel()[idx_flat],
            depth,
            ellipsoid=ellipsoid,
        )  # uint64, ordre NESTED
    else:
        ipix_lookup = healpix_geo.ring.lonlat_to_healpix(
            lon_rotated.ravel()[idx_flat],
            lat_grid.ravel()[idx_flat],
            depth,
            ellipsoid=ellipsoid,
        )  # uint64, ordre ring

    # ------------------------------------------------------------------
    # 4. Lecture des valeurs dans la carte
    # ------------------------------------------------------------------
    data_flat = np.full(width_px * height_px, np.nan)
    data_flat[idx_flat] = hpx_map[ipix_lookup]
    data_img = data_flat.reshape(height_px, width_px)

    # ------------------------------------------------------------------
    # 5. Normalisation couleur
    # ------------------------------------------------------------------
    valid_vals = data_flat[idx_flat]
    valid_vals = valid_vals[np.isfinite(valid_vals)]
    if norm is None:
        if vmin is None:
            vmin = float(np.percentile(valid_vals, 2)) if len(valid_vals) else 0.0
        if vmax is None:
            vmax = float(np.percentile(valid_vals, 98)) if len(valid_vals) else 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # ------------------------------------------------------------------
    # 6. Image RGBA
    # ------------------------------------------------------------------
    if isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap).copy()
    else:
        cmap_obj = cmap.copy() if hasattr(cmap, "copy") else cmap
    cmap_obj.set_bad(bad_color, alpha=1.0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    rgba = sm.to_rgba(data_img, bytes=False)   # (H, W, 4)

    rgba[~inside] = mcolors.to_rgba(bgcolor)   # fond hors ellipse

    if flip == "geo":
        rgba = rgba[:, ::-1, :]

    # ------------------------------------------------------------------
    # 7. Sélection / création de l'Axes
    #
    #   Priorité : sub  >  hold=True  >  hold=False (nouvelle figure)
    # ------------------------------------------------------------------
    if sub is not None:
        nrows, ncols, sidx = sub
        ax = plt.gcf().add_subplot(nrows, ncols, sidx, facecolor=bgcolor)
        fig = ax.get_figure()
    elif hold:
        ax = plt.gca()
        ax.set_facecolor(bgcolor)
        fig = ax.get_figure()
    else:
        fig = plt.figure(figsize=figsize, facecolor=bgcolor)
        ax = fig.add_subplot(111, facecolor=bgcolor)

    fig.patch.set_facecolor(bgcolor)

    # ------------------------------------------------------------------
    # 8. Rendu
    # ------------------------------------------------------------------
    ax.imshow(
        rgba, extent=[-2, 2, -1, 1], origin="upper",
        interpolation="nearest", aspect="equal",
    )
    ax.add_patch(Ellipse(
        (0, 0), width=4, height=2,
        edgecolor="white", facecolor="none", linewidth=1.0, zorder=3,
    ))
    if graticule:
        _draw_graticule(ax, graticule_step_deg=graticule_step, rot=rot)
        _draw_graticule_labels(ax, step=graticule_step, rot=rot)

    if title:
        ax.set_title(title, color="white", fontsize=11, pad=6)

    ax.set_xlim(-2.05, 2.05)
    ax.set_ylim(-1.12, 1.05)
    ax.axis("off")

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

    # Pas de return — comme healpy.mollview


# ---------------------------------------------------------------------------
# mollgnomview — zoom gnomonique local
# ---------------------------------------------------------------------------

def mollgnomview(
    hpx_map: np.ndarray,
    lon_center: float,
    lat_center: float,
    *,
    nest: bool = False,
    fov_deg: float = 10.0,
    title: str = "",
    cmap: str | mcolors.Colormap = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    ellipsoid: str = "sphere",
    unit: str = "",
    width_px: int = 800,
    height_px: int = 800,
    figsize: tuple[float, float] = (7, 7),
    colorbar: bool = True,
    hold: bool = False,
    sub: Optional[tuple[int, int, int]] = None,
) -> None:
    """
    Zoom local en projection gnomonique sur une carte HEALPix.

    Équivalent à ``healpy.gnomview``. La profondeur est déduite de la taille
    de la carte. Même convention hold/sub que mollview.

    Parameters
    ----------
    hpx_map : np.ndarray
        Carte HEALPix (RING par défaut, NESTED si ``nest=True``).
    lon_center, lat_center : float
        Centre de la vue en degrés.
    nest, fov_deg, title, cmap, vmin, vmax, ellipsoid, unit,
    width_px, height_px, figsize, colorbar, hold, sub :
        Voir ``mollview``.
    """
    hpx_map = np.asarray(hpx_map, dtype=np.float64)
    depth = _depth_from_npix(hpx_map.size)

    half = np.tan(np.radians(fov_deg / 2.0))
    xs = np.linspace(-half, half, width_px)
    ys = np.linspace(half, -half, height_px)
    xg, yg = np.meshgrid(xs, ys)

    lon_c = np.radians(lon_center)
    lat_c = np.radians(lat_center)
    rho   = np.sqrt(xg**2 + yg**2)
    c     = np.arctan(rho)
    cos_c, sin_c = np.cos(c), np.sin(c)

    lat_rad = np.arcsin(
        cos_c * np.sin(lat_c)
        + np.where(rho > 1e-12, yg * sin_c * np.cos(lat_c) / rho, 0.0)
    )
    lon_rad = lon_c + np.arctan2(
        xg * sin_c,
        rho * np.cos(lat_c) * cos_c - yg * np.sin(lat_c) * sin_c,
    )

    ipix_nested = lonlat_to_healpix(
        np.degrees(lon_rad).ravel(),
        np.degrees(lat_rad).ravel(),
        depth, ellipsoid=ellipsoid,
    )
    ipix_lookup = ipix_nested if nest else _nested_to_ring(ipix_nested, depth)
    data_img = hpx_map[ipix_lookup].reshape(height_px, width_px)

    valid_vals = data_img[np.isfinite(data_img)]
    if vmin is None:
        vmin = float(np.percentile(valid_vals, 2)) if len(valid_vals) else 0.0
    if vmax is None:
        vmax = float(np.percentile(valid_vals, 98)) if len(valid_vals) else 1.0
    norm_obj = mcolors.Normalize(vmin=vmin, vmax=vmax)

    if sub is not None:
        nrows, ncols, sidx = sub
        ax = plt.gcf().add_subplot(nrows, ncols, sidx, facecolor="black")
        fig = ax.get_figure()
    elif hold:
        ax = plt.gca()
        ax.set_facecolor("black")
        fig = ax.get_figure()
    else:
        fig, ax = plt.subplots(figsize=figsize, facecolor="black")

    fig.patch.set_facecolor("black")
    im = ax.imshow(
        data_img, origin="upper", cmap=cmap, norm=norm_obj,
        extent=[-fov_deg / 2, fov_deg / 2, -fov_deg / 2, fov_deg / 2],
        interpolation="nearest",
    )
    ax.set_facecolor("black")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("white")
    ax.set_xlabel("Δlon [deg]", color="white")
    ax.set_ylabel("Δlat [deg]", color="white")
    if title:
        ax.set_title(title, color="white", fontsize=11)
    if colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(colors="white")
        cbar.outline.set_edgecolor("white")
        if unit:
            cbar.set_label(unit, color="white")

