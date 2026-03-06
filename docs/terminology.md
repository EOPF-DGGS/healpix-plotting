# Terminology

## Data formats

```{glossary}
HEALPix
  Hierarchical Equal Area isoLatitude Pixelization. Divides the sphere into equal-area cells identified by ``(level, cell_id)``. healpix-plotting uses `healpix-geo` for all HEALPix ↔ lon/lat conversions.

Cell ID
  Integer identifier for a single HEALPix cell. Its meaning depends on the indexing scheme (*nested*, *ring*, or *zuniq*).

Level
  HEALPix resolution parameter. At level ``k`` there are ``12 × 4^k`` cells. Range: 0 (12 cells) to 29 (~3 × 10^17 cells).

Nested scheme
  HEALPix indexing where child cells have IDs derived from their parent. Efficient for hierarchical operations.

Ring scheme
  HEALPix indexing ordered by iso-latitude rings. Required by some harmonic-transform libraries (e.g. HEALPy).

Zuniq scheme
  Compact unique identifier encoding both level and cell ID in a single integer.

WGS84
  World Geodetic System 1984 — the standard GPS ellipsoid (EPSG:4326).
```

## Resampling

```{glossary}
Sampling grid
  The regular pixel grid that HEALPix data is resampled onto before rendering. Described by a shape (pixels), a spatial extent (degrees), and optionally an affine transform.

Nearest-neighbour resampling
  Each output pixel takes the value of the closest HEALPix cell centre. Fast, preserves sharp boundaries.

Bilinear resampling
  Each output pixel is a distance-weighted blend of neighbouring cell values. Slower, produces smoother output for continuous fields.

Aggregation
  When multiple cell IDs map to the same cell, their values are combined with a reduce function (mean, max, min, sum, …) before resampling.

Background value
  Fill value used for output pixels that fall outside the spatial coverage of the input ``cell_ids``.
```

## Projections

```{glossary}
Cartopy CRS
  Coordinate Reference System object from the `cartopy` library.

Mollweide
  Equal-area pseudocylindrical map projection.

PlateCarree
  Simple equirectangular projection (lon/lat as x/y).
```
