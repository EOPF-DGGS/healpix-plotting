# healpix-plot

`healpix-plot` prioritises **getting a usable figure quickly** over perfectly accurate cell-geometry rendering. It rasterises HEALPix data via nearest-neighbour resampling onto a regular lon/lat grid and renders the result with Cartopy's `imshow`.

Unlike astronomy-focused HEALPix tools, this library is built with **Earth observation and geoscience in mind**: the underlying coordinate operations are provided by [healpix-geo](https://healpix-geo.readthedocs.io/en/latest/), which supports geodetically correct reference ellipsoids such as WGS84.

```{toctree}
---
maxdepth: 2
caption: User guide
hidden: true
---
installation
user-guide/index
tutorials/index
```

```{toctree}
---
maxdepth: 2
caption: Reference
hidden: true
---
api
terminology
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Tutorials
:link: tutorials/index
:link-type: doc
Step-by-step notebooks: from your first global map to RGB composites and regional subsets.
:::

:::{grid-item-card} API Reference
:link: api
:link-type: doc
Complete documentation for every public function and class.
:::

::::
