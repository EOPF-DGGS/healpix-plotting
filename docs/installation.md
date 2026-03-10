# Installation

This guide explains you how to install healpix-plot on your system.

## Requirements

- Python **≥ 3.13**

## Install

### Via conda

::::{tab-set}

:::{tab-item} conda

```bash
conda install -c conda-forge healpix-plot
```

:::

:::{tab-item} mamba

```bash
mamba install -c conda-forge healpix-plot
```

:::

:::{tab-item} pixi

```bash
pixi add healpix-plot
```

:::

::::

### Via pip

::::{tab-set}

:::{tab-item} pip

```bash
pip install healpix-plot
```

:::

:::{tab-item} uv

```bash
uv add healpix-plot
```

:::

::::

## Verify

```python
import healpix_plot

print(healpix_plot.__version__)
```
