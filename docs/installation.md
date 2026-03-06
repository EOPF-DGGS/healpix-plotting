# Installation

This guide explains you how to install healpix-plotting on your system.

## Requirements

- Python **≥ 3.13**

## Install

### Via conda

::::{tab-set}

:::{tab-item} conda

```bash
conda install -c conda-forge healpix-plotting
```

:::

:::{tab-item} mamba

```bash
mamba install -c conda-forge healpix-plotting
```

:::

:::{tab-item} pixi

```bash
pixi add healpix-plotting
```

:::

::::

### Via pip

::::{tab-set}

:::{tab-item} pip

```bash
pip install healpix-plotting
```

:::

:::{tab-item} uv

```bash
uv add healpix-plotting
```

:::

::::

## Verify

```python
import healpix_plotting

print(healpix_plotting.__version__)
```
