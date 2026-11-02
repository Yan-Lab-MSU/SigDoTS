# Package index

## Run SigDoTS

[`sigdots()`](../reference/sigdots.md) runs the whole analysis and
returns the result;
[`sigdots_control()`](../reference/sigdots_control.md) holds the
settings that rarely need changing. Each page below explains the
function, every argument with its allowed values and default, and what
is returned.

- [`sigdots()`](sigdots.md) : Significance test for spatial domains
- [`sigdots_control()`](sigdots_control.md) : Settings of sigdots() that
  rarely need changing

## Prepare the input

[`sigdots()`](../reference/sigdots.md) takes a normalised expression
matrix with coordinates, a Seurat object or a SpatialExperiment object
directly. From raw counts,
[`sigdots_prepare()`](../reference/sigdots_prepare.md) makes the matrix.

- [`sigdots_prepare()`](sigdots_prepare.md) : Expression matrix for
  sigdots() from a counts matrix

## Read the result

What [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) show and how to read
the split tree and the tests table; the domain of every spot as a table
or written back into your object; the plots.

- [`print(`*`<sigdots>`*`)`](summary.sigdots.md)
  [`summary(`*`<sigdots>`*`)`](summary.sigdots.md)
  [`print(`*`<summary.sigdots>`*`)`](summary.sigdots.md) : Print and
  summarise a SigDoTS result
- [`domains()`](domains.md) : The domain of every spot
- [`add_domains()`](add_domains.md) : Add the domains to a Seurat or
  SpatialExperiment object
- [`plot(`*`<sigdots>`*`)`](plot.sigdots.md) : Plot a SigDoTS result

## Example data

- [`sigdots_example`](sigdots_example.md) : A small simulated dataset
  with three spatial domains

## Overview

- [`SigDoTS`](SigDoTS-package.md)
  [`SigDoTS-package`](SigDoTS-package.md) : SigDoTS: significance
  testing for spatial domains
