# roahd 1.4.2

## Upgrades

## Minor updates

* Switch from Travis to Github Actions for continuous integration.
* Setup automatic `R CMD check` on Windows, macOS and Linux for both the latest
release and the development version of R.
* Setup automatic deployment of a [website](https://astamm.github.io/roahd/) for
the package that references a package introduction, its help and all vignettes.
* Setup automatic computation of [test coverage](https://codecov.io) and
report to both the [Github page](https://github.com/astamm/roahd) and
[website](https://astamm.github.io/roahd/) of the package.
* Added CRAN status badge to `README`.
* New package maintainer.

## Fixes

### Major fixes

* Updated all `matrix` class checks for compliance with R-4.0 in which the
`matrix` class inherits from the `array` class.

### Minor fixes

* Fixed typos in doc, vignette and README.
* Fixed bug in `fbplot()` display.

# roahd 1.4.1

## Fixes

### Minor fixes

* Fixed dependency error on a new version of **scales** package that breaks the
use of multivariate `fbplot` in the corner-case of zero outliers.

# roahd 1.4

## Upgrades

### Major upgrades

* Extended Spearman's correlation coefficient computation for multivariate
datasets with more than two components.
* Added bootstrap-based computation of Spearman's correlation coefficient bias
and standard deviation.
* Added methods to provide bootstrap-based confidence intervals on Spearman's
coefficients for two univariate functional datasets or a multivariate functional
dataset.
* Added a bootstrap-based test on Spearman's correlation coefficient for two
multivariate functional datasets.
* Added an outliergram version (without graphical display of original data) of
multivariate functional datasets.
* Added example multivariate functional datasets of ECG signals.

### Minor updates

* Added two convenience functions to append compatible functional datasets
(univariate or multivariate).
* Added a `[` operator overload for multivariate functional dataset
representation `mfData`.

## Fixes

### Major fixes

* Fixed bug in `cor_spearman()` function. Now the standard Spearman correlation
is not computed on ranks of MHI/MEI, but on MHI/MEI itself. The difference is
very small, but allows for full reproducibility of the results in the original
paper.

### Minor fixes

* Fixed typos in doc.
* Standardized formulas for the application of F inflations in outliergram and
boxplot.

# roahd 1.2

## Fixes

### Major fixes

* Removed check for uniformity in the grid of `fData()` and `mfData()`
constructor.
* Added the possibility to subset `fData` in time with logical vectors.
* Fixes in methods `BD`, `BD_relative`, `HI` and `EI`: the previous
computational technique was based on arguments from the popular reference "Exact
fast computation of band depth for large functional datasets: How quickly can
one million curves be ranked?" by Sun, Genton and Nychka, which in the case of
BD, and HI/EI are wrong. Now the implementation exploited sticks to the
definition, at the cost of a higher computational burden (and thus, time to
complete the computation).

# roahd 1.1

## Fixes

### Major fixes 

* Modified the check of the grid provided to build fData objects. Since support is provided only for evenly spaced grids, a check is needed before building an `fData` object. Before it was:

```r
all(abs(diff(unique(diff(grid)))) < 1e-14)
```

Now it is:

```r
max(diff(unique(diff(grid)))) / diff(range(grid)) < 1e-13
```

which is much more robust in practical cases.

* Extended `README.md`.
* Added `cov_fun` method to compute covariance and cross-covariance functions
for either univariate or multivariate functional data. Implemented the `S3` class
`Cov` and plotting specialization `plot.Cov`, wrapping `graphics::image`.

### Minor fixes 

* Fixed typos in documentation.
* Fixed typos in vignette.
* Added [Travis](https://travis-ci.org/ntarabelloni/roahd) and
[Codecov](https://codecov.io/gh/ntarabelloni/roahd) support.
* Modified the default parameter value for `trial_size` in `fbplot` from
`Data$N` to `8 * Data$N`.
* Added check to `fbplot` and `outliergram` that raises warnings when parameters
different than those supported are provided through `adjust` argument.

# roahd 1.0

* Initial release of the package.
