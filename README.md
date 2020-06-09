
# roahd

<!-- badges: start -->
[![check-standard](https://github.com/astamm/roahd/workflows/R-CMD-check/badge.svg)](https://github.com/astamm/roahd/actions)
[![test-coverage](https://github.com/astamm/roahd/workflows/test-coverage/badge.svg)](https://github.com/astamm/roahd/actions)
[![codecov](https://codecov.io/gh/astamm/roahd/branch/master/graph/badge.svg)](https://codecov.io/gh/astamm/roahd)
[![pkgdown](https://github.com/astamm/roahd/workflows/pkgdown/badge.svg)](https://github.com/astamm/roahd/actions)
<!-- badges: end -->

Package __roahd__ (_Robust Analysis of High-dimensional Data_) allows to use
a set of statistical tools for the _exploration_ and _robustification_ of
univariate and multivariate __functional datasets__ through the use of depth-based
statistical methods.


In the implementation of functions special attention was put to their efficiency,
so that they can be profitably used also for the analysis of high-dimensional
datasets.

(_For a full-featured description of the package, please turn to the Vignette_)

## `fData` and `mfData` objects

A simple `S3` representation of functional data object, `fData`,
allows to encapsulate the important features of univariate functional datasets (like the
grid of the dependent variable, the pointwise observations etc.):

```r
# Grid representing the dependent variable
grid = seq( 0, 1, length.out = 100 )

# Pointwise-measurements of the functional dataset
Data = matrix( c( sin( 2 * pi * grid ),
                  cos ( 2 * pi * grid ),
                  sin( 2 * pi * grid + pi / 4 ) ), ncol = 100, byrow = TRUE )

# S3 object encapsulating the univariate functional dataset            
fD = fData( grid, Data )

# S3 representation of a multivariate functional dataset
mfD = mfData( grid, list( 'comp1' = Data, 'comp2' = Data ) )
```
Also, this allows to exploit simple calls to customised functions which
simplify the exploratory analysis:

```r
# Algebra of fData objects
fD + 1 : 100
fD * 4

fD_1 + fD_2

# Subsetting fData objects (providing other fData objects)
fD[ 1, ]
fD[ 1, 2 : 4]

# Smaple mean and (depth-based) median(s)
mean( fD )
mean( fD[ 1, 10 : 20 ] )
median_fData( fD, type = 'MBD' )

# Plotting functions
plot( fD )
plot( mean( fD ), add = TRUE )

plot( fD[ 2:3, :] )
```


## Robust methods for functional data analysis

A part of the package is specifically devoted to the computation of depths and
other statistical indexes for functional data:

  - Band Dephts and Modified Band Depths,
  - Modified band depths for multivariate functional data,
  - Epigraph and Hypograph indexes,
  - Spearman and Kendall's correlation indexes for functional data.
  - Confidence intervals and tests on Spearman's correlation coefficients for univariate andmultivariate functional data.

These also are the core of the visualization/robustification tools like
functional boxplot (`fbplot`) and outliergram (`outliergram`), allowing
the visualization and identification of amplitude/shape outliers.

Thanks to the functions for the simulation of synthetic functional datasets,
both `fbplot` and `outliergram` procedures can be auto-tuned to the dataset
at hand, in order to control the true positive outliers rate.
