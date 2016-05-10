# Changelog



Here's a list of what is changed in this update of __roahd__:

### Major fixes 

1) Modified the check of the grid provided to build fData objects. 
Since support is provided only for evenly spaced grids, a check is needed before building an `fData` object.
Before it was:

```r
all( abs( diff( unique( diff( grid  )  )  )  ) < 1e-14  )
```

Now it is:

```r
max( diff( unique( diff( grid  )  )  )  ) / diff( range( grid  )  ) < 1e-13
```

which is much more robust in practical cases.

2) Extended README.md

3) Added `cov_fun` method to compute covariance and cross-covariance functions
for either univariate or multivariate functional data. Implemented the `S3` class
`Cov` and plotting specialisation `plot.Cov`, wrapping `graphics::image`.



### Minor fixes 

1) Fixed typos in documentation

2) Fixed typos in vignette

3) Added [Travis](https://travis-ci.org/ntarabelloni/roahd) and [Codecov](https://codecov.io/gh/ntarabelloni/roahd) support

4) Modified the default parameter value for `trial_size` in `fbplot` from `Data$N` 
to `8 * Data$N`.

5) Added check to `fbplot` and `outliergram` that raises warnings when parameters
different than those supported are provided through `adjust` argument.
