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


### Minor fixes 

1) Fixed typos in documentation

2) Fixed typos in vignette

3) Added [Travis](https://travis-ci.org/ntarabelloni/roahd) and [Codecov](https://codecov.io/gh/ntarabelloni/roahd) support
