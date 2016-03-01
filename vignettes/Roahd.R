## ------------------------------------------------------------------------
# The number of observations in the functional dataset.
N = 3
# The number of points in the 1D grid where the functional data are measured.
P = 1e2
# The previous two variable names are used consistently throughout the tutorial
# and the package's documentation to indicate the number of observations and the
# grid size.

# The grid over which the functional dataset is defined
grid = seq( 0, 1, length.out = P )

Data = matrix( c( sin( 2 * pi * grid  ),
                  cos( 2 * pi * grid ),
                  4 * grid * ( 1 - grid )  ),
              nrow = N, ncol = P, byrow = TRUE )

fD = fData( grid, Data )

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

