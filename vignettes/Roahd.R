## ----fData, fig.align='center', fig.height=4, fig.width=4, cache=TRUE----

library( roahd )

# The number of observations in the functional dataset.
N = 5
# The number of points in the 1D grid where the functional data are measured.
P = 1e2
# The previous two variable names are used consistently throughout the tutorial
# and the package's documentation to indicate the number of observations and the
# grid size.

# The grid over which the functional dataset is defined
grid = seq( 0, 1, length.out = P )

# Creating the values of the functional dataset
Data = matrix( c( sin( 2 * pi * grid  ),
                  cos( 2 * pi * grid ),
                  4 * grid * ( 1 - grid ),
                  tan( grid ),
                  log( grid ) ),
              nrow = N, ncol = P, byrow = TRUE )

# Building an fData object
# The constructor takes a grid and a matrix-like structure for data values
# (see help for more details on how to use the constructor)
fD = fData( grid, Data )

# Inspecting the structure of an fData object
str( fD )

plot( fD, main = 'Univariate FD', xlab = 'time [s]', ylab = 'values', lwd = 2 )

## ----mfData, cache = TRUE, fig.width = 7, fig.height = 4, fig.align = 'center'----
# Creating some values for first component of the dataset
Data_1 = t( sapply( runif( 10, 0, 4 ), 
                    function( phase ) sin( 2 * pi * grid + phase ) ) )

# Creating some values of functions for  
Data_2 = t( sapply( runif( 10, 0, 4 ), 
                    function( phase ) log( grid + phase ) ) )

# Building an fData object
# The constructor takes a grid and a list of matrix-like structures for data values,
# each one representing the data values of a single component of the dataset 
# (i.e. D_{,,k}, k = 1, ... L ).
# (see help for more details on how to use the constructor)
mfD = mfData( grid, list( Data_1, Data_2 ) )

str( mfD )

# Each component of the mfData object is a fData object 
sapply( mfD$fDList, class )


plot( mfD, lwd = 2, main = 'Multivariate FD',
      xlab = 'time', ylab = list( 'Values 1', 'Values 2' ))

## ----mfDataComponents, eval = TRUE, fig.width = 4, fig.height = 4, fig.align = 'center'----
plot( mfD$fDList[[ 1 ]], main = 'First component',
      xlab = 'time', ylab = 'Values', lwd = 2 )

## ----fDataToMfData, eval = FALSE-----------------------------------------
#  fD_1 = fData( grid, Data_1 )
#  
#  fD_2 = fData( grid, Data_2 )
#  
#  mfD = as.mfData( list( fD_1, fD_2 ) )
#  
#  mfD = as.mfData( lapply( 1 : 10, function( i )( fD_1 ) ) )

## ----fData_subset_1, eval = FALSE----------------------------------------
#  
#  # Subsetting fData and returning result in matrix form
#  fD[ 1 , 1, as_fData = FALSE ]
#  
#  fD[ 1, , as_fData = FALSE ]
#  
#  fD[ 2, 10 : 20, as_fData = FALSE ]
#  
#  fD[ , 10, as_fData = FALSE ]

## ----fData_subset_2, fig.align='center', fig.width=7, fig.height=4, cache=TRUE----
# As default behaviour the subset is returned in fData form
par( mfrow = c(1,2) )
plot( fD, main = 'Original dataset', lwd = 2 )
plot( fD[ , 1 : 20 ], main = 'Zooming in', lwd = 2 )

## ----fData_algebra_1, eval = FALSE---------------------------------------
#  fD + fD
#  
#  fD + matrix( 1, nrow = N, ncol = P )
#  
#  fD + array( 2, dim = c( N, P ) )
#  
#  fD + 1 : P

## ----fData_algebra_2, eval = FALSE---------------------------------------
#  fD * 2
#  
#  fD / 3
#  
#  fD * ( 1 : N )
#  
#  fD / ( 1 : N )

## ----BD_MBD_1, eval = FALSE, cache = TRUE--------------------------------
#  
#    BD( fD )
#    BD( fD$values )
#  
#    MBD( fD )
#    MBD( fD$values )
#  
#    MBD( fD, manage_ties = TRUE )
#    MBD( fD$values, manage_ties = TRUE )

## ----HRD, eval = FALSE, cache = TRUE-------------------------------------
#  
#  HRD( fD )
#  HRD( fD$values )
#  
#  MHRD( fD )
#  MHRD( fD$values )
#  

## ----multiMBD, eval = FALSE, cache = TRUE--------------------------------
#  
#    multiBD( mfD, weights = 'uniform' )
#    multiMBD( mfD, weights = 'uniform', manage_ties = FALSE )
#  
#    multiBD( mfD, weights = c( 0.6, 0.4) )
#    multiMBD( mfD, weights = c( 0.7, 0.3 ), manage_ties = FALSE )
#  
#    multiBD( list( fD_1$values, fD_2$values ), weights = c( 0.6, 0.4) )
#    multiMBD( list( fD_1$values, fD_2$values ), weights = c( 0.7, 0.3 ), manage_ties = FALSE )

## ----mean_1, eval = FALSE, cache = TRUE----------------------------------
#  
#    mean( fD )
#    mean( mfD)
#  
#    median_fData( fD, type = 'MBD' )
#    median_fData( fD, type = 'MHRD' )
#  

## ----mean_2, eval = TRUE, fig.height=3, fig.width=6, fig.align='center'----

  par( mfrow = c(1,2) )
  plot( fD, main = 'Mean', lwd = 2 )
  plot( mean( fD ), add = TRUE, lwd = 2, col = 'black', lty = 2 )
  
  plot( fD, main = 'Median', lwd = 2 )
  plot( median_fData( fD, type = 'MBD' ), add = TRUE, lwd = 2, lty = 2, col = 'black' )

## ----indexes, eval = FALSE-----------------------------------------------
#  
#    # Calling on fData objects
#    EI( fD )
#    MEI( fD )
#  
#    HI( fD )
#    MHI( fD )
#  
#    # Calling on 2D matrix type objects
#    EI( fD$values )
#    EI( matrix( rnorm( 20 ), nrow = 4, ncol = 5 ) )
#  

## ----cor_kendall_1, eval = TRUE, cache = TRUE----------------------------
  N = 10
  P = 1e3

  grid = seq( 0, 1, length.out = P )
  
  Data_1 = t( sapply( 1 : N, function( i )( sin( 2 * pi * grid ) + i ) ) )
  # Monotone nonlinear transformation of data
  Data_2 = Data_1^3
  
  mfD = mfData( grid, list( Data_1, Data_2 ) )

## ----cor_kendall_2, eval = TRUE, cache = TRUE, collapse = TRUE, fig.align = 'center', fig.height = 4, fig.width = 6----
  plot( mfD, main = list( 'Comp. 1', 'Comp. 2') )  

  # Kendall correlation of monotonically dependent data is exactly 1
  cor_kendall( mfD, ordering = 'max' )
  cor_kendall( mfD, ordering = 'area' )

## ----Spearman, eval = TRUE, cache = TRUE, collapse = TRUE----------------
# Spearman correlation of monotonically dependent data is exactly 1
  cor_spearman( mfD, ordering = 'MEI' )
  cor_spearman( mfD, ordering = 'MHI' )

## ----simulation_univariate, execute = FALSE, collapse = TRUE, fig.align = 'center', fig.width = 7, fig.height = 4----
  N = 50
  P = 1e3

  grid = seq( 0, 1, length.out = P )
  
  Cov = exp_cov_function( grid, alpha = 0.2, beta  = 0.3 )
  
  Data = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
  
  fD = fData( grid, Data )
  
  plot( fD, main = 'Gaussian fData', xlab = 'grid', lwd = 2)

## ----simulation_multivariate, execute = FALSE, collapse = TRUE, fig.align = 'center', fig.width = 7, fig.height = 4----
  N = 10
  P = 1e3

  grid = seq( 0, 1, length.out = P )
  
  Cov_1 = exp_cov_function( grid, alpha = 0.1, beta  = 0.5 )
  Cov_2 = exp_cov_function( grid, alpha = 0.5, beta = 0.1)
  
  centerline = matrix( c( sin( 2 * pi * grid ),
                          cos( 2 * pi * grid ) ), nrow = 2, byrow = TRUE )
  
  Data = generate_gauss_mfdata( N, 2, centerline, 0.8, list( Cov, Cov ) )
  
  mfD = mfData( grid, Data )
  
  plot( mfD, main = list( 'Comp.1', 'Comp. 2'), xlab = 'grid', lwd = 2)

## ----fbplot_fData, fig.align='center', fig.height=5, fig.width=7, cache=TRUE, results='hide'----

  set.seed(1618)

  N = 1e2
  P = 1e2
  
  grid = seq( 0, 1, length.out = P )
  
  Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
  
  Data = generate_gauss_fdata( N, sin( 2 * pi * grid ), Cov )
  
  fD = fData( grid, Data )
  
  Data = generate_gauss_mfdata( N, 2, matrix( sin( 2 * pi * grid ), nrow = 2, ncol = P, byrow = TRUE ), 0.6, listCov = list( Cov, Cov ) )

  mfD = mfData( grid, Data )
  
  fbplot( fD, main = 'Fbplot', Fvalue = 3.5 )

  fbplot( mfD, main = list( 'Comp. 1', 'Comp. 2' ), Fvalue = 3.5 )

## ----fbplot_fData_adjust, eval = TRUE, collapse = TRUE, fig.align = 'center', fig.width = 7, fig.height = 4, cache = TRUE----
  fbplot( fD, adjust = list( N_trials = 20, trial_size = N, TPR = 0.007, F_min = 0.1, F_max = 20 ), xlab = 'grid', ylab = 'values', main = 'Adjusted functional boxplot' )

## ----outliergram, cache = FALSE, fig.align = 'center', fig.width = 7, fig.height = 4----
set.seed( 1618 )

N = 100
P = 200
N_extra = 4

grid = seq( 0, 1, length.out = P )

Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.5 )

Data = generate_gauss_fdata( N, sin( 4 * pi * grid ), Cov )

Data_extra = array( 0, dim = c( N_extra, P ) )

Data_extra[ 1, ] = generate_gauss_fdata( 1, sin( 4 * pi * grid + pi / 2 ), Cov )

Data_extra[ 2, ] = generate_gauss_fdata( 1, sin( 4 * pi * grid - pi / 2 ), Cov )

Data_extra[ 3, ] = generate_gauss_fdata( 1, sin( 4 * pi * grid + pi/ 3 ), Cov )

Data_extra[ 4, ] = generate_gauss_fdata( 1, sin( 4 * pi * grid - pi / 3), Cov )

fD = fData( grid, rbind( Data, Data_extra ) )

outliergram( fD, display = TRUE )
outliergram( fD, Fvalue = 5, display = TRUE )

## ----outliergram_adjusted, fig.alicn = 'center', fig.width = 7, fig.height = 4----
outliergram( fD, adjust = list( N_trials = 5, trial_size = 5 * nrow( Data ), TPR = 0.01, VERBOSE = FALSE ), display = TRUE )

