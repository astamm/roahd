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
# The constructor takes a grid and a list of matrix-like structures for data values
# (see help for more details on how to use the constructor)
mfD = mfData( grid, list( Data_1, Data_2 ) )

str( mfD )

# Each component of the mfData object is a fData object 
sapply( mfD$fDList, class )


plot( mfD, lwd = 2, main = 'Multivariate FD',
      xlab = 'time', ylab = list( 'Values 1', 'Values 2' ))


## ----mfDataComponents, fig.width = 4, fig.height = 4, fig.align = 'center'----

plot( mfD$fDList[[ 1 ]], main = 'First component', 
      xlab = 'time', ylab = 'Values', lwd = 2 )

# S3 method for mean of univariate functional data
plot( mean( mfD$fDList[[ 1 ]] ), lwd = 2, col = 'black', add = TRUE )

## ----fDataToMfData, eval = FALSE-----------------------------------------
#  
#  fD_1 = fData( grid, Data_1 )
#  
#  fD_2 = fData( grid, Data_2 )
#  
#  mfD = as.mfData( list( fD_1, fD_2 ) )
#  
#  mfD = as.mfData( lapply( 1 : 10, function( i )( fD_1 ) ) )

## ----fData_subset_1, eval = FALSE----------------------------------------
#  
#  # Subseting fData and returning result in matrix form
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

## ----BD_MBD_1, eval = TRUE, cache = TRUE---------------------------------

  BD( fD )
  
  BD( fD$values )
  
  MBD( fD )
  
  MBD( fD )
  
  MBD( fD, manage_ties = TRUE )
  
  MBD( fD$values )

