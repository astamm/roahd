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
#  # Subseting fData and returning result in matrix form
#  fD[ 1 , 1, as_fData = FALSE ]
#  
#  fD[ 1, , as_fData = FALSE ]
#  
#  fD[ 2, 10 : 20, as_fData = FALSE ]
#  
#  fD[ , 10, as_fData = FALSE ]

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

## ----outliergram, cache = FALSE, fig.align = 'center', fig.width = 7, fig.height = 4----
set.seed( 1618 )

N = 200
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
outliergram( fD, Fvalue = 10, display = TRUE )

## ----outliergram_adjusted, fig.alicn = 'center', fig.width = 7, fig.height = 4----
outliergram( fD, adjust = list( N_trials = 10, trial_size = 5 * nrow( Data ), TPR = 0.01, VERBOSE = FALSE ), display = TRUE )

