

# TESTING FDATA -----------------------------------------------------------

N = 1e2

P = 1e3
t0 = 0
t1 = 1

time_grid = seq( t0, t1, length.out = P )

Cov = exp_cov_function( time_grid, alpha = 0.3, beta = 0.4 )

Data = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )

test_that( 'Creation of fData object',
           expect_silent( fData( time_grid, Data ) ) )

fD = fData( time_grid, Data )

test_that( 'Plot of fData object',
           expect_silent( plot_fData( fData( time_grid, Data ),
                                      xlab = 'time', ylab = 'values',
                                      main = 'A functional dataset' ) ) )


# TESTING STATISTICS OPERATIONS -------------------------------------------

test_that( 'Mean of fData object',
           expect_equal( as.numeric( mean_fData( fD )$values ),
                         colMeans( fD$values ) ) )

test_that( 'Median of fData obejct - MBD',
           expect_equal( as.numeric( median_fData( fD )$value ),
                         fD$values[ which.max( MBD( fD$values ) ), ]) )

test_that( 'Median of fData obejct - MHRD',
           expect_equal( as.numeric( median_fData( fD, type = 'MHRD' )$value ),
                         fD$values[ which.max( MHRD( fD$values ) ), ]) )


# meanfD = mean( fD )
# medfD = median.fData( fD )

# quartz()
# plot( fD, xlab = 'time', ylab = 'value', main = 'Basic stastics of fData objects')
# plot( meanfD, add = T, lwd = 2, lty = 2, col = 'darkblue' )
# plot( medfD, add = T, lwd = 2, lty = 2, col = 'darkred' )


# TESTING ALGEBRAIC OPERATIONS --------------------------------------------

fD = fData( seq( 0, 1, length.out = 10 ), values = matrix( seq( 1, 10 ),
                                                           nrow = 21, ncol = 10,
                                                           byrow = TRUE ) )

test_that( 'Sum of fData and raw vector',
           expect_equal( sum( ( fD + 1 : 10 )$values -
                                matrix( 2 * seq( 1, 10 ),  nrow = 21, ncol = 10,
                                        byrow = TRUE ) ), 0 ) )

test_that( 'Sum of fData and raw vector',
           expect_equal( sum( ( fD - 1 : 10 )$values ), 0 ) )

test_that( 'Sum of fData and array',
           expect_equal( sum( ( fD + array( 1, dim = c( 1, 10 ) ) )$values -
                                matrix( seq( 2, 11 ),
                                        nrow = 21, ncol = 10, byrow = TRUE ) ),
                         0 ) )

fD.2 = fData( seq( 0, 1, length.out = 11 ), values = matrix( seq( 1, 11 ),
                                                             nrow = 21, ncol = 11,
                                                             byrow = TRUE ) )

test_that( 'Sum of two compliant fData objects',
           expect_equal( sum( ( fD + fD -
                                  matrix( 2 * seq( 1, 10 ),
                                          nrow = 21, ncol = 10,
                                          byrow = TRUE ) )$values ), 0  ) )

test_that( 'Sum of two noncompliant fData objects',
           expect_error( fD + fD.2, regexp = 'Error.*') )

test_that( 'Product of fData object with scalar',
           expect_equal( fD * 2,  fD + fD ) )

test_that( 'Division of fData object by scalar',
           expect_equal( ( fD * 4 ) / 2, fD + fD ) )


# TESTING FUNCTIONAL DATA SUBSETTING --------------------------------------

test_that( 'Functional data subsetting - case 1',
           expect_identical( fD[ 1, ],
                             fData( seq( 0, 1, length.out = 10 ), 1 : 10 ) ) )

test_that( 'Functional data subsetting - case 2',
           expect_identical( fD[ , 1:2 ],
                             fData( seq( 0, 1, length.out = 10 )[1:2],
                                    matrix( 1 : 2,
                                            nrow = 21, ncol = 2,
                                            byrow = T ) ) ) )

test_that( 'Functional data subsetting - case 2',
           expect_identical( fD[ 1:2, 1:2, as_fData = FALSE ],
                             matrix( seq(1,2),
                                     nrow = 2, ncol = 2,
                                     byrow = TRUE ) ) )


# TESTING MFDATA ----------------------------------------------------------

N = 1e2

P = 1e3

t0 = 0
t1 = 1

time_grid = seq( t0, t1, length.out = P )

Cov = exp_cov_function( time_grid, alpha = 0.3, beta = 0.4 )

Data_1 = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )
Data_2 = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )

test_that( 'Creation of mfData object',
           expect_silent( mfData( time_grid, list( Data_1, Data_2 ) ) ) )

test_that( 'Plot of mfData object',
           expect_silent( plot_mfData( mfData( time_grid, list( Data_1, Data_2 ) ),
                                       xlab = 'time', ylab = list( 'values',
                                                                   'values' ),
                                       main = list( 'First Component',
                                                    'Second Component' ) ) ) )
