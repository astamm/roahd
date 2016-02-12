
# TESTING OUTLIERGRAM -----------------------------------------------------

set.seed( 1618 )

N = 200
N_outliers = 4
P = 200

time_grid = seq( 0, 1, length.out = P )


Cov = exp_cov_function( time_grid, alpha = 0.2, beta = 0.8 )

Data = generate_gauss_fdata( N,
                             centerline = sin( 4 * pi * time_grid ),
                             Cov = Cov )


Data_out = array( 0, dim = c( N_outliers, P ) )


Data_out[ 1, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * time_grid + pi / 2 ),
                                       Cov = Cov )

Data_out[ 2, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * time_grid - pi / 2 ),
                                       Cov = Cov )

Data_out[ 3, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * time_grid + pi/ 3 ),
                                       Cov = Cov )

Data_out[ 4, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * time_grid - pi / 3),
                                       Cov = Cov )
Data = rbind( Data, Data_out )

test_that( 'Outliergram - no adjustment',
           expect_equal( outliergram( time_grid, Data, display = FALSE ),
                         c( 3, 20, 30, 31, 101, 112, 175, 199, 201, 202, 203, 204 ) ) )

test_that( 'Outliergram - with adjustment',
           expect_equal( outliergram( time_grid, Data,
                                      adjust = list( N_trials = 10,
                                                     trial_size = 5 * nrow( Data ),
                                                     TPR = 0.01,
                                                     VERBOSE = FALSE ),
                                      display = FALSE ),
                         c( 30,  31, 101, 175, 201, 202, 203, 204 ) ) )

