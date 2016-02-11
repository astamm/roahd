

# TESTING THE FUNCTIONAL BOXPLOT ------------------------------------------

time_grid = seq( 0, 1, length.out = 1e2 )

D = matrix( c( sin( 2 * pi * time_grid ) + 0,
               sin( 2 * pi * time_grid ) + 1,
               sin( 2 * pi * time_grid ) + 2,
               sin( 2 * pi * time_grid ) + 3,
               sin( 2 * pi * time_grid ) + 4,
               sin( 2 * pi * time_grid ) + 5,
               sin( 2 * pi * time_grid ) + 6,
               sin( 2 * pi * time_grid ) + 7,
               sin( 2 * pi * time_grid ) + 8,
               sin( 2 * pi * time_grid ) + 9,
               sin( 2 * pi * time_grid ) + 10,
               sin( 2 * pi * time_grid ) - 1,
               sin( 2 * pi * time_grid ) - 2,
               sin( 2 * pi * time_grid ) - 3,
               sin( 2 * pi * time_grid ) - 4,
               sin( 2 * pi * time_grid ) - 5,
               sin( 2 * pi * time_grid ) - 6,
               sin( 2 * pi * time_grid ) - 7,
               sin( 2 * pi * time_grid ) - 8,
               sin( 2 * pi * time_grid ) - 9,
               sin( 2 * pi * time_grid ) - 10),
            nrow = 21, ncol = length( time_grid ), byrow = T )

fD = fData( time_grid, D )

test_that( 'Functional boxplot - simple',
           expect_silent( fbplot( fD, xlab = 'time', ylab = 'value',
                                  main = 'My Functional Boxplot' ) ) )

# TESTING THE ADJUSTED FUNCTIONAL BOXPLOT ---------------------------------

time_grid = seq( 0, 1, length.out = 1e2 )

N = 1e2

Data = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ),
                             Cov = exp_cov_function( time_grid,
                                                     alpha = 0.3,
                                                     beta  = 0.4 ) )
fD = fData( time_grid, Data )

test_that( 'Functional boxplot - adjusted',
           expect_silent(
             fbplot( fD, adjust = list( N_trials = 10,
                                        trial_size =  10 * N,
                                        VERBOSE = FALSE ),
                     xlab = 'time', ylab = 'Values',
                     main = 'My adjusted functional boxplot' ) ) )
