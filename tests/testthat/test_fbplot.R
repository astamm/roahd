


# TESTING THE FUNCTIONAL BOXPLOT FOR UNVIARIATE FUNCTIONAL DATA -----------

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

test_that( 'Functional boxplot for univariate data - simple - 1',
           expect_silent( fbplot( fD, Fvalue = 10 ) ) )

test_that( 'Functional boxplot for univariate data - simple - 2',
           expect_silent(
             fbplot( fD, xlab = 'time', ylab = 'value',
                     main = 'My Functional Boxplot' )
           ) )


# TESTING THE ADJUSTED FUNCTIONAL BOXPLOT OF UNIVARIATE DATA --------------

set.seed( 161803 )

time_grid = seq( 0, 1, length.out = 1e2 )

N = 5e2

Data = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ),
                             Cov = exp_cov_function( time_grid,
                                                     alpha = 0.3,
                                                     beta  = 0.4 ) )
fD = fData( time_grid, Data )

test_that( 'Functional boxplot for univariate data - adjusted',
           expect_silent(
             fbplot( fD, adjust = list( N_trials = 10,
                                        trial_size = N,
                                        VERBOSE = FALSE ),
                     xlab = 'time', ylab = 'Values',
                     main = 'My adjusted functional boxplot' )
           ) )


# TESTING THE FUNCTIONAL BOXPLOT FOR MULTIVARIATE FUNCTIONAL DATA ---------

# 1 )

time_grid = seq( 0, 1, length.out = 1e2 )

D = matrix( c( sin( 2 * pi * time_grid ) - 10,
               sin( 2 * pi * time_grid ) - 9,
               sin( 2 * pi * time_grid ) - 8,
               sin( 2 * pi * time_grid ) - 7,
               sin( 2 * pi * time_grid ) - 6,
               sin( 2 * pi * time_grid ) - 5,
               sin( 2 * pi * time_grid ) - 4,
               sin( 2 * pi * time_grid ) - 3,
               sin( 2 * pi * time_grid ) - 2,
               sin( 2 * pi * time_grid ) - 1,
               sin( 2 * pi * time_grid ) + 0,
               sin( 2 * pi * time_grid ) + 1,
               sin( 2 * pi * time_grid ) + 2,
               sin( 2 * pi * time_grid ) + 3,
               sin( 2 * pi * time_grid ) + 4,
               sin( 2 * pi * time_grid ) + 5,
               sin( 2 * pi * time_grid ) + 6,
               sin( 2 * pi * time_grid ) + 7,
               sin( 2 * pi * time_grid ) + 8,
               sin( 2 * pi * time_grid ) + 9,
               sin( 2 * pi * time_grid ) + 10 ),
            nrow = 21, ncol = length( time_grid ), byrow = T )

mfD = mfData( time_grid, list( D, D * abs( 1 : 21 - 11 ) / 5 ) )

test_that( 'Functional boxplot for multivariate data - simple - 1',
           expect_silent( fbplot( mfD, Fvalue = 3 ) ) )

test_that( 'Functional boxplot for multivariate data - simple - 2',
           expect_error( fbplot( mfD, adjust = list( N_trials = 10 ) ) ) )

# 2 )


time_grid = seq( 0, 1, length.out = 1e2 )

N = 5e2

Data_1 = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ),
                             Cov = exp_cov_function( time_grid,
                                                     alpha = 0.3,
                                                     beta  = 0.4 ) )
Data_2 = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ),
                             Cov = exp_cov_function( time_grid,
                                                     alpha = 0.3,
                                                     beta = 0.4 ) )

mfD = mfData( time_grid, list( Data_1, Data_2 ) )

# test_that( 'Functional boxplot for multivariate data - simple - 3',
           # expect_silent(
             fbplot( mfD, Fvalue = 3 )
             # ) )

