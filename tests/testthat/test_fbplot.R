

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

fbplot( time_grid, Data = D,  )


# TESTING THE ADJUSTED FUNCTIONAL BOXPLOT ---------------------------------

time_grid = seq( 0, 1, length.out = 1e2 )

N = 1e2

# C( s, t ) = \alpha \exp( - beta | s - t | )
# Amplitude factor
alpha = 0.3

# Correlation decay factor
beta = 0.4

Cov = outer( time_grid, time_grid, function( s, t )( alpha * exp( - beta * abs( s - t ) ) ) )

Data = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )

fbplot( time_grid, Data, adjust = list( N_trials = 2,
                                        VERBOSE = TRUE ),
        xlab = 'time [ms]', ylab = 'Data', main = 'Functional boxplot' )
