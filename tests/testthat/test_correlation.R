
# TESTING MAX-MIN FUNCTIONS -----------------------------------------------

P = 1e3

time_grid = seq( 0, 1, length.out = P )

h = time_grid[ 2 ] - time_grid[ 1 ]

Data = matrix( c( 1 * time_grid,
                  2 *  time_grid,
                  3 * ( 0.5 - abs( time_grid - 0.5 ) ) ),
               nrow = 3, ncol = P, byrow = TRUE )

fD = fData( time_grid, Data )

# quartz()
# plot( fD, lwd = 2 )

test_that( 'Max function for functional data, which = TRUE, grid',
           expect_equal( max( fD, which = TRUE )$grid,
                             c( 1, 1, 0 + 499 * h ) ) )

test_that( 'Max function for functional data, which = TRUE, value',
           expect_equal( max( fD, which = TRUE )$value,
                             c( 1, 2, 3 * ( 0.5 - abs( 0.5 - 499 * h) ) ) ) )

test_that( 'Min function for functional data, which = TRUE, grid',
           expect_equal( min( fD, which = TRUE )$grid,
                         c( 0, 0, 0 ) ) )

test_that( 'Min function for functional data, which = TRUE, value',
           expect_equal( min( fD, which = TRUE )$value,
                         c( 0, 0, 0 ) ) )

test_that( 'Max function for functional data, which = FALSE',
           expect_equal( max( fD, which = FALSE ),
                         c( 1, 2, 3 * ( 0.5 - abs( 0.5 - 499 * h) ) ) ) )

test_that( 'Min function for functional data, which = FALSE',
           expect_equal( min( fD, which = FALSE ),
                         c( 0, 0, 0 ) ) )
