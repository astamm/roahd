
# TESTING MAX-MIN FUNCTIONS -----------------------------------------------

P = 1e4

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
                             c( 1, 1, 0 + 4999 * h ) ) )

test_that( 'Max function for functional data, which = TRUE, value',
           expect_equal( max( fD, which = TRUE )$value,
                             c( 1, 2, 3 * ( 0.5 - abs( 0.5 - 4999 * h) ) ) ) )

test_that( 'Min function for functional data, which = TRUE, grid',
           expect_equal( min( fD, which = TRUE )$grid,
                         c( 0, 0, 0 ) ) )

test_that( 'Min function for functional data, which = TRUE, value',
           expect_equal( min( fD, which = TRUE )$value,
                         c( 0, 0, 0 ) ) )

test_that( 'Max function for functional data, which = FALSE',
           expect_equal( max( fD, which = FALSE ),
                         c( 1, 2, 3 * ( 0.5 - abs( 0.5 - 4999 * h) ) ) ) )

test_that( 'Min function for functional data, which = FALSE',
           expect_equal( min( fD, which = FALSE ),
                         c( 0, 0, 0 ) ) )

test_that( 'Area under the curve',
           expect_equal( area_under_curve( fD ),
                         c( 0.5, 1, 0.75 ) ) )

# ORDERING ----------------------------------------------------------------

P = 1e3

time_grid = seq( 0, 1, length.out = P )

h = time_grid[ 2 ] - time_grid[ 1 ]

Data_1 = matrix( c( 1 * time_grid,
                    2 *  time_grid ),
                 nrow = 2, ncol = P, byrow = TRUE )

Data_2 = matrix( 3 * ( 0.5 - abs( time_grid - 0.5 ) ),
                 nrow = 1, byrow = TRUE )

Data_3 = rbind( Data_1, Data_1 )


fD_1 = fData( time_grid, Data_1 )
fD_2 = fData( time_grid, Data_2 )
fD_3 = fData( time_grid, Data_3 )

# MAX ORDERING
test_that( 'Max_ordering - case 1',
           expect_equal( max_ordered( fD_1, fD_2 ),
                         c( TRUE, FALSE ) ) )

test_that( 'Max_ordering - case 2',
           expect_equal( max_ordered( fD_2, fD_1 ),
                         c( FALSE, TRUE ) ) )

test_that( 'Max_ordering - case 3',
           expect_error( max_ordered( fD_1, fD_3 ) ) )

test_that( 'Max_ordering - case 4',
           expect_error( max_ordered( fD_3, fD_1 ) ) )

test_that( 'Max_ordering - case 5',
           expect_equal( max_ordered( fD_2, fD_3 ),
                         c( FALSE, TRUE, FALSE, TRUE ) ) )

test_that( 'Max_ordering - case 6',
           expect_equal( max_ordered( fD_3, fD_2 ),
                         c( TRUE, FALSE, TRUE, FALSE ) ) )

# AREA ORDERING
test_that( 'Area ordering - case 1',
           expect_equal( area_ordered( fD_1, fD_2 ),
                         c( TRUE, FALSE ) ) )

test_that( 'Area ordering - case 2',
           expect_equal( area_ordered( fD_2, fD_1 ),
                         c( FALSE, TRUE ) ) )

test_that( 'Max_ordering - case 3',
           expect_error( area_ordered( fD_1, fD_3 ) ) )

test_that( 'Max_ordering - case 4',
           expect_error( area_ordered( fD_3, fD_1 ) ) )

test_that( 'Max_ordering - case 5',
           expect_equal( area_ordered( fD_2, fD_3 ),
                         c( FALSE, TRUE, FALSE, TRUE ) ) )

test_that( 'Max_ordering - case 6',
           expect_equal( max_ordered( fD_3, fD_2 ),
                         c( TRUE, FALSE, TRUE, FALSE ) ) )

