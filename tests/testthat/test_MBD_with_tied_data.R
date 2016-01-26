
# TESTING MBDs IN PRESENCE OF TIED DATA -----------------------------------

D = matrix( c( c( 1, 0.5, 0.25, 0.1, 0.05   ),
               c( 1, 0.75, 0.25, 0.2, 0.1   ),
               c( 1, 0.7, 0.25, 0.25, 0.15  ),
               c( 1, 0.9, 0.35, 0.3, 0.25   ),
               c( 1, 0.6, 0.25, 0.2, 0.2    ) ,
               c( 0.9, 0.8, 0.25, 0.1, 0.08 ),
               c( 1, 0.4, 0.3, 0.2, 0.1     ),
               c( 1, 0.4, 0.3, 0.2, 0.1 )   ),
            ncol = 5, nrow = 8, byrow = T )

# quartz()
# matplot( t( D ), lty = 1, type = 'l' )

# Direct computation of MBDs
N = nrow( D )
P = ncol( D )

depths = rep( 0, N )

for( i in 1 : N )
{
  for( j in 1 : ( N - 1 ) )
  {
    for( k in ( j + 1 ) : N )
    {
      for( r in 1 : P )
      {
        if( ( D[ j, r ] - D[ i, r ] ) * ( D[ k, r ] - D[ i, r ] ) <= 0  )
        {
          depths[ i ] = depths[ i ] + 1
        }
      }
    }
  }
}
depths = depths / ( N * ( N - 1 ) / 2 * P )


test_that( "Correct behaviour of MBD in presence if tied-data",
           expect_equal( depths, MBD( D ) ) )


