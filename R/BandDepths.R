

# BD_relative = function( Data_target, Data_reference, strict_inclusion = FALSE, ties_method = 'average' )
# {
#   # Data_target:       Is the dataset containing signals whose depths are to be computed
#   #                    with respect to the reference population.
#   #
#   # Data_reference:    Is the dataset containing the reference population used to measure
#   #                    the centrality through depths. In particular, the depths of these
#   #                    functions will not be computed.
#   #
#   # strict_inclusion:  Specifies whether in the computation of the bands one should also
#   #                    count those bands created by the function itself and all the other
#   #                    elements of the dataset (i.e. if one wants to adhere to Genton's
#   #                    formula strictly).
#   #
#   # ties_method:        A flag specifying the way you want to break possible ties the value
#   #                     is just passed to `rank` hence must be an admissible value for that
#   #                     function.
#
#   # Observations
#   N = nrow( Data_target )
#
#   # Manipulation of Data_target to obtain a matrix representation such that,
#   # in case of just one observation we have a matrix like structure with time
#   # on the columns and observations on the rows
#   if( is.null( dim( Data_target ) )|
#       is.array( Data_target ) &
#       length( dim( Data_target ) ) == 1 )
#   {
#     Data_target = t( as.matrix( Data_target ) )
#
#   } else if( is.matrix( Data_target ) ) {
#
#     if( ncol( Data_target ) == 1 ){
#       Data_target = t( Data_target )
#     }
#   } else {
#     stop( 'Error: unsupported value provided to MBD_relative' )
#   }
#
#   if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with not
#                                        compliant dimensions to MBD_relative')
#
#   # Number of target_functions
#   N_target = nrow( Data_target )
#
#   Depths = rep( 0, N_target );
#
#   for( iObs in 1 : N_target )
#   {
#     rk = t( apply( cbind( Data_reference, Data_target[ iObs, ] ),
#                    2,
#                    function( v )( rank( v, ties.method = ties_method ) ) ) );
#
#     N_a = N - max( rk[ 1, ] )
#
#     N_b = min( rk[ 1, ] ) - 1
#
#     if( ! strict_inclusion )
#     {
#       Depths[ iObs ] = N_a * N_b
#     } else {
#       Depths[ i]
#     }
#
#   }
#
#
# }

BD = function( Data )
{
  # Number of rows
  N = nrow( Data )

  # Number of columns
  P = ncol( Data )

  # Compute ranks of matrix-like representation of data with `min' tie-breaking rule
  rk_min = apply( D, 2, function( v )( rank( v, ties.method = 'min' ) ) )

  # Compute ranks of matrix-like representation of data with `max' tie-breaking rule
  rk_max = apply( D, 2, function( v )( rank( v, ties.method = 'max' ) ) )

  # Actual number of function values strictly above
  N_a = N - apply( rk_max, 1, max )

  # Actual number of function values strictly below
  N_b = apply( rk_min, 1, min ) - 1

  Depths = apply( N_a * N_b + ( N - 1 ),
                  1,
                  sum ) / ( P * N * ( N - 1 ) / 2 )

  return( Depths )
}


MBD = function( Data )
{
  # Number of rows
  N = nrow( Data )

  # Number of time points in the discrete grid
  P = ncol( Data )

  # Compute ranks of matrix-like representation of data with `min' tie-breaking rule
  rk_min = apply( Data, 2, function( v )( rank( v, ties.method = 'min' ) ) )

  # Compute ranks of matrix-like representation of data with `max' tie-breaking rule
  rk_max = apply( Data, 2, function( v )( rank( v, ties.method = 'max' ) ) )

  # Times each function value is repeated in the dataset, for each time point
  # ( matrix is N x P)
  Repetitions = rk_max - rk_min + 1

  # Actual number of function values strictly above
  N_a = N - rk_max

  # Actual number of function values strictly below
  N_b = rk_min - 1

  # added_bands_function = function( N, K )
  # {
  #   sapply( K, function( k )( 0.5 * N * ( N - 1 ) - 0.5 * ( N - k ) * ( N - k - 1 )  ) )
  # }

  # Now, owing to repetitions we gain more bands containing each functional observation.
  # In standard cases (no coincident values) we have N - 1, in general we have N - 1 +
  # N - 2 + ... + N - K, where K is the number of repetitions of the function value in the
  # dataset, time by time. If you split the sum and use Gauss formula twice, you
  # can get the following expression
  added_bands_function = function( N, K )
  {
    sapply( K,
            function( k )( N * k - 0.5 * ( k * k + k  )  ) )
  }

  added_bands = t( apply( Repetitions,
                          1,
                          function( reps )( added_bands_function( N, reps ) ) ) )

  Depths = apply( N_a * N_b + added_bands,
                  1,
                  sum ) / ( P * N * ( N - 1 ) / 2 )

  return( Depths )
}

# MBD_relative = function( Data_target, Data_reference, strict_inclusion = FALSE, ties_method = 'average' )
# {
#   # Data_target:       Is the dataset containing signals whose depths are to be computed
#   #                    with respect to the reference population.
#   #
#   # Data_reference:    Is the dataset containing the reference population used to measure
#   #                    the centrality through depths. In particular, the depths of these
#   #                    functions will not be computed.
#   #
#   # strict_inclusion:  Specifies whether in the computation of the bands one should also
#   #                    count those bands created by the function itself and all the other
#   #                    elements of the dataset (i.e. if one wants to adhere to Genton's
#   #                    formula strictly).
#   #
#   # ties_method:        A flag specifying the way you want to break possible ties the value
#   #                     is just passed to `rank` hence must be an admissible value for that
#   #                     function.
#
#   # Observations
#   N = nrow( Data_target )
#
#   # Time points
#   P = ncol( Data_reference )
#
#   # Manipulation of Data_target to obtain a matrix representation such that,
#   # in case of just one observation we have a matrix like structure with time
#   # on the columns and observations on the rows
#   if( is.null( dim( Data_target ) )|
#       is.array( Data_target ) &
#       length( dim( Data_target ) ) == 1 )
#   {
#     Data_target = t( as.matrix( Data_target ) )
#
#   } else if( is.matrix( Data_target ) ) {
#
#     if( ncol( Data_target ) == 1 ){
#       Data_target = t( Data_target )
#     }
#   } else {
#     stop( 'Error: unsupported value provided to MBD_relative' )
#   }
#
#   if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with not
#                                        compliant dimensions to MBD_relative')
#
#   # Number of target_functions
#   N_target = nrow( Data_target )
#
#   Depths = rep( 0, N_target );
#
#   # Select formula depending on strict_inclusion policy
#   if( ! strict_inclusion )
#   {
#     # Here there can be room for some efficienty improvement
#
#     # First I select a time point, then proceed to compute depth of each element
#     # in the target group
#     for( iTimePoint in 1 : P )
#     {
#       # Loop over target group of functions
#       for( jObs in 1 : N_target )
#       {
#         m.leq = sum( Data_reference[ , iTimePoint ] <= Data_target[ jObs, iTimePoint ] )
#         m.geq = sum( Data_reference[ , iTimePoint ] >= Data_target[ jObs, iTimePoint ] )
#
#         exceed = m.leq + m.geq - N
#
#         if( m.leq < m.geq ){
#
#           m.geq = m.geq - exceed
#         }
#         else{
#
#           m.leq = m.leq - exceed
#         }
#
#         # Update Depths
#         Depths[ jObs ] = Depths[ jObs ] + m.leq * m.geq
#       }
#     }
#   }
#   else
#   {
#     # First I select a time point, then proceed to compute depth of each element
#     # in the target group
#     for( iTimePoint in 1 : P )
#     {
#       # Loop over target group of functions
#       for( jObs in 1 : N_target )
#       {
#         m.leq = sum( Data_reference[ , iTimePoint ] < Data_target[ jObs, iTimePoint ] )
#         m.geq = sum( Data_reference[ , iTimePoint ] > Data_target[ jObs, iTimePoint ] )
#
#         Depths[ jObs ] = Depths[ jObs ] + m.leq * m.geq
#       }
#     }
#   }
#
#   Depths = Depths / ( N * ( N - 1 ) / 2 * P  )
#
#   return( Depths )
# }
