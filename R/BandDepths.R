
MBD = function( Data, strict_inclusion = FALSE, ties_method = 'average' )
{
  # Data:                 Is the dataset containing signals whose depths are to be computed.
  #
  # strict_inclusion:  Specifies whether in the computation of the bands one should also
  #                    count those bands created by the function itself and all the other
  #                    elements of the dataset (i.e. if one wants to adhere to Genton's
  #                    formula strictly).
  #
  # ties_method:        A flag specifying the way you want to break possible ties; the value
  #                     is just passed to `rank` hence must be an admissible value for that
  #                     function.

  # Observations
  N = nrow( Data )

  # Time points
  P = ncol( Data )

  # Vector of final depths
  Depths = rep( 0, N );

  # Actual computation of modified band depths
  if( ! strict_inclusion )
  {
    for( iTimePoint in 1 : P )
    {
      rk = rank( Data[ , iTimePoint ], ties.method = ties_method );

      for( jObs in 1 : N )
      {
        N.up = rk[ jObs ] - 1;

        N.down = N - rk[ jObs ];

        Depths[ jObs ] = Depths[ jObs ] + N.up * N.down + N - 1;
        # Depths are increased by N-1, accounting for the bands given by the current
        # signal and any other signal in the dataset.
        # In the case of vertically shiftedß signals the most outlying observations
        # would have a non-zero depth, and precisely (N-1)/(N choose 2) ).
      }
    }
  }
  else
  {
    for( iTimePoint in 1 : P )
    {
      rk = rank( Data[, iTimePoint ] );

      for( jObs in 1 : N )
      {
        N.up = rk[ jObs ] - 1;

        N.down = N - rk[ jObs ];

        Depths[ jObs ] = Depths[ jObs ] + N.up * N.down;
      }
    }
  }

  Depths = Depths / ( N * ( N - 1 ) / 2 * P );

  return( Depths )
}

MBD_relative = function( Data_target, Data_reference, strict_inclusion = FALSE, ties_method = 'average' )
{
  # Data_target:       Is the dataset containing signals whose depths are to be computed
  #                    with respect to the reference population.
  #
  # Data_reference:    Is the dataset containing the reference population used to measure
  #                    the centrality through depths. In particular, the depths of these
  #                    functions will not be computed.
  #
  # strict_inclusion:  Specifies whether in the computation of the bands one should also
  #                    count those bands created by the function itself and all the other
  #                    elements of the dataset (i.e. if one wants to adhere to Genton's
  #                    formula strictly).
  #
  # ties_method:        A flag specifying the way you want to break possible ties; the value
  #                     is just passed to `rank` hence must be an admissible value for that
  #                     function.

  # Observations
  N = nrow( Data_target )

  # Time points
  P = ncol( Data_reference )

  # Manipulation of Data_target to obtain a matrix representation such that,
  # in case of just one observation we have a matrix like structure with time
  # on the columns and observations on the rows
  if( is.null( dim( Data_target ) )|
      is.array( Data_target ) &
      length( dim( Data_target ) ) == 1 )
  {
    Data_target = t( as.matrix( Data_target ) )

  } else if( is.matrix( Data_target ) ) {

    if( ncol( Data_target ) == 1 ){
      Data_target = t( Data_target )
    }
  } else {
    stop( 'Error: unsupported value provided to MBD_relative' )
  }

  if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with not
                                       compliant dimensions to MBD_relative')

  # Number of target_functions
  N_target = nrow( Data_target )

  Depths = rep( 0, N_target );

  if( ! strict_inclusion )
  {
    # Here there can be room for some efficienty improvement
    for( iTimePoint in 1 : P )
    {
      for( jObs in 1 : N_target )
      {
        m.leq = sum( Data_reference[ , iTimePoint ] <= Data_target[ jObs, iTimePoint ] )
        m.geq = sum( Data_reference[ , iTimePoint ] >= Data_target[ jObs, iTimePoint ] )

        exceed = m.leq + m.geq   - N

        if( m.leq < m.geq ){

          m.geq = m.geq - exceed
        }
        else{

          m.leq = m.leq - exceed
        }

        Depths[ jObs ] = Depths[ jObs ] + m.leq * m.geq;
      }
    }
  }
  else
  {
    for( iTimePoint in 1 : P )
    {
      for( jObs in 1 : N_target )
      {
        m.leq = sum( Data_reference[ , iTimePoint ] < Data_target[ jObs, iTimePoint ] )
        m.geq = sum( Data_reference[ , iTimePoint ] > Data_target[ jObs, iTimePoint ] )

        Depths[ jObs ] = Depths[ jObs ] + m.leq * m.geq;
      }
    }
  }

  Depths = Depths / ( N * ( N - 1 ) / 2 * P  );

  return( Depths )
}
