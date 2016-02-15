#' Band Depth for univariate functional data
#'
#' \code{BD} computes the Band Depth of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
BD = function( Data )
{
  # Number of rows
  N = nrow( Data )

  # Compute ranks of matrix-like representation of data with `min' tie-breaking rule
  rk = apply( Data, 2, function( v )( rank( v, ties.method = 'average' ) ) )

  # Actual number of function values strictly above
  N_a = N - apply( rk, 1, max )

  # Actual number of function values strictly below
  N_b = apply( rk, 1, min ) - 1

  Depths = ( N_a * N_b + ( N - 1 ) ) / ( N * ( N - 1 ) / 2 )

  return( Depths )
}

#' Relative Band Depth of functions in a functional dataset
#' \code{BD_relative} computes the Band Depth (BD) of a functional dataset (Data_target) with respect to another (Data_reference)
#'
#' @param Data_target is the dataset containing signals whose depths are to be computed with respect to the reference population.
#' @param Data_reference is the dataset containing the reference population used to measure the centrality through depths. In particular, the depths of these functions will not be computed.
#'
BD_relative = function( Data_target, Data_reference )
{
  # Observations
  N = nrow( Data_reference )

  # Number of columns
  P = ncol( Data_reference )

  Data_target = toRowMatrixForm( Data_target )

  if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with not
                                       compliant dimensions to BD_relative')

  # Number of target_functions
  N_target = nrow( Data_target )

  Depths = rep( 0, N_target );

  for( iObs in 1 : N_target )
  {
    rk = apply( rbind( Data_target[ iObs, ], Data_reference ),
                2, rank, ties.method = 'average' )

    N_a = N + 1 - max( rk[ 1, ] )

    N_b = min( rk[ 1, ] ) - 1

    Depths[ iObs ] = ( N_a * N_b ) / ( N * ( N - 1 ) / 2 )
  }

  return( Depths )
}


#' Modified Band Depth for univariate funcitonal data
#'
#' \code{MBD} Computes the Modified Band Depth (MBD) of elements of a functional dataset
#'
#' @param Data is the dataset of functional observations whose MBD are to be computed. Its rows
#' should be the functional observations, while the columns must be time-point evaluations of such
#' observations
#' @param manage_ties is a flag specifying if the method should accomodate for the presence of ties
#' in the dataset.
MBD = function( Data, manage_ties = FALSE )
{
  # Number of rows
  N = nrow( Data )

  # Number of time points in the discrete grid
  P = ncol( Data )

  if( manage_ties )
  {
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

    # if( any( Repetitions > 1 ) )
    # {
      # Now, owing to repetitions we gain more bands containing each functional observation.
      # In standard cases (no coincident values) we have N - 1, in general we have N - 1 +
      # N - 2 + ... + N - K, where K is the number of repetitions of the function value in the
      # dataset, time by time. If you split the sum and use Gauss formula, you
      # can get the following expression

      added_bands = N * Repetitions - 0.5 * ( Repetitions * Repetitions + Repetitions )

      Depths = rowSums( N_a * N_b  + added_bands ) / ( P * ( N - 1 ) * N / 2 )

    } else {

      rk =  apply( Data, 2, function( v )( rank( v ) )  )

      Depths = ( rowSums( ( N - rk ) * ( rk - 1 ) ) / P + N - 1 ) / ( N * ( N - 1 ) / 2 )

    }

  return( Depths )
}

#' Relative Modified Band Depth of functions in a functional dataset
#' \code{MBD_relative} computes the Modified Band Depth (MBD) of a functional dataset (Data_target) with respect to another (Data_reference)
#'
#' @param Data_target is the dataset containing signals whose depths are to be computed with respect to the reference population.
#' @param Data_reference is the dataset containing the reference population used to measure the centrality through depths. In particular, the depths of these functions will not be computed.
MBD_relative = function( Data_target, Data_reference )
{
  # Observations
  N = nrow( Data_reference )

  # Number of columns
  P = ncol( Data_reference )

  Data_target = toRowMatrixForm( Data_target )

  if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with not
                                       compliant dimensions to MBD_relative')

  # Number of target_functions
  N_target = nrow( Data_target )

  Depths = rep( 0, N_target );

  for( iObs in 1 : N_target )
  {
    rk = apply( rbind( Data_target[ iObs, ], Data_reference ),
                2, rank, ties.method = 'average' );

    N_a = N + 1 - rk[ 1, ]

    N_b = rk[ 1, ] - 1

    Depths[ iObs ] = sum( N_a * N_b ) / (  P *  N * ( N - 1 ) / 2 )
  }

  return( Depths )
}
