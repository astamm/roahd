#' Epigraph Index of functions in a functional dataset
#'
#' \code{EI} computes the Epigraph Index of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
EI = function( Data )
{
  # Number of observations
  N = nrow( Data )

  # Number of time points
  P = ncol( Data )

  rk = apply( Data, 2, function( v )( rank( v, ties.method = 'min' ) ) )

  N_a = N - apply( rk, 1, max ) + 1

  EI = 1 - N_a / N

  return( EI )
}


#' Modified Epigraph Index of functions in a functional dataset
#'
#' \code{MEI} computes the Modified Epigraph Index of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
MEI = function( Data )
{
  # Number of observations
  N = nrow( Data )

  # Number of time points
  P = ncol( Data )

  # Matrix of ranks, with `min' policy for breaking ties
  rk = apply( Data, 2, function( v ) ( rank( v, ties.method = 'min' ) ) )

  # Number of curves equal or above, time by time
  N_a = N - rk + 1

  MEI = 1 - rowSums( N_a ) / ( N * P )

  return( MEI )
}



#' Hypograph Index of functions in a functional dataset
#'
#' \code{EI} computes the Modified Index of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
HI = function( Data )
{
  # Number of observations
  N = nrow( Data )

  rk = apply( Data, 2, function( v )( rank( v, ties.method = 'max' ) ) )

  N_b = apply( rk, 1, min )

  HI = N_b / N

  return( HI )
}


#' Modified Hypograph Index of functions in a functional dataset
#'
#' \code{MEI} computes the Modified Hypograph Index of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
MHI = function( Data )
{
  # Number of observations
  N = nrow( Data )

  # Number of time points
  P = ncol( Data )

  # Matrix of ranks, with `max' policy for breaking ties
  rk = apply( Data, 2, function( v ) ( rank( v, ties.method = 'max' ) ) )

  # Number of curves equal or below, time by time
  # N_b = rk

  MHI = rowSums( rk ) / ( N * P )

  return( MHI )
}


#' Half-region depth
#'
#' \code{HRD} computes the Half Region Depth of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
HRD = function( Data )
{
  ei = EI( Data )

  hi = HI( Data )

  return( mapply( min, 1 - ei, hi ) )
}


#' Modified Half-region depth
#'
#' \code{MHRD} computes the Modified Half Region Depth of a functional dataset
#'
#' @param Data a matrix-like dataset of functional data, with observations as rows and time points as columns
MHRD = function( Data )
{
  mei = MEI( Data )

  mhi = MHI( Data )

  return( mapply( min, 1 - mei, mhi ) )

}
