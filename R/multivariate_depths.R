
#' Modified Band Depth for multivariate functional data
#'
#' \code{MBD} Computes the Band Depth (MBD) of elements of a multivariate
#' functional dataset
#'
#' @param listData specifies the multivariate functional dataset. It is a list
#' of 2-dimensional matirices containing as rows the observations and as columns
#' the measurements of the functional data over the continuous variable (e.g.
#' time)
#' @param weights either a set of weights (of the same length of \code{listData}
#' ) or the string \code{uniform} specifying that a set of uniform weights
#' (1 / L, where L is the number of dimensions of the functional dataset and
#' thus the length of \code{listData}) is to be used.
#' @param manage_ties a logical flag specifying whether the check for ties and
#' the relative treatment is to be carried out while computing the MBDs in each
#' dimension. It is directly passed to \code{MBD}.
multiMBD = function( listData, weights = 'uniform', manage_ties = FALSE )
{
  L = length( listData )

  N = nrow( listData[[ 1 ]] )
  P = ncol( listData[[ 2 ]] )

  if( ! all( sapply( listData, nrow ) - N == 0 ) |
      ! any( sapply( listData, ncol ) - P == 0 ) )
  {
    stop( ' Error in multiMBD: you provided a list with mismatching univariate
          functional datasets' )
  }

  if( is.character( weights ) )
  {
    if( weights == 'uniform' )
    {
      weights = rep( 1 / L, L )
    } else {
      stop( ' Error in multiMBD: misspecified weights characetr information')
    }
  } else {
    if( length( weights ) != L )
    {
      stop( ' Error in multiMBD: you provided more/fewer weights than
            required' )
    } else if( sum( weights ) != 1 ){
      stop( ' Error in multiMBD: sum of weights is not 1')
    }
    weights = as.numeric( weights )
  }

  return( as.numeric( sapply( listData, MBD, manage_ties ) %*% weights ) )
}

#' Band Depth for multivariate functional data
#'
#' \code{BD} Computes the Band Depth (MBD) of elements of a multivariate
#' functional dataset
#'
#' @param listData specifies the multivariate functional dataset. It is a list
#' of 2-dimensional matirices containing as rows the observations and as columns
#' the measurements of the functional data over the continuous variable (e.g.
#' time)
#' @param weights either a set of weights (of the same length of \code{listData}
#' ) or the string \code{uniform} specifying that a set of uniform weights
#' (1 / L, where L is the number of dimensions of the functional dataset and
#' thus the length of \code{listData}) is to be used.
multiBD = function( listData, weights = 'uniform' )
{
  L = length( listData )

  N = nrow( listData[[ 1 ]] )
  P = ncol( listData[[ 2 ]] )

  if( ! all( sapply( listData, nrow ) - N == 0 ) |
      ! any( sapply( listData, ncol ) - P == 0 ) )
  {
    stop( ' Error in multiBD: you provided a list with mismatching univariate
          functional datasets' )
  }

  if( is.character( weights ) )
  {
    if( weights == 'uniform' )
    {
      weights = rep( 1 / L, L )
    } else {
      stop( ' Error in multiBD: misspecified weights characetr information')
    }
  } else {
    if( length( weights ) != L )
    {
      stop( ' Error in multiBD: you provided more/fewer weights than
            required' )
    } else if( sum( weights ) != 1 ){
      stop( ' Error in multiMBD: sum of weights is not 1')
    }
    weights = as.numeric( weights )
    }

  return( as.numeric( sapply( listData, BD ) %*% weights ) )
  }
