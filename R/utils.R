

#' Conversion to row-matrix form of generic vectors/arrays/matrices
#'
#' \code{toRowMatrixForm} manipulates D to obtain a matrix representation such that,
#' in case of just one observation we have a row vector representation, and in general
#' we have a consistent output with matrix-structure
#'
#' @param D a generic array or atomic vecor to be converted in row-matrix format.
toRowMatrixForm = function( D )
{
  if( is.null( dim( D ) ) |
      is.array( D ) &
      length( dim( D ) ) == 1 )
  {
    D = t( as.matrix( D ) )

  } else if( is.matrix( D ) ) {

    if( ncol( D ) == 1 ){
      D = t( D )
    }
  } else {
    stop( 'Error: unsupported value provided to toRowMatrixForm' )
  }

  return( D )
}

#' Function to setup alpha value for a set of colors
#'
#' \code{set_alpha} manipulates a vector of color representations in order to
#' setup the alpha value, and get the desired transparency level
#'
#' @param col a vector of colors
#' @param alpha the value(s) of alpha for (each of) the colors
#'
set_alpha = function( col, alpha )
{
  alpha = alpha * 255

  rgb_colors = rbind( col2rgb( col ), alpha = alpha )

  return( apply( rgb_colors, 2, function( x )( rgb( x[ 1 ],
                                                    x[ 2 ],
                                                    x[ 3 ],
                                                    x[ 4 ],
                                                    maxColorValue = 255 ) ) ) )
}

#' Function to obtain an exponential covariance function over a grid
#'
#'  \eqn{C( s, t ) = \alpha e^{ - \beta | s - t | }}
#'
#' @param time_grid a vector of time points
#' @param alpha the alpha parameter in the exponential covariance formula
#' @param beta the beta parameter in the exponential covariance formula
#'
exp_cov_function = function( time_grid, alpha, beta )
{
  return(  outer( time_grid, time_grid,
                  function( s, t )( alpha * exp( - beta * abs( s - t ) ) ) ) )
}
