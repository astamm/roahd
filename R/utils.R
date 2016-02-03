

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
