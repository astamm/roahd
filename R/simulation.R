#' Functional boxplot
#'
#' \code{fbplot} displays the functional boxplot of a dataset of functional data.
#'
#' @param M the number of functional observations to generate
#' @param centerline the centerline of the distribution (either mean or median)
#' @param Cov covairance operator, in form of a P x P matrix (where P is the number of time points of the discrete grid over which functional data are observed)
#' @param CholCov the Cholesky factor of the P x P Covariance matrix
#'
generate_gauss_fdata = function( M, centerline, Cov = NULL, CholCov = NULL )
{
  centerline = as.vector( centerline )

  if( is.null( Cov ) & is.null( CholCov ) ){
    stop( 'Error: You have to provide at least either covariance matrix or
          its cholesky factor to .generate_gauss_fdata\n')
  } else if( ! is.null( CholCov ) ) {

    P = ncol( CholCov )

    if( length( centerline ) != nrow( CholCov ) | nrow( CholCov ) != P  ){
      stop( 'Error: You provided mismatching centerline and covaraince matrix
            Cholesky factor to .generate_gauss_fdata\n')
    }
    } else if( ! is.null( Cov ) ){

      P = ncol( Cov )

      if( length( centerline ) != nrow( Cov ) | nrow( Cov ) != P  )
      {
        stop( 'Error: You provided mismatching centerline and covaraince matrix to
              .generate_gauss_fdata\n')
      }

      CholCov = chol( Cov )
      }


  temp  = matrix( rnorm( M * P ),
                  nrow = M,
                  ncol = P ) %*% CholCov + centerline

  return( t( t( matrix( rnorm( M * P ),
                  nrow = M,
                  ncol = P ) %*% CholCov ) + centerline ) )
}
