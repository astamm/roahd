#' Functional boxplot
#'
#' \code{generate_gauss_data} displays the functional boxplot of a dataset of
#' functional data.
#'
#' @param M the number of functional observations to generate
#' @param centerline the centerline of the distribution (either mean or median)
#' @param Cov covairance operator, in form of a P x P matrix (where P is the
#' number of time points of the discrete grid over which functional data are
#' observed)
#' @param CholCov the Cholesky factor of the P x P Covariance matrix
#'
generate_gauss_fdata = function( M, centerline = NULL,
                                 Cov = NULL, CholCov = NULL )
{
  if( is.null( Cov ) & is.null( CholCov ) ){
    stop( 'Error: You have to provide at least either covariance matrix or
          its cholesky factor to .generate_gauss_fdata\n')
  } else if( ! is.null( CholCov ) ) {

    P = ncol( CholCov )

    if( length( centerline ) != nrow( CholCov ) | nrow( CholCov ) != P  ){
      stop( 'Error: You provided mismatching centerline and covaraince matrix
            Cholesky factor to generate_gauss_fdata\n')
    }
    } else if( ! is.null( Cov ) ){

      P = ncol( Cov )

      if( length( centerline ) != nrow( Cov ) | nrow( Cov ) != P  )
      {
        stop( 'Error: You provided mismatching centerline and covariance matrix
to generate_gauss_fdata\n')
      }

      CholCov = chol( Cov )
      }

  return( t( t( matrix( rnorm( M * P ),
                  nrow = M,
                  ncol = P ) %*% CholCov ) + centerline ) )
}

generate_gauss_mfdata = function( M, L, centerline, correlations,
                                  listCov = NULL, listCholCov = NULL )
{
  if( length( correlations ) != 0.5 * ( L ) * ( L - 1 ) )
  {
    stop( 'Error in generate_gauss_mfdata: you have to provide all the
          correlations among functional components' )
  }

  if( nrow( centerline ) != L )
  {
    stop( 'Error in generate_gauss_mfdata: you have to provide a centerline for
each dimension' )
  }

  if( is.null( listCov ) & is.null( listCholCov ) ){
    stop( 'Error: You have to provide at least either covariance matrices or
          their cholesky factors to generate_gauss_mfdata')
  } else if( ! is.null( listCholCov ) )
  {
      if( length( listCholCov ) != L )
      {
        stop( 'Error: You have to provide a covariance Cholesky factor for each
              dimension' )
      }

      P = ncol( listCholCov[[ 1 ]] )

      if( ncol( centerline ) != P | any( sapply( listCholCov, nrow ) != P ) |
          any( sapply( listCholCov, ncol ) != P ) )
      {
        stop( 'Error: You provided mismatching centerline and covariance
matrices Cholesky factors to generate_gauss_mfdata')
      }
    } else if( ! is.null( listCov ) )
    {

      P = ncol( listCov[[ 1 ]] )

        if( ncol( centerline ) != P | any( sapply( listCov, nrow ) != P ) |
            any( sapply( listCov, ncol ) != P ) )
        {
          stop( 'Error: You provided mismatching centerline and covariance
matrices to generate_gauss_mfdata')
        }
      listCholCov = lapply( listCov, chol )
    }

  # Generating the matrix of correlations among dimensions
  R = matrix( 1, ncol = L, nrow = L )

  R[ upper.tri( R ) ] = as.numeric( correlations )

  R[ lower.tri( R ) ] = as.numeric( correlations )

  R_chol = chol( R )

  # Generating gaussian data with correlations among dimensions
  Data = matrix( rnorm( M * L * P ), ncol = L, nrow = M * P )

  Data = Data %*% R_chol

  return( eval( parse( text =
                         paste( 'list( ',
                                paste( 't( t( matrix( Data[ , ', 1 : L,
                                       ' ], nrow = M, ncol = P ) %*% listCholCov[[ ',
                                       1 : L, ' ]] ) + as.numeric( centerline[ ',
                                       1 : L, ', ] ) )',
                                       sep = '', collapse  = ', ' ),
                                ' )', sep = '' ) ) ) )
}
