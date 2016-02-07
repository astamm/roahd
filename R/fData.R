
#'
#' \code{fData} generic class for (univariate) functional data.
#'
#' @param grid the (evenly spaced) grid over which the functional data is defined
#' @param values the time-by-time values of the functional dataset, provided either in a vector or matrix-like representation, with N rows (observations) and P columns (time points)
#'
fData = function( grid, values )
{
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in fData: you provided an unevenly spaced grid')

  h = grid[ 2 ] - grid[ 1 ]

  P = length( grid )

  values = toRowMatrixForm( values )

  return( structure( list( t0 = grid[ 1 ],
                           tP = grid[ P ],
                           h = h,
                           P = P,
                           N = nrow( values ),
                           values = values ),
                     class = c( 'fData' ) ) )
}


plot.fData = function( fData, lty = 1, col = NULL, ... )
{
  if( is.null( col ) )
  {
    col = set_alpha( scales::hue_pal( )( fData$N ),
                     0.8 )
  }

  matplot( seq( fData$t0, fData$tP, length.out = fData$P ),
           t( fData$values ), lty = lty, type = 'l',
           col =  col, ...  )

}

#'
#' \code{mfData} generic class for (multivariate) functional data.
#'
#' @param grid the (evenly spaced) grid over which the functional data is defined
#' @param Data_list a list containing the time-by-time values of functional dataset, where each node of the list (possibly named) represents a particular dimension of the multivariate dataset.
#'
mfData = function( grid, Data_list )
{
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in mfData: you provided an unevenly spaced grid')

  dimMatrix = sapply( Data_list, dim )

  if( unique( apply( dimMatrix, 1, function( x )( length( unique( x ) ) ) ) ) != 1 )
  {
    stop( ' Error in mfData: you provided mismatching datasets as Data_list')
  }

  L = length( Data_list )

  if( is.null( names( Data_list ) ) )
  {
    nms = 1 : L
  } else {
    nms = names( Data_list )
  }

  fDList = NULL
  message( ' * * * Think about this representation' )
  for( nmL in nms )
  {
    fDList[[ nmL ]] = fData( grid, toRowMatrixForm( Data_list[[ nmL ]] ) )
  }

  return( structure( list( L = L,
                           N = fDList[[ 1 ]]$N,
                           P = fDList[[ 1 ]]$P,
                           fDList = fDList ),
                     class = c( 'mfData' ) ) )
}


plot.mfData = function( mfData, lty = 1, col = NULL, ... )
{
  if( is.null( col ) )
  {
    col = set_alpha( scales::hue_pal( )( mfData$N ),
                     0.8 )
  }

  mfrow_rows = ceiling( mfData$L / 2 )
  mfrow_cols = 2

  par( mfrow = c( mfrow_rows, mfrow_cols ) )

  invisible( sapply( mfData$fDList, plot, lty, col, ... ) )
}
