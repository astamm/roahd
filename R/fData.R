
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


# Constructor of the class "EvenTimeStr"
mfData = function( time_grid, h = NULL, P = NULL )
{
  ! missing( grid ) || stop( ' Error in EvenTimeStr: missing grid in EvenTimeStr' )

  message( ' * * * Handle better the global constant here ' )
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in EvenTimeStr: uneven grid passed in EvenTimeStr')



  structure( list( grid = grid,
                   h = h,
                   P = P,
                   t0 = grid[ 1 ],
                   tP = grid[ P ]),
             class = c( 'EvenTimeStr', 'TimeStr' ) )
}
