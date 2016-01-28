#' Functional boxplot
#'
#' \code{fbplot} displays the functional boxplot of a dataset of functional data.
#'
#' @param time_grid optional, a vector corresponding to the time grid over which the functions are evaluated (it must be
#' evenly spaced)
#' @param Data the dataset of functional data, having observations on the rows and time points on columns
#' @param Depths either a vector containing the depths for each element of the dataset, or a string containing the name of the method you want to use to compute it. In this case the name of the method must be included in the caller's environment
#' @param Fvalue the inflation factor of the functional boxplot
#' @param display either a logical value indicating wether you want the functional boxplot to be displayed, or the number of the graphical devices where you want the boxplot to be plotted.
#'
fbplot = function( time_grid = NULL, Data, Depths = 'MBD', Fvalue = 1.5, display = TRUE, ... )
{
  library(RColorBrewer)

  # Number of observations
  N = nrow( Data )

  if( is.null( time_grid ) )
  {
    time_grid = 1 : ncol( Data)
  }

  stopifnot( length( time_grid ) == ncol( Data ) )

  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths = eval( parse( text = paste( Depths, '( Data )', sep = '' ) ) )
  } else {
    stopifnot( length( Depths ) == N )
  }

  Data_mean = apply( Data, 2, mean )

  id.central.region = which( Depths >= quantile( Depths, prob = 0.5 ) )

  max.envelope.central = apply( Data[ id.central.region, ], 2, max )
  min.envelope.central = apply( Data[ id.central.region, ], 2, min )

  fence.upper = ( max.envelope.central - Data_mean ) * Fvalue + Data_mean
  fence.lower = ( min.envelope.central - Data_mean ) * Fvalue + Data_mean

  id_out = which( apply( Data, 1, function(x)( any( x > fence.upper  ) | any ( x < fence.lower ) ) ) )
  # prob.out = length( id_out ) / nrow( Data )

  if( is.numeric( display ) )
  {
    dev.set( display )
  }

  if( ! display == FALSE )
  {
    # Creating color palettes
    colors.default = colorRampPalette( RColorBrewer::brewer.pal(9,"Blues") )( nrow( Data ) - length( id_out ) )
    colors.out     = colorRampPalette( RColorBrewer::brewer.pal(9,"Reds")[5:9] )( length( id_out ) )

    # Plotting non-outlying data
    matplot( time_grid, t( Data[ - id_out, ] ), lty = 1, type = 'l',
             col = colors.default, ylim = range( Data ), ... )

    # Plotting outlying data
    matplot( time_grid, t( Data[ id_out, ] ), lty = 1, type = 'l', col = colors.out, lwd = 1, add = T )

    # Highligting the central envelope
    lines( time_grid, max.envelope.central, lty = 1, col = 'yellow', lwd = 3 )
    lines( time_grid, min.envelope.central, lty = 1, col = 'yellow', lwd = 3 )

    # Filling in the central envelope
    rgb.temp = col2rgb( 'yellow' )
    polygon( c(time_grid, rev( time_grid) ), c( min.envelope.central, rev( max.envelope.central ) ),
             col = rgb( rgb.temp[1], rgb.temp[2], rgb.temp[3], alpha = 100, maxColorValue = 255 ), border = NA)

    # Plotting the sample median
    lines( time_grid, Data[ which.max( Depths ), ], lty = 1, type = 'l', col = 'white', lwd = 3)

    # Plotting fences
    id.uppermost = which.min( apply( - t( t( Data[ - id_out, ]) + fence.upper ), 1, min )  )
    id.lowermost = which.min( apply(   t( t( Data[ - id_out, ]) - fence.lower ), 1, min )  )

    ## actual fences ( F * envelope)
    lines( time_grid, fence.upper, lty = 2, col = 'red', lwd = 2 )
    lines( time_grid, fence.lower, lty = 2, col = 'red', lwd = 2 )

    ## closest data to fences
    lines( time_grid, Data[ -id_out, ][ id.uppermost, ], lty = 1, col = 'darkorange1', lwd = 3 )
    lines( time_grid, Data[ -id_out, ][ id.lowermost, ], lty = 1, col = 'darkorange1', lwd = 3 )

    ## vertical segments
    half.time_grid = which.min( abs( time_grid - 0.5 ) )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( max.envelope.central[ half.time_grid ], Data[ -id_out, half.time_grid ][ id.uppermost ] ),
           lty = 1, col = 'orange2', lwd = 3 )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( min.envelope.central[ half.time_grid ], Data[ -id_out, half.time_grid ][ id.lowermost ] ),
           lty = 1, col = 'orange2', lwd = 3 )
  }

  return( list( Depth = Depths,
                id_outliers = id_out ) )
}
