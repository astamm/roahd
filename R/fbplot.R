#' Functional boxplot
#'
#' \code{fbplot} displays the functional boxplot of a dataset of functional data.
#'
#' @param time_grid optional, a vector corresponding to the time grid over which the functions are evaluated (it must be
#' evenly spaced)
#' @param Data the dataset of functional data, having observations on the rows and time points on columns
#' @param Depths either a vector containing the depths for each element of the dataset, or a string containing the name of the method you want to use to compute it. In this case the name of the method must be included in the caller's environment
#' @param adjust either FALSE if you would like the default value for the inflation factor, F = 1.5, to be used, or a list specifying the parameters required by the adjustment.
#' @param display either a logical value indicating wether you want the functional boxplot to be displayed, or the number of the graphical devices where you want the boxplot to be plotted.
#'
fbplot = function( time_grid = NULL, Data, Depths = 'MBD',
                   adjust = FALSE, display = TRUE, ..., verbose = FALSE )
{
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
    Depths_spec = Depths
    Depths = eval( parse( text = paste( Depths, '( Data )', sep = '' ) ) )
  } else {
    stopifnot( length( Depths ) == N )
  }

  if( ! is.list( adjust ) )
  {
    # Plain functional boxplot with default F value: F = 1.5

    out = .fbplot( time_grid, Data, Depths, Fvalue = 1.5 )

  } else {

    library( robustbase )

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         5 * N,
                         adjust$trial_size )

    TPR = ifelse( is.null( adjust$TPR ),
                  2 * pnorm( 4 * qnorm( 0.25 ) ),
                  adjust$FPR )

    F_min = ifelse( is.null( adjust$F_min ),
                    0.5,
                    adjust$F_min )

    F_max= ifelse( is.null( adjust$F_max ),
                   5,
                   adjust$F_max )

    tol = ifelse( is.null( adjust$tol ),
                  1e-3,
                  adjust$tol )

    maxiter = ifelse( is.null( adjust$maxiter ),
                      100,
                      adjust$maxiter )

    VERBOSE = ifelse( is.null( adjust$VERBOSE ),
                      FALSE,
                      adjust$VERBOSE )

    # Estimation of robust covaraince matrix
    Cov = covOGK( Data, sigmamu = s_Qn )$cov

    # Cholesky factor
    CholCov <- chol( Cov )

    # Centerline of the dataset
    centerline = Data[ which.max( Depths ), ]

    Fvalues = rep( 0, N_trials )

    cost_functional = function( F_curr )( length( .fbplot( time_grid,
                                                           Data_gauss,
                                                           Depths = Depths_spec,
                                                           Fvalue = F_curr )$ID_out ) / trial_size - TPR )

    for( iTrial in 1 : N_trials )
    {
      if( VERBOSE > 0 )
      {
        cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
      }

      Data_gauss = generate_gauss_fdata( trial_size, centerline, CholCov = CholCov )

      if( VERBOSE > 0 )
      {
        cat( ' * * * * beginning optimisation\n' )
      }

      opt = uniroot( cost_functional,
                     interval = c( F_min, F_max ),
                     tol = tol,
                     maxiter = maxiter )
      if( VERBOSE > 0 )
      {
        cat( ' * * * * optimisation finished.\n')
      }

      Fvalues[ iTrial ] = opt$root
    }

    Fvalue = mean( Fvalues )

    out = .fbplot( time_grid, Data, Depths, Fvalue = Fvalue  )
  }

  ID_out = out$ID_out

  # Plotting part
  if( is.numeric( display ) )
  {
    dev.set( display )
  }

  if( ! display == FALSE )
  {
    # Creating color palettes
    colors.default = colorRampPalette( RColorBrewer::brewer.pal(9,"Blues") )( nrow( Data ) - length( ID_out ) )
    colors.out     = colorRampPalette( RColorBrewer::brewer.pal(9,"Reds")[5:9] )( length( ID_out ) )

    # Plotting non-outlying data
    matplot( time_grid, t( Data[ - ID_out, ] ), lty = 1, type = 'l',
             col = colors.default, ylim = range( Data ), ... )

    # Plotting outlying data
    matplot( time_grid, t( Data[ ID_out, ] ), lty = 1, type = 'l', col = colors.out, lwd = 1, add = T )

    # Highligting the central envelope
    lines( time_grid, out$max_envelope_central, lty = 1, col = 'yellow', lwd = 3 )
    lines( time_grid, out$min_envelope_central, lty = 1, col = 'yellow', lwd = 3 )

    # Filling in the central envelope
    rgb.temp = col2rgb( 'yellow' )
    polygon( c(time_grid, rev( time_grid) ), c( out$min_envelope_central, rev( out$max_envelope_central ) ),
             col = rgb( rgb.temp[1], rgb.temp[2], rgb.temp[3], alpha = 100, maxColorValue = 255 ), border = NA)

    # Plotting the sample median
    lines( time_grid, Data[ which.max( Depths ), ], lty = 1, type = 'l', col = 'white', lwd = 3)

    # Plotting fences
    id_uppermost = which.min( apply( - t( t( Data[ - ID_out, ]) + out$fence_upper ), 1, min )  )
    id_lowermost = which.min( apply(   t( t( Data[ - ID_out, ]) - out$fence_lower ), 1, min )  )

    ## actual fences ( F * envelope)
    lines( time_grid, out$fence_upper, lty = 2, col = 'red', lwd = 2 )
    lines( time_grid, out$fence_lower, lty = 2, col = 'red', lwd = 2 )

    ## closest data to fences
    lines( time_grid, Data[ -ID_out, ][ id_uppermost, ], lty = 1, col = 'darkorange1', lwd = 3 )
    lines( time_grid, Data[ -ID_out, ][ id_lowermost, ], lty = 1, col = 'darkorange1', lwd = 3 )

    ## vertical segments
    half.time_grid = which.min( abs( time_grid - 0.5 ) )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$max_envelope_central[ half.time_grid ], Data[ -ID_out, half.time_grid ][ id_uppermost ] ),
           lty = 1, col = 'orange2', lwd = 3 )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$min_envelope_central[ half.time_grid ], Data[ -ID_out, half.time_grid ][ id_lowermost ] ),
           lty = 1, col = 'orange2', lwd = 3 )
  }

  return( list( Depth = Depths,
                ID_outliers = ID_out ) )
}


.fbplot = function( time_grid, Data, Depths = 'MBD', Fvalue = 1.5 )
{
  # Number of observations
  N = nrow( Data )

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

  id_central_region = which( Depths >= quantile( Depths, prob = 0.5 ) )

  max_envelope_central = apply( Data[ id_central_region, ], 2, max )
  min_envelope_central = apply( Data[ id_central_region, ], 2, min )

  fence_upper = ( max_envelope_central - Data_mean ) * Fvalue + Data_mean
  fence_lower = ( min_envelope_central - Data_mean ) * Fvalue + Data_mean

  ID_outlying = which( apply( Data, 1, function(x)( any( x > fence_upper  ) | any ( x < fence_lower ) ) ) )

  return( list( ID_out = ID_outlying,
                min_envelope_central = min_envelope_central,
                max_envelope_central = max_envelope_central,
                fence_lower = fence_lower,
                fence_upper = fence_upper ) )
}
