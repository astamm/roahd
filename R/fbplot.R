#' Functional boxplot
#'
#' \code{fbplot} displays the functional boxplot of a dataset of functional data.
#'
#' @param fData the univariate functional dataset whose functional boxplot is
#' desired
#' @param Depths either a vector containing the depths for each element of the
#' dataset, or a string containing the name of the method you want to use to
#' compute it. In this case the name of the method must be included in the
#' caller's environment
#' @param adjust either FALSE if you would like the default value for the
#' inflation factor, F = 1.5, to be used, or a list specifying the parameters
#' required by the adjustment.
#' @param display either a logical value indicating wether you want the
#' functional boxplot to be displayed, or the number of the graphical devices
#' where you want the boxplot to be plotted.
#' @param ... additional graphical parameters to be used in plotting functions
#' @param verbose logical flag to enable verbosity in the adjustment process
#'
fbplot = function( fData, Depths = 'MBD',
                   adjust = FALSE, display = TRUE, ..., verbose = FALSE )
{
  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec == 'MBD' )
    {
      Depths = MBD( fData$values, manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( fData$values )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == fData$N )
  }

  if( ! is.list( adjust ) )
  {
    # Plain functional boxplot with default F value: F = 1.5

    out = .fbplot( time_grid, fData$values, Depths, Fvalue = 1.5 )

  } else {

    library( robustbase )

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         5 * fData$N,
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
    Cov = robustbase::covOGK( fData$values, sigmamu = robustbase::s_Qn )$cov

    # Cholesky factor
    CholCov <- chol( Cov )

    # Centerline of the dataset
    centerline = fData$values[ which.max( Depths ), ]

    Fvalues = rep( 0, N_trials )

    cost_functional = function( F_curr )( length( .fbplot( time_grid,
                                                           Data_gauss,
                                                           Depths = Depths_spec,
                                                           Fvalue = F_curr )$ID_out ) /
                                            trial_size - TPR )

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

    out = .fbplot( time_grid, fData$values, Depths, Fvalue = Fvalue  )
  }

  ID_out = out$ID_out

  # Plotting part
  if( is.numeric( display ) )
  {
    dev.set( display )
  }

  if( ! display == FALSE )
  {
    library(scales)

    # Creating color palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( fData$N - length( ID_out ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( ID_out ) )
    col_envelope = set_alpha( 'blue', alpha = 0.4 )
    col_center = set_alpha( 'blue', alpha = 1 )
    col_fence_structure = 'darkblue'

    time_grid = seq( fData$t0, fData$tP, length.out = fData$P )

    # Plotting non-outlying data
    matplot( time_grid,
             t( fData$values[ - ID_out, ] ), lty = 1, type = 'l',
             col = col_non_outlying, ylim = range( rbind( fData$values,
                                                          out$fence_upper,
                                                          out$fence_lower ) ), ... )

    # Plotting outlying data
    matplot( time_grid, t( fData$values[ ID_out, ] ), lty = 1, type = 'l',
             col = col_outlying, lwd = 3, add = T )


    # Filling in the central envelope

    polygon( c(time_grid, rev( time_grid) ),
             c( out$min_envelope_central, rev( out$max_envelope_central ) ),
             col = col_envelope, border = NA)
    lines( time_grid, out$max_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    lines( time_grid, out$min_envelope_central, lty = 1, col = col_envelope, lwd = 3 )

    # Plotting the sample median
    lines( time_grid, fData$values[ which.max( Depths ), ], lty = 1, type = 'l',
           col = col_center, lwd = 3)

    # Plotting fences ( F * envelope)
    # lines( time_grid, out$fence_upper, lty = 2, col = 'red', lwd = 2 )
    # lines( time_grid, out$fence_lower, lty = 2, col = 'red', lwd = 2 )

    # Plotting closest data to fences
    id_uppermost = which.min( apply( - t( t( fData$values[ - ID_out, ]) +
                                            out$fence_upper ), 1, min )  )
    id_lowermost = which.min( apply(   t( t( fData$values[ - ID_out, ]) -
                                            out$fence_lower ), 1, min )  )

    lines( time_grid, fData$values[ -ID_out, ][ id_uppermost, ], lty = 1,
           col = col_fence_structure, lwd = 3 )
    lines( time_grid, fData$values[ -ID_out, ][ id_lowermost, ], lty = 1,
           col = col_fence_structure, lwd = 3 )

    # Plotting vertical whiskers
    half.time_grid = which.min( abs( time_grid - 0.5 ) )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$max_envelope_central[ half.time_grid ],
              fData$values[ -ID_out, half.time_grid ][ id_uppermost ] ),
           lty = 1, col = col_fence_structure, lwd = 3 )

    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$min_envelope_central[ half.time_grid ],
              fData$values[ -ID_out, half.time_grid ][ id_lowermost ] ),
           lty = 1, col = col_fence_structure, lwd = 3 )
  }

  return( list( Depth = Depths,
                ID_outliers = ID_out ) )
}


.fbplot = function( time_grid, Data, Depths = 'MBD', Fvalue = 1.5 )
{
  # Number of observations
  DataN = nrow( Data )

  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths = eval( parse( text = paste( Depths, '( Data )', sep = '' ) ) )
  } else {
    stopifnot( length( Depths ) == DataN )
  }

  Data_center = Data[ which.max( Depths ), ]

  id_central_region = which( Depths >= quantile( Depths, prob = 0.5 ) )

  max_envelope_central = apply( Data[ id_central_region, ], 2, max )
  min_envelope_central = apply( Data[ id_central_region, ], 2, min )

  fence_upper = ( max_envelope_central - Data_center ) * Fvalue + Data_center
  fence_lower = ( min_envelope_central - Data_center ) * Fvalue + Data_center

  ID_outlying = which( apply( Data, 1, function(x)( any( x > fence_upper  ) |
                                                      any ( x < fence_lower ) ) ) )

  return( list( ID_out = ID_outlying,
                min_envelope_central = min_envelope_central,
                max_envelope_central = max_envelope_central,
                fence_lower = fence_lower,
                fence_upper = fence_upper ) )
}
