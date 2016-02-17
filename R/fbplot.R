#' Functional boxplot
#'
#' \code{fbplot} displays the functional boxplot of a dataset of functional data.
#'
#' @param Data the univariate or multivariate functional dataset whose functional
#' boxplot is desired
#' @param Depths either a vector containing the depths for each element of the
#' dataset, or a string containing the name of the method you want to use to
#' compute it. In this case the name of the method must be included in the
#' caller's environment
#' @param Fvalue the value of the inflation factor F, default is F = 1.5
#' @param adjust either FALSE if you would like the default value for the
#' inflation factor, F = 1.5, to be used, or a list specifying the parameters
#' required by the adjustment.
#' @param display either a logical value indicating wether you want the
#' functional boxplot to be displayed, or the number of the graphical devices
#' where you want the boxplot to be plotted.
#' @param xlab the label to use on the x axis when displaying the functional
#' boxplot
#' @param ylab the label(s) to use on the y axis when displaying the functional
#' boxplot
#' @param main the main title(s) to used when displaying the functional boxplot
#' @param ... additional graphical parameters to be used in plotting functions
#' @param verbose logical flag to enable verbosity in the adjustment process
#'
fbplot = function( Data,
                   Depths = 'MBD',
                   Fvalue = 1.5,
                   adjust = FALSE,
                   display = TRUE,
                   xlab = NULL,
                   ylab = NULL,
                   main = NULL,
                   ...,
                   verbose = FALSE )
{
  UseMethod( 'fbplot', Data )
}

#' Functional boxplot for univariate functional data
#'
#' \code{fbplot.fData} displays the functional boxplot of a dataset of
#' functional data.
#'
#' @param Data the univariate functional dataset whose functional boxplot is
#' desired
#' @param Depths either a vector containing the depths for each element of the
#' dataset, or a string containing the name of the method you want to use to
#' compute it. In this case the name of the method must be included in the
#' caller's environment
#' @param Fvalue the value of the inflation factor F, default is F = 1.5,
#' @param adjust either FALSE if you would like the default value for the
#' inflation factor, F = 1.5, to be used, or a list specifying the parameters
#' required by the adjustment.
#' @param display either a logical value indicating wether you want the
#' functional boxplot to be displayed, or the number of the graphical devices
#' where you want the boxplot to be plotted.
#' @param ... additional graphical parameters to be used in plotting functions
#' @param verbose logical flag to enable verbosity in the adjustment process
#'
fbplot.fData = function( Data,
                         Depths = 'MBD',
                         Fvalue = 1.5,
                         adjust = FALSE,
                         display = TRUE,
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         ...,
                         verbose = FALSE )
{
  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec == 'MBD' )
    {
      Depths = MBD( Data$values, manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( Data$values )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == Data$N )
  }

  if( ! is.list( adjust ) )
  {
    # Plain functional boxplot with default F value: F = 1.5

    out = .fbplot_fData( time_grid, Data$values, Depths, Fvalue )

  } else {

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         5 * Data$N,
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
    Cov = robustbase::covOGK( Data$values, sigmamu = robustbase::s_Qn )$cov

    # Cholesky factor
    CholCov <- chol( Cov )

    # Centerline of the dataset
    centerline = Data$values[ which.max( Depths ), ]

    Fvalues = rep( 0, N_trials )

    cost_functional = function( F_curr )( length(
      .fbplot_fData( time_grid,
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

    out = .fbplot_fData( time_grid, Data$values, Depths, Fvalue = Fvalue  )
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
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( Data$N - length( ID_out ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( ID_out ) )
    col_envelope = set_alpha( 'blue', alpha = 0.4 )
    col_center = set_alpha( 'blue', alpha = 1 )
    col_fence_structure = 'darkblue'

    time_grid = seq( Data$t0, Data$tP, length.out = Data$P )

    xlab = ifelse( is.null( xlab ), '', xlab )
    ylab = ifelse( is.null( ylab ), '', ylab )
    main = ifelse( is.null( main ), '', main )

    if( length( ID_out ) > 0 )
    {
      # Plotting non-outlying data
      matplot( time_grid,
               t( Data$values[ - ID_out, ] ), lty = 1, type = 'l',
               col = col_non_outlying,
               ylim = range( Data$values ),
               xlab = xlab, ylab = ylab, main = main, ... )

      # Computing maximum and minimum envelope
      max_envelope_limit = apply( Data$values[ - ID_out, ], 2, max )
      min_envelope_limit = apply( Data$values[ - ID_out, ], 2, min )
    } else {
      # Plotting all data
      matplot( time_grid,
               t( Data$values ), lty = 1, type = 'l',
               col = col_non_outlying,
               ylim = range( Data$values ),
               xlab = xlab, ylab = ylab, main = main, ... )

      # Computing maximum and minimum envelope
      max_envelope_limit = apply( Data$values, 2, max )
      min_envelope_limit = apply( Data$values, 2, min )
    }


    # Filling in the central envelope

    polygon( c(time_grid, rev( time_grid) ),
             c( out$min_envelope_central, rev( out$max_envelope_central ) ),
             col = col_envelope, border = NA)
    lines( time_grid, out$max_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    lines( time_grid, out$min_envelope_central, lty = 1, col = col_envelope, lwd = 3 )

    # Plotting the sample median
    lines( time_grid, Data$values[ which.max( Depths ), ], lty = 1, type = 'l',
           col = col_center, lwd = 3)

    lines( time_grid, max_envelope_limit, lty = 1,
    col = col_fence_structure, lwd = 3 )
    lines( time_grid, min_envelope_limit, lty = 1,
    col = col_fence_structure, lwd = 3 )

    # Plotting vertical whiskers
    half.time_grid = which.min( abs( time_grid - 0.5 ) )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$max_envelope_central[ half.time_grid ],
              max_envelope_limit[ half.time_grid ] ),
           lty = 1, col = col_fence_structure, lwd = 3 )

    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$min_envelope_central[ half.time_grid ],
              min_envelope_limit[ half.time_grid ] ),
           lty = 1, col = col_fence_structure, lwd = 3 )

    # Plotting outlying data
    if( length( ID_out ) > 0 )
    {
      matplot( time_grid, t( Data$values[ ID_out, ] ), lty = 1, type = 'l',
               col = col_outlying, lwd = 3, add = T )
    }
  }

  return( list( Depth = Depths,
                ID_outliers = ID_out ) )
}


.fbplot_fData = function( time_grid, Data, Depths = 'MBD', Fvalue = 1.5 )
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


fbplot.mfData = function( Data,
                          Depths = list( def = 'MBD',
                                         weights = 'uniform' ),
                          Fvalue = 1.5,
                          adjust = FALSE,
                          display = TRUE,
                          xlab = NULL,
                          ylab = NULL,
                          main = NULL,
                          ...,
                          verbose = FALSE )
{
  if( is.list( adjust ) )
  {
    stop( ' Error in fbplot.mfData: for now the adjustment support is not
provided in the multivariate version of the functional boxplot' )
  }

  listOfValues = toListOfValues( Data )

  # Checking if depths have already been provided or must be computed
  if( is.list( Depths ) )
  {
    if( length( Depths ) != 2 & ! identical( names( Depths ),
                                             c( 'def', 'weights' ) ) )
    {
      stop( " Error in fbplot.mfData: you have to provide both a specification
            for the depth definition and one for the set of weights.")
    }

    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec[[ 'def' ]] == 'MBD' )
    {
      Depths = multiMBD( listOfValues,
                         weights = Depths_spec[[ 'weights' ]],
                         manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( listOfValues )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == Data$N )
  }

  out = .fbplot_mfData( time_grid, listOfValues, Depths, Fvalue  )

  ID_out = out$ID_out

  if( is.numeric( display ) )
  {
    dev.set( display )
  }

  if( ! display == FALSE )
  {

    # Subdividing the graphical window
    mfrow_rows = ceiling( Data$L / 2 )
    mfrow_cols = 2

    par( mfrow = c( mfrow_rows, mfrow_cols ) )

    # Creating color palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( Data$N - length( ID_out ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( ID_out ) )
    col_envelope = set_alpha( 'blue', alpha = 0.4 )
    col_center = set_alpha( 'blue', alpha = 1 )
    col_fence_structure = 'darkblue'

    time_grid = seq( Data$t0, Data$tP, length.out = Data$P )

    xlab = ifelse( is.null( xlab ), '', xlab )
    ylab = ifelse( is.null( ylab ), '', ylab )
    main = ifelse( is.null( main ), '', main )

    for( iL in 1 : Data$L )
    {
      Data_curr = Data$fDList[[ iL ]]$values

      if( length( ID_out ) > 0 )
      {
        # Plotting non-outlying data
        matplot( time_grid,
                 t( Data_curr[ - ID_out, ] ), lty = 1, type = 'l',
                 col = col_non_outlying,
                 ylim = range( Data_curr ),
                 xlab = xlab, ylab = ylab, main = main, ... )

        # Computing maximum and minimum envelope
        max_envelope_limit = apply( Data_curr[ - ID_out, ], 2, max )
        min_envelope_limit = apply( Data_curr[ - ID_out, ], 2, min )
      } else {
        # Plotting all data
        matplot( time_grid,
                 t( Data_curr ), lty = 1, type = 'l',
                 col = col_non_outlying,
                 ylim = range( Data_curr ),
                 xlab = xlab, ylab = ylab, main = main, ... )

        # Computing maximum and minimum envelope
        max_envelope_limit = apply( Data_curr, 2, max )
        min_envelope_limit = apply( Data_curr, 2, min )
      }

      # Filling in the central envelope

      polygon( c(time_grid, rev( time_grid) ),
               c( as.numeric( out$min_envelope_central[ iL, ] ),
                  rev( as.numeric( out$max_envelope_central[ iL, ] ) ) ),
               col = col_envelope, border = NA)
      lines( time_grid, as.numeric( out$max_envelope_central[ iL, ] ),
             lty = 1, col = col_envelope, lwd = 3 )
      lines( time_grid, as.numeric( out$min_envelope_central[ iL, ] ),
             lty = 1, col = col_envelope, lwd = 3 )

      # Plotting the sample median
      lines( time_grid, Data_curr[ which.max( Depths ), ], lty = 1, type = 'l',
             col = col_center, lwd = 3)

      lines( time_grid, max_envelope_limit, lty = 1,
             col = col_fence_structure, lwd = 3 )
      lines( time_grid, min_envelope_limit, lty = 1,
             col = col_fence_structure, lwd = 3 )

      # Plotting vertical whiskers
      half.time_grid = which.min( abs( time_grid - 0.5 ) )
      lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
             c( out$max_envelope_central[ iL, half.time_grid ],
                max_envelope_limit[ half.time_grid ] ),
             lty = 1, col = col_fence_structure, lwd = 3 )

      lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
             c( out$min_envelope_central[ iL, half.time_grid ],
                min_envelope_limit[ half.time_grid ] ),
             lty = 1, col = col_fence_structure, lwd = 3 )


      # Plotting outlying data
      if( length( ID_out ) > 0 )
      {
        matplot( time_grid, t( Data_curr[ ID_out, ] ), lty = 1, type = 'l',
                 col = col_outlying, lwd = 3, add = T )
      }
    }
  }

  return( list( Depth = Depths,
                ID_outliers = ID_out ) )
}


.fbplot_mfData = function( time_grid,
                           listOfValues,
                           Depths = list( def = 'MBD',
                                          weights = 'uniform' ),
                           Fvalue = 1.5 )
{

  L = length( listOfValues )
  N = nrow( listOfValues[[ 1 ]] )
  P = ncol( listOfValues[[ 2 ]] )

  # Checking if depths have already been provided or must be computed
  if( is.list( Depths ) )
  {
    if( length( Depths ) != 2 & ! identical( names( Depths ),
                                             c( 'def', 'weights' ) ) )
    {
      stop( " Error in .fbplot_mfData: you have to provide both a specification
            for the depth definition and one for the set of weights.")
    }

    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec[[ 'def' ]] == 'MBD' )
    {
      Depths = multiMBD( listOfValues,
                         weights = Depths_spec[[ 'weights' ]],
                         manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( listOfValues )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == N )
  }

  # Nice hack to extract selected row from each matrix of the list and
  # concatenate the result into a row-major matrix
  Data_center = t( sapply( listOfValues, `[`,  which.max( Depths ), 1 : P ) )

  id_central_region = which( Depths >= quantile( Depths, prob = 0.5 ) )

  max_envelope_central = t( sapply( 1 : L, function( i ) (
    apply( listOfValues[[ i ]][ id_central_region, ], 2, max ) ) ) )

  min_envelope_central = t( sapply( 1 : L, function( i ) (
    apply( listOfValues[[ i ]][ id_central_region, ], 2, min ) ) ) )

  fence_upper = ( max_envelope_central - Data_center ) * Fvalue + Data_center
  fence_lower = ( min_envelope_central - Data_center ) * Fvalue + Data_center

  ID_outlying = unique( unlist( sapply( 1 : L, function( iL ) ( which(
    apply( listOfValues[[ iL ]], 1,
           function( x ) ( any( x > as.numeric( fence_upper[ iL, ] ) ) |
                             any( x < as.numeric( fence_lower[ iL, ] ) ) ) ) )
  ) ) ) )

  return( list( ID_out = ID_outlying,
                min_envelope_central = min_envelope_central,
                max_envelope_central = max_envelope_central,
                fence_lower = fence_lower,
                fence_upper = fence_upper ) )
}
