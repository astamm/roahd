#' Outliergram for multivariate functional datasets
#'
#' This function performs the outliergram of a multivariate functional dataset.
#'
#' The method applies the extension of the univariate outliergram to the case of multivariate
#' functional datasets. Differently from the function for the univariate case, only the
#' outliergram plot is displayed.
#'
#' @section Adjustment:
#'
#' Differently from the case of univariate functional data, in this case the function does not apply
#' an automatic tuning of the F parameter, since the related procedure would become computationally
#' too heavy for general datasets. If a good value of F is sought, it is recommended to run several
#' trials of the outliergram and manually select the best value.
#'
#' @param mfData the multivariate functional dataset whose outliergram has to be
#' determined;
#' @param MBD_data a vector containing the MBD for each element of the dataset; If missing, MBDs are computed with the specified choice of weights;
#' @param weights: the weights choice to be used to compute multivariate MBDs and MEIs;
#' @param MEI_data a vector containing the MEI for each element of the dataset.
#' If not not provided, MEIs are computed;
#' @param Fvalue the \eqn{F} value to be used in the procedure that finds the
#' shape outliers by looking at the lower parabolic limit in the outliergram.
#' Default is \code{1.5};
#' @param display either a logical value indicating wether you want the
#' outliergram to be displayed, or the number of the graphical device
#' where you want the outliergram to be displayed;
#' @param xlab the label to use on the x axis in the outliergram plot;
#' @param ylab the label to use on the x axis in the outliergram plot;
#' @param main the title to use in the outliergram;
#'
#' @references
#'
#' Ieva, F. & Paganoni, A.M. Stat Papers (2017). https://doi.org/10.1007/s00362-017-0953-1.
#'
#' @seealso \code{\link{outliergram}}, \code{\link{mfData}}, \code{\link{MBD}},
#' \code{\link{MEI}}
#'
#' @examples
#'
#' N = 2e2
#' P = 1e2
#'
#' t0 = 0
#' t1 = 1
#'
#' set.seed(1)
#'
#' # Defining the measurement grid
#' grid = seq( t0, t1, length.out = P )
#'
#' # Generating an exponential covariance matrix to be used in the simulation of
#' # the functional datasets (see the related help for details)
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.2)
#'
#' # Simulating the measurements of two univariate functional datasets with
#' # required center and covariance function
#' f1 = function(x) x * ( 1 - x )
#' f2 = function(x) x^3
#' Data = generate_gauss_mfdata( N, L = 2,
#'                               centerline = matrix(c(sin(2 * pi * grid),
#'                                                     cos(2 * pi * grid)), nrow=2, byrow=TRUE),
#'                               listCov = list(C, C), correlations = 0.1 )
#'
#' # Building the mfData object
#' mfD = mfData( grid, Data )
#'
#'
#' dev.new()
#' out = multivariate_outliergram(mfD, Fvalue = 13, shift=FALSE)
#' col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
#'                                     l = 60 )( N - length( out$ID_outliers ) )
#' col_non_outlying = set_alpha( col_non_outlying, 0.5 )
#' col_outlying = scales::hue_pal( h = c( - 90, 180  ),
#'                                 c = 150 )( length( out$ID_outliers ) )
#' colors = rep('black', N)
#' colors[out$ID_outliers] = col_outlying
#' colors[colors == 'black'] = col_non_outlying
#'
#' lwd = rep(1, N)
#' lwd[out$ID_outliers] = 2
#'
#' dev.new()
#' plot(mfD, col=colors, lwd=lwd)
#'
#' @export
multivariate_outliergram = function( mfData,
                          MBD_data = NULL,
                          MEI_data = NULL,
                          weights='uniform',
                          p_check = 0.05, q_low = 0, q_high = 1,
                          Fvalue = 1.5, shift = TRUE,
                          display = TRUE, xlab = NULL, ylab = NULL, main = NULL )
{
  N = mfData$N

  out = .outliergram_mfData( mfData,
                             MBD_data = MBD_data, MEI_data = MEI_data,
                             p_check = p_check, q_low = q_low, q_high = q_high,
                             Fvalue = Fvalue,
                             shift = shift )

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Plot
  if( display )
  {

    if( is.null( xlab ) )
    {
      xlab = 'MEI'
    }

    if( is.null( ylab ) )
    {
      ylab = 'MBD'
    }

    if( is.null( main ) )
    {
      main = 'Outliergram'
    }

    # Setting up palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( length( out$ID_NO ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( out$ID_SO ) )

    dev.cur()
    # Plotting outliergram
    ## Upper parabolic limit
    grid_1D = seq( 0, 1, length.out = 100 )

    plot( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2,
          lty = 2, type = 'l', col = 'darkblue', lwd = 2,
          ylim = c( 0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 ),
          xlab = xlab,
          ylab = ylab,
          main = main )

    if( length( out$ID_SO ) > 0 )
    {
      points( out$MEI_data[ - out$ID_SO ], out$MBD_data[ - out$ID_SO ],
              pch = 16, col = col_non_outlying )
      points( out$MEI_data[ out$ID_SO ], out$MBD_data[ out$ID_SO ],
              pch = 16, cex = 1.5, col = col_outlying )
      for( idOut in out$ID_SO )
      {
        text( out$MEI_data[ idOut ],
              out$MBD_data[ idOut ] + 0.5 / 30,
              idOut,
              col = col_outlying[ match( idOut, out$ID_SO ) ] )
      }
    } else {
      points( out$MEI_data, out$MBD_data,
              pch = 16, col = col_non_outlying )
    }

    # lower parabolic limit
    if( Fvalue == 1.5 )
    {
      lines( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
               out$Q_d3 - 1.5 * out$IQR_d,
             lty = 2, lwd = 2, col = 'lightblue' )
    }
    else
    {
      lines( grid_1D,  a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
               Fvalue * out$Q_d1,
             lty = 2, lwd = 2, col = 'lightblue' )
    }
  }
  return( list( Fvalue = Fvalue,
                d = out$d,
                ID_outliers = out$ID_SO ) )
}

.outliergram_mfData = function( mfData, MBD_data = NULL, MEI_data = NULL, weights='uniform',
                                p_check = 0.05, q_low = 0, q_high = 1,
                                Fvalue = 1.5, shift = TRUE )
{
  N = mfData$N

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Computing MBD
  if( is.null( MBD_data ) ){
    MBD_data = multiMBD( toListOfValues(mfData), weights = weights  )
  }

  # Computing MEI
  if( is.null( MEI_data ) ){
    MEI_data = multiMEI( toListOfValues(mfData), weights = weights )
  }

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves

  if( Fvalue == 1.5 )
  {
    ID_shape_outlier = which( d >= Q_d3 + 1.5 * IQR_d )
  } else {
    ID_shape_outlier = which( d >= Fvalue * Q_d1 )
  }

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : N, ID_shape_outlier )

  if( shift )
  {
    stop('Unsupported feature, yet.')
  }

  return( list( ID_SO = ID_shape_outlier,
                ID_NO = ID_non_outlying,
                MEI_data = MEI_data,
                MBD_data = MBD_data,
                Q_d3 = Q_d3,
                Q_d1 = Q_d1,
                IQR_d = IQR_d,
                d = d) )

}

