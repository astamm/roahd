#' Outliergram for univariate functional datasets
#'
#' This function performs the outliergram of a univariate functional dataset,
#' possibly with an adjustment of the true positive rate of outliers discovered
#' under assumption of gaussianity.
#'
#' @section Adjustment:
#'
#' When the adjustment option is selected, the value of \eqn{F} is optimized for
#' the univariate functional dataset provided with \code{fData}. In practice,
#' a number \code{adjust$N_trials} of times a synthetic population
#' (of size \code{adjust$trial_size} with the same covariance (robustly
#' estimated from data) and centerline as \code{fData} is simulated without
#' outliers and each time an optimized value \eqn{F_i} is computed so that a
#' given proportion (\code{adjust$TPR}) of observations is flagged as outliers.
#' The final value of \code{F} for the outliergram is determined as an average
#' of \eqn{F_1, F_2, \ldots, F_{N_{trials}}}. At each time step the optimization
#' problem is solved using \code{stats::uniroot} (Brent's method).
#'
#' @param fData the univariate functional dataset whose outliergram has to be
#' determined.
#' @param MBD_data a vector containing the MBD for each element of the dataset.
#' If missing, MBDs are computed.
#' @param MEI_data a vector containing the MEI for each element of the dataset.
#' If not not provided, MEIs are computed.
#' @param p_check percentage of observations with either low or high MEI to be
#' checked for outliers in the secondary step (shift towards the center of the
#' dataset).
#' @param Fvalue the \eqn{F} value to be used in the procedure that finds the
#' shape outliers by looking at the lower parabolic limit in the outliergram.
#' Default is \code{1.5}. You can also leave the default value and, by providing
#' the parameter \code{adjust}, specify that you want \code{Fvalue} to be
#' adjusted for the dataset provided in \code{fData}.
#' @param adjust either \code{FALSE} if you would like the default value for the
#' inflation factor, \eqn{F = 1.5}, to be used, or a list specifying the
#' parameters required by the adjustment.
#'  \itemize{
#'  \item{"\code{N_trials}"}{: the number of repetitions of the adjustment
#'  procedure based on the simulation of a gaussian population of functional
#'  data, each one producing an adjusted value of \eqn{F}, which will lead
#'  to the averaged adjusted value \eqn{\bar{F}}. Default is 20;}
#'  \item{"\code{trial_size}"}{: the number of elements in the gaussian
#'  population of functional data that will be simulated at each repetition of
#'  the adjustment procedure. Default is \code{5 * fData$N};}
#'  \item{"\code{TPR}"}{: the True Positive Rate of outliers, i.e. the proportion
#'  of observations in a dataset without shape outliers that have to be considered
#'  outliers. Default is \code{2 * pnorm( 4 * qnorm( 0.25 ) )};}
#'  \item{"\code{F_min}"}{: the minimum value of \eqn{F}, defining the left
#'  boundary for the optimization problem aimed at finding, for a given dataset
#'  of simulated gaussian data associated to \code{fData}, the optimal value of
#'  \eqn{F}. Default is 0.5;}
#'  \item{"\code{F_max}"}{: the maximum value of \eqn{F}, defining the right
#'  boundary for the optimization problem aimed at finding, for a given dataset
#'  of simulated gaussian data associated to \code{fData}, the optimal value of
#'  \eqn{F}. Default is 20;}
#'  \item{"\code{tol}"}{: the tolerance to be used in the optimization problem
#'  aimed at finding, for a given dataset of simulated gaussian data associated
#'  to \code{fData}, the optimal value of \eqn{F}. Default is \code{1e-3};}
#'  \item{"\code{maxiter}"}{: the maximum number of iterations to solve the
#'  optimization problem aimed at finding, for a given dataset of simulated
#'  gaussian data associated to \code{fData}, the optimal value of \eqn{F}.
#'  Default is \code{100};}
#'  \item{"\code{VERBOSE}"}{: a parameter controlling the verbosity of the
#'  adjustment process;}
#'  }
#' @param display either a logical value indicating whether you want the
#' outliergram to be displayed, or the number of the graphical device
#' where you want the outliergram to be displayed.
#' @param xlab a list of two labels to use on the x axis when displaying the
#' functional dataset and the outliergram
#' @param ylab a list of two labels to use on the y axis when displaying the
#' functional dataset and the outliergram;
#' @param main a list of two titles to be used on the plot of the functional
#' dataset and the outliergram;
#' @param ... additional graphical parameters to be used \emph{only} in the plot
#' of the functional dataset
#'
#' @return
#'
#' Even when used graphically to plot the outliergram, the function returns a
#' list containing:
#' \itemize{
#' \item{\code{Fvalue}}{: the value of the parameter F used;}
#' \item{\code{d}}{: the vector of values of the parameter \eqn{d} for each observation
#' (distance to the parabolic border of the outliergram);}
#' \item{\code{ID_outliers}}{: the vector of observations id corresponding to outliers.}}
#'
#' @references
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and visualization
#' for functional data: the outliergram, \emph{Biostatistics}, 15(4), 603-619.
#'
#' @seealso \code{\link{fData}}, \code{\link{MEI}}, \code{\link{MBD}},
#' \code{\link{fbplot}}
#'
#' @examples
#'
#'
#' set.seed( 1618 )
#'
#' N = 200
#' P = 200
#' N_extra = 4
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.8 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 4 * pi * grid ),
#'                              Cov = Cov )
#'
#' Data_extra = array( 0, dim = c( N_extra, P ) )
#'
#' Data_extra[ 1, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid + pi / 2 ),
#'                                          Cov = Cov )
#'
#' Data_extra[ 2, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid - pi / 2 ),
#'                                          Cov = Cov )
#'
#' Data_extra[ 3, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid + pi/ 3 ),
#'                                          Cov = Cov )
#'
#' Data_extra[ 4, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid - pi / 3),
#'                                          Cov = Cov )
#' Data = rbind( Data, Data_extra )
#'
#' fD = fData( grid, Data )
#'
#' outliergram( fD, display = TRUE )
#'
#' outliergram( fD, Fvalue = 2.5, display = TRUE )
#' \dontrun{
#' outliergram( fD,
#'              adjust = list( N_trials = 10,
#'                             trial_size = 5 * nrow( Data ),
#'                             TPR = 0.01,
#'                             VERBOSE = FALSE ),
#'              display = TRUE )
#' }
#'
#' @importFrom grDevices dev.set dev.cur
#' @importFrom stats cor pnorm rnorm qnorm uniroot
#' @importFrom graphics text lines polygon plot points matplot par
#' @importFrom dplyr filter group_by summarize
#' @importFrom magrittr %>%
#'
#' @export
outliergram = function( fData, MBD_data = NULL, MEI_data = NULL, p_check = 0.05,
                        Fvalue = 1.5,
                        adjust = FALSE, display = TRUE,
                        xlab = NULL, ylab = NULL, main = NULL, ... )
{
  N = fData$N

  grid = seq( fData$t0,
              fData$tP,
              length.out = fData$P )

  if( ! is.list( adjust ) )
  {
    # Plain outliergram with default F value: F = 1.5
    stopifnot( is.numeric(Fvalue) )

    out = .outliergram( fData,
                        MBD_data = MBD_data, MEI_data = MEI_data,
                        p_check = p_check,
                        Fvalue = Fvalue,
                        shift = TRUE )
  } else {

    nodenames = c( 'N_trials', 'trial_size', 'TPR', 'F_min', 'F_max',
                   'tol', 'maxiter', 'VERBOSE' )
    unused = setdiff( names( adjust ), nodenames )

    # Checking for unused parameters
    if( length( unused ) > 0 )
    {
      for( i in unused )
        warning( 'Warning: unused parameter ', i, ' in adjust argument of outliergram' )
    }

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         5 * fData$N,
                         adjust$trial_size )

    TPR = ifelse( is.null( adjust$TPR ),
                  2 * pnorm( 4 * qnorm( 0.25 ) ),
                  adjust$TPR )

    F_min = ifelse( is.null( adjust$F_min ),
                    0.5,
                    adjust$F_min )

    F_max= ifelse( is.null( adjust$F_max ),
                    20,
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

    Cov = robustbase::covOGK( fData$values, sigmamu = robustbase::s_Qn )$cov

    CholCov <- chol( Cov )

    if( is.null ( MBD_data ) )
    {
      MBD_data = MBD( fData$values, manage_ties = TRUE )
    }

    centerline = fData$values[ which.max( MBD_data ), ]

    Fvalues = rep( 0, N_trials )

    obj_function = function( F_curr )( length(
      .outliergram( fData_gauss,
                    MBD_data = NULL,
                    MEI_data = NULL,
                    Fvalue = F_curr,
                    shift = FALSE )$ID_SO ) / trial_size - TPR )

    for( iTrial in 1 : N_trials )
    {
      if( VERBOSE > 0 )
      {
        cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
      }

      fData_gauss = fData( grid,
                           generate_gauss_fdata( N = trial_size,
                                                 centerline = centerline,
                                                 CholCov = CholCov ) )

      if( VERBOSE > 0 )
      {
        cat( ' * * * * beginning optimization\n' )
      }

      opt = uniroot( obj_function,
                     interval = c( F_min, F_max ),
                     tol = tol,
                     maxiter = maxiter )

      Fvalues[ iTrial ] = opt$root
    }

    Fvalue = mean( Fvalues )

    out = .outliergram( fData, MBD_data, MEI_data, p_check, Fvalue = Fvalue, shift = TRUE  )
  }

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Plot
  if( display )
  {

    if( is.null( xlab ) )
    {
      xlab = list( '', 'MEI' )
    }

    if( is.null( ylab ) )
    {
      ylab = list( '', 'MBD' )
    }

    if( is.null( main ) )
    {
      main = list( '', 'Outliergram' )
    }

    # Setting up palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( length( out$ID_NO ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    if( length( out$ID_SO ) > 0 )
    {
      col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                      c = 150 )( length( out$ID_SO ) )

    } else
    {
      col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                      c = 150 )( 1 )

    }

    dev.cur()
    par( mfrow = c( 1, 2 ) )

    # Plotting functional data
    if( length( out$ID_SO ) > 0 )
    {
      matplot( grid, t( fData$values[ - out$ID_SO, ] ), type = 'l', lty = 1,
               ylim = range( fData$values ),
               col = col_non_outlying,
               xlab = xlab[[1]],
               ylab = ylab[[1]],
               main = main[[1]],
               ... )
      matplot( grid, t( toRowMatrixForm( fData$values[ out$ID_SO, ] ) ),
               type = 'l', lty = 1, lwd = 3, ylim = range( fData$values ),
               col = col_outlying, add = TRUE )
    } else {
      matplot( grid, t( fData$values ), type = 'l', lty = 1,
               ylim = range( fData$values ),
               col = col_non_outlying,
               xlab = xlab[[1]],
               ylab = ylab[[1]],
               main = main[[1]],
               ... )
    }


    # Adding text labels with curve ID
    w_spacing = diff( range( grid ) ) / ( 2 * length( out$ID_SO ) )

    for( iOut in seq_along( out$ID_SO ) )
    {
      text( grid[ 1 ] + ( 2 * iOut - 1 ) * w_spacing,
            fData$values[ out$ID_SO[ iOut ],
                  which.min( abs( grid - grid[ 1 ] -
                                    ( 2 * iOut - 1 ) * w_spacing ) ) ] +
              diff( range( fData$values[ out$ID_SO[ iOut ]  ] ) ) / 30,
            out$ID_SO[ iOut ],
            col = col_outlying[ iOut ] )
    }

    # Plotting outliergram

    # Upper parabolic limit
    grid_1D = seq( 0, 1, length.out = 100 )

    plot( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2,
          lty = 2, type = 'l', col = 'darkblue', lwd = 2,
          ylim = c( 0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 ),
          xlab = xlab[[2]],
          ylab = ylab[[2]],
          main = main[[2]] )

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
    lines( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
             out$Q_d3 - Fvalue * out$IQR_d,
           lty = 2, lwd = 2, col = 'lightblue' )
  }

  return( list( Fvalue = Fvalue,
                d = out$d,
                ID_outliers = out$ID_SO ) )
}

.outliergram = function( fData, MBD_data = NULL, MEI_data = NULL,
                         p_check = 0.05,
                         Fvalue = 1.5, shift = TRUE )
{
  N = fData$N

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Computing MBD
  if( is.null( MBD_data ) ){

    MBD_data = MBD( fData$values )
  }

  # Computing MEI
  if( is.null( MEI_data ) )
  {
    MEI_data = MEI( fData$values )
  }

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves
  ID_shape_outlier = which( d >= Q_d3 + Fvalue * IQR_d )

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : nrow( fData$values ), ID_shape_outlier )

  if( shift )
  {
    # Low MEI curves will be checked for upward shift
      ID_non_outlying_Low_MEI = ID_non_outlying[
        which( MEI_data[ - ID_shape_outlier ] <=
                 quantile( MEI_data,
                           probs = p_check ) ) ]

      # High MEI curves will be checked for downward shift
      ID_non_outlying_High_MEI = ID_non_outlying[
        which( MEI_data[ - ID_shape_outlier ] >=
                 quantile( MEI_data, probs = 1 - p_check ) ) ]


      # Manage high MEI data
      lst = manage_high_MEI_data(fData = fData,
                                 ID_non_outlying_High_MEI = ID_non_outlying_High_MEI,
                                 ID_shape_outlier = ID_shape_outlier,
                                 ID_non_outlying = ID_non_outlying,
                                 Q_d1 = Q_d1,
                                 Q_d3 = Q_d3,
                                 IQR_d = IQR_d,
                                 Fvalue = Fvalue)

      ID_non_outlying = lst[['ID_non_outlying']]
      ID_shape_outlier = lst[['ID_shape_outlier']]

      # Manage low MEI data
      lst = manage_low_MEI_data(fData = fData,
                                ID_non_outlying_Low_MEI = ID_non_outlying_Low_MEI,
                                ID_shape_outlier = ID_shape_outlier,
                                ID_non_outlying = ID_non_outlying,
                                Q_d1 = Q_d1,
                                Q_d3 = Q_d3,
                                IQR_d = IQR_d,
                                Fvalue = Fvalue)

      ID_non_outlying = lst[['ID_non_outlying']]
      ID_shape_outlier = lst[['ID_shape_outlier']]


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
#' @param MBD_data a vector containing the MBD for each element of the dataset; If missing, MBDs
#' are computed with the specified choice of weights;
#' @param weights the weights choice to be used to compute multivariate MBDs and MEIs;
#' @param MEI_data a vector containing the MEI for each element of the dataset.
#' If not not provided, MEIs are computed;
#' @param p_check percentage of observations with either low or high MEI to be
#' checked for outliers in the secondary step (shift towards the center of the
#' dataset).
#' @param Fvalue the \eqn{F} value to be used in the procedure that finds the
#' shape outliers by looking at the lower parabolic limit in the outliergram.
#' Default is \code{1.5};
#' @param shift whether to apply the shifting algorithm to properly manage observations having low
#' or high MEI. Default is TRUE.
#' @param display either a logical value indicating whether you want the
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
#' out = multivariate_outliergram(mfD, Fvalue = 2., shift=TRUE)
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
#' @importFrom dplyr filter group_by summarize
#' @importFrom magrittr %>%
#'
#' @export
multivariate_outliergram = function( mfData,
                                     MBD_data = NULL,
                                     MEI_data = NULL,
                                     weights='uniform',
                                     p_check = 0.05,
                                     Fvalue = 1.5, shift = TRUE,
                                     display = TRUE, xlab = NULL, ylab = NULL, main = NULL )
{
  N = mfData$N

  out = .outliergram_mfData( mfData,
                             MBD_data = MBD_data, MEI_data = MEI_data,
                             p_check = p_check,
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

    if( length( out$ID_SO ) > 0 )
    {
      col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                      c = 150 )( length( out$ID_SO ) )

    } else {
      col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                      c = 150 )( 1 )
    }

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
      lines( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
               out$Q_d3 - Fvalue * out$IQR_d,
             lty = 2, lwd = 2, col = 'lightblue' )
  }
  return( list( Fvalue = Fvalue,
                d = out$d,
                ID_outliers = out$ID_SO ) )
}

.outliergram_mfData = function( mfData, MBD_data = NULL, MEI_data = NULL, weights='uniform',
                                p_check = 0.05,
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
  ID_shape_outlier = which( d >= Q_d3 + Fvalue * IQR_d )

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : N, ID_shape_outlier )

  if( shift )
  {
    for( iDim in seq_along(mfData$fDList))
    {
      # Low MEI curves will be checked for upward shift
      ID_non_outlying_Low_MEI = ID_non_outlying[
        which( MEI_data[ - ID_shape_outlier ] <=
                 quantile( MEI_data,
                           probs = p_check ) ) ]

      # High MEI curves will be checked for downward shift
      ID_non_outlying_High_MEI = ID_non_outlying[
        which( MEI_data[ - ID_shape_outlier ] >=
                 quantile( MEI_data, probs = 1 - p_check ) ) ]


      # Manage high MEI data
      lst = manage_high_MEI_data(fData = mfData$fDList[[iDim]],
                                 ID_non_outlying_High_MEI = ID_non_outlying_High_MEI,
                                 ID_shape_outlier = ID_shape_outlier,
                                 ID_non_outlying = ID_non_outlying,
                                 Q_d1 = Q_d1,
                                 Q_d3 = Q_d3,
                                 IQR_d = IQR_d,
                                 Fvalue = Fvalue)

      ID_non_outlying = lst[['ID_non_outlying']]
      ID_shape_outlier = lst[['ID_shape_outlier']]

      # Manage low MEI data
      lst = manage_low_MEI_data(fData = mfData$fDList[[iDim]],
                                ID_non_outlying_Low_MEI = ID_non_outlying_Low_MEI,
                                ID_shape_outlier = ID_shape_outlier,
                                ID_non_outlying = ID_non_outlying,
                                Q_d1 = Q_d1,
                                Q_d3 = Q_d3,
                                IQR_d = IQR_d,
                                Fvalue = Fvalue)

      ID_non_outlying = lst[['ID_non_outlying']]
      ID_shape_outlier = lst[['ID_shape_outlier']]
    }
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


manage_high_MEI_data = function(fData,
                                ID_non_outlying_High_MEI,
                                ID_shape_outlier,
                                ID_non_outlying,
                                Q_d1,
                                Q_d3,
                                IQR_d,
                                Fvalue=1.5)
{
  N = fData$N

  a_0_2 = -2 / ( N * ( N - 1 ) )
  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  stopifnot( length(intersect(ID_shape_outlier, ID_non_outlying)) == 0 )

  # Managing High MEI data
  obs = min_diffs = NULL
  mins = data.frame(min_diffs = apply(fData$values, 2, compute_min_diff_min),
                    obs = apply(fData$values, 2, which.min))

  min_diff_min = mins %>%
    dplyr::filter(obs %in% ID_non_outlying_High_MEI ) %>%
    dplyr::group_by(obs) %>%
    dplyr::summarize(min_diff_min = min(min_diffs))

  ID_to_check = as.list(min_diff_min)[['obs']]

  aux_function_MBD = function( ID )(
    MBD( rbind( fData$values[ - ID, ],
                Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )

  aux_function_MEI = function( ID )(
    MEI( rbind( fData$values[ - ID, ],
                Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )


  if( nrow( min_diff_min ) > 0 )
  {
    Data_tilde = toRowMatrixForm( fData$values[ ID_to_check, ] +
                                    as.list(min_diff_min)[['min_diff_min']])

    MBD_curr = sapply( ID_to_check, aux_function_MBD )

    MEI_curr = sapply( ID_to_check, aux_function_MEI )

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    ID_out_extra = ID_to_check[ which( d_curr >= Q_d3 + Fvalue * IQR_d ) ]

    ID_shape_outlier = c( ID_shape_outlier, ID_out_extra )
    ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra )
  }

  return( list(ID_shape_outlier = ID_shape_outlier,
               ID_non_outlying = ID_non_outlying) )

}

manage_low_MEI_data = function(fData,
                               ID_non_outlying_Low_MEI,
                               ID_shape_outlier,
                               ID_non_outlying,
                               Q_d1,
                               Q_d3,
                               IQR_d,
                               Fvalue = 1.5)
{
  N = fData$N

  a_0_2 = -2 / ( N * ( N - 1 ) )
  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  stopifnot( length(intersect(ID_shape_outlier, ID_non_outlying)) == 0 )

  # Managing Low MEI data
  obs = max_diffs = NULL
  maxs = data.frame(max_diffs = apply(fData$values, 2, compute_max_diff_max),
                    obs = apply(fData$values, 2, which.max))

  max_diff_max = maxs %>%
    dplyr::filter(obs %in% ID_non_outlying_Low_MEI ) %>%
    dplyr::group_by(obs) %>%
    dplyr::summarize(max_diff_max = max(max_diffs))

  ID_to_check = as.list(max_diff_max)[['obs']]

  aux_function_MBD = function( ID )(
    MBD( rbind( fData$values[ - ID, ],
                Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )

  aux_function_MEI = function( ID )(
    MEI( rbind( fData$values[ - ID, ],
                Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )

  if( nrow( max_diff_max ) > 0 )
  {
    Data_tilde = toRowMatrixForm( fData$values[ ID_to_check, ] +
                                    as.list(max_diff_max)[['max_diff_max']])

    MBD_curr = sapply( ID_to_check, aux_function_MBD )

    MEI_curr = sapply( ID_to_check, aux_function_MEI )

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    ID_out_extra = ID_to_check[ which( d_curr >= Q_d3 + Fvalue * IQR_d ) ]

    ID_shape_outlier = c( ID_shape_outlier, ID_out_extra )
    ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra )
  }

  return( list(ID_shape_outlier = ID_shape_outlier,
               ID_non_outlying = ID_non_outlying) )

}


compute_min_diff_min = function(x)
{
  return(diff(sort(x, partial=c(1,2), decreasing=FALSE)[c(1,2)]))
}

compute_max_diff_max = function(x)
{
  return(diff(-sort(-x, partial=c(1,2), decreasing=FALSE)[c(1,2)]))
}

