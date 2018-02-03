#' Maxima of a univariate functional dataset
#'
#' This function computes the maximum value of each element of a
#' univariate functional dataset, optionally returing also the value of the
#' grid where they are fulfilled.
#'
#' @param fData the functional dataset containing elements whose maxima have to
#' be computed, in form of \code{fData} object.
#' @param ... additional parameters.
#' @param which logical flag specifying whether the grid values where maxima are
#' fulfilled have to be returned too.
#'
#' @return If \code{which = FALSE}, the function returns a vector containing the
#' maxima for each element of the functional dataset; if \code{which = TRUE},
#' the function returns a \code{data.frame} whose field \code{value} contains
#' the values of maxima, and \code{grid} contains the grid points where maxima
#' are reached.
#'
#' @examples
#'
#' P = 1e3
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' Data = matrix( c( 1 * grid,
#'                   2 *  grid,
#'                   3 * ( 0.5 - abs( grid - 0.5 ) ) ),
#'                nrow = 3, ncol = P, byrow = TRUE )
#'
#' fD = fData( grid, Data )
#'
#' maxima( fD, which = TRUE )
#'
#' @seealso \code{\link{minima}}
#'
#' @export
#'
maxima = function( fData, ..., which = FALSE )
{
  if( which )
  {
    return( data.frame( grid = fData$t0 + ( apply( fData$values, 1,
                                                   which.max ) - 1 ) * fData$h,
                        value = apply( fData$values, 1, max ) ) )
  } else {
    return( values = apply( fData$values, 1, max ) )
  }

}


#' Minima of a univariate functional dataset
#'
#' This function computes computes the minimum value of each element of a
#' univariate functional dataset, optionally returing also the value of the
#' grid where they are fulfilled.
#'
#' @param fData the functional dataset containing elements whose minima have to
#' be computed, in form of \code{fData} object.
#' @param ... additional parameters.
#' @param which logical flag specifying whether the grid values where minima are
#' fulfilled have to be returned too.
#'
#' @return If \code{which = FALSE}, the function returns a vector containing the
#' minima for each element of the functional dataset; if \code{which = TRUE},
#' the function returns a \code{data.frame} whose field \code{value} contains
#' the values of minima, and \code{grid} contains the grid points where minima
#' are reached.
#'
#' @examples
#'
#' P = 1e3
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' Data = matrix( c( 1 * grid,
#'                   2 *  grid,
#'                   3 * ( 0.5 - abs( grid - 0.5 ) ) ),
#'                nrow = 3, ncol = P, byrow = TRUE )
#'
#' fD = fData( grid, Data )
#'
#' minima( fD, which = TRUE )
#'
#' @seealso \code{\link{maxima}}
#'
#' @export
#'
minima = function( fData, ..., which = FALSE )
{
  if( which )
  {
    return( data.frame( grid = fData$t0 + ( apply( fData$values, 1, which.min ) - 1 ) * fData$h,
                        value = apply( fData$values, 1, min ) ) )
  } else {
    return( values = apply( fData$values, 1, min ) )
  }

}


#'  Maximum order relation between univariate functional data
#'
#'  This function implements an order relation between univariate functional
#'  data based on the maximum relation, that is to say a pre-order relation
#'  obtained by comparing the maxima of two different functional data.
#'
#'  Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)}
#'  and another functional dataset \eqn{Y_1(t),} \eqn{Y_2(t), \ldots, Y_M(t)}
#'  defined over the same compact interval \eqn{I=[a,b]}, the function computes
#'  the maxima in both the datasets, and checks whether the first ones are lower
#'  or equal than the second ones.
#'
#'  By default the function tries to compare each \eqn{X_i(t)} with the
#'  corresponding \eqn{Y_i(t)}, thus assuming \eqn{N=M}, but when either \eqn{N=1}
#'  or \eqn{M=1}, the comparison is carried out cycling over the dataset with
#'  fewer elements. In all the other cases (\eqn{N\neq M,} and either
#'  \eqn{N \neq 1} or \eqn{M \neq 1}) the function stops.
#'
#' @param fData the first univariate functional dataset containing elements to
#' be compared, in form of \code{fData} object.
#' @param gData the second univariate functional dataset containing elements to
#' be compared, in form of \code{fData} object.
#'
#' @return
#' The function returns a logical vector of length \eqn{\max(N,M)} containing the
#' value of the predicate for all the corresponding elements.
#'
#' @references
#'
#' Valencia, D., Romo, J. and Lillo, R. (2015). A Kendall correlation
#' coefficient for functional dependence,
#' \emph{Universidad Carlos III de Madrid technical report},
#' \code{http://EconPapers.repec.org/RePEc:cte:wsrepe:ws133228}.
#'
#'
#'
#' @seealso \code{\link{maxima}}, \code{\link{minima}}, \code{\link{fData}},
#' \code{\link{area_ordered}}
#'
#' @examples
#'
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' Data_1 = matrix( c( 1 * grid,
#'                     2 *  grid ),
#'                  nrow = 2, ncol = P, byrow = TRUE )
#'
#' Data_2 = matrix( 3 * ( 0.5 - abs( grid - 0.5 ) ),
#'                  nrow = 1, byrow = TRUE )
#'
#' Data_3 = rbind( Data_1, Data_1 )
#'
#'
#' fD_1 = fData( grid, Data_1 )
#' fD_2 = fData( grid, Data_2 )
#' fD_3 = fData( grid, Data_3 )
#'
#' max_ordered( fD_1, fD_2 )
#'
#' max_ordered( fD_2, fD_3 )
#'
#' @export
max_ordered = function( fData, gData )
{
  if( fData$P != gData$P ||
      fData$h != gData$h ||
      fData$t0 != gData$t0 ||
      fData$tP != gData$tP  )
  {
    stop( ' Error in max_ordered: provided fData objects have mismatching time grids')
  }

  if( fData$N != gData$N )
  {
    fData$N == 1 || gData$N == 1 ||
      stop( ' Error  in max_ordered: you must provide equally sized fData
            objects or at least one of them with 1 observation.')
  }

  return( maxima( fData ) - maxima( gData ) <= 0 )
}

#'
#' Area under curve of elements of univariate functional data
#'
#' This method computes the (signed) area under the curve of elements of a
#' univariate functional dataset, namely, their integral.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]} and observed on an evenly sapced
#' 1D grid \eqn{[a = t_0, t_1, \ldots, t_{P-1} = b \subset I}, the function
#' computes:}
#'
#' \deqn{ \sum_{i=1}^{P-2} \frac{X(t_{i+1}) - X(t_{i-1})}{2} h \approx
#' \int_a^b X(t) dt,}
#'
#' where \eqn{h = t_1 - t_0}.
#'
#' @param fData the functional dataset containing elements whose areas under the
#' curve have to be computed, in form of \code{fData} object.
#'
#' @return The function returns a numeric vector containing the values of areas
#' under the curve for all the elements of the functional dataset \code{fData}.
#'
#' @seealso \code{\link{area_ordered}}, \code{\link{fData}}
#'
#' @examples
#'
#' P = 1e3
#' grid = seq( 0, 1, length.out = P )
#'
#' fD = fData( grid,
#'             matrix( c( sin( 2 * pi * grid ),
#'                        cos( 2 * pi * grid ),
#'                        4 * grid * ( 1 - grid ) ),
#'                     nrow = 3, ncol = P, byrow = TRUE ) )
#' plot( fD )
#'
#' area_under_curve( fD )
#'
#' @export
#'
area_under_curve = function( fData)
{
  if( fData$N > 1 )
  {
    return( rowSums( ( fData$values[ , - 1 ] +
                         fData$values[ , - fData$P ] ) / 2 * fData$h ) )
  } else {
    return( sum( ( fData$values[ , - 1 ] +
                         fData$values[ , - fData$P ] ) / 2 * fData$h ) )
  }
}

#'  Area-under-curve order relation between univariate functional data
#'
#'  This function implements an order relation between univariate functional
#'  data based on the area-under-curve relation, that is to say a pre-order
#'  relation obtained by comparing the area-under-curve of two
#'  different functional data.
#'
#'  Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)}
#'  and another functional dataset \eqn{Y_1(t),} \eqn{Y_2(t), \ldots, Y_M(t)}
#'  defined over the same compact interval \eqn{I=[a,b]}, the function computes
#'  the area-under-curve (namely, the integral) in both the datasets, and checks
#'  whether the first ones are lower or equal than the second ones.
#'
#'  By default the function tries to compare each \eqn{X_i(t)} with the
#'  corresponding \eqn{Y_i(t)}, thus assuming \eqn{N=M}, but when either \eqn{N=1}
#'  or \eqn{M=1}, the comparison is carried out cycling over the dataset with
#'  fewer elements. In all the other cases (\eqn{N\neq M,} and either
#'  \eqn{N \neq 1} or \eqn{M \neq 1}) the function stops.
#'
#' @param fData the first univariate functional dataset containing elements to
#' be compared, in form of \code{fData} object.
#' @param gData the second univariate functional dataset containing elements to
#' be compared , in form of \code{fData} object.
#'
#' @return
#' The function returns a logical vector of length \eqn{\max(N,M)} containing the
#' value of the predicate for all the corresponding elements.
#'
#' @references
#'
#' Valencia, D., Romo, J. and Lillo, R. (2015). A Kendall correlation
#' coefficient for functional dependence,
#' \emph{Universidad Carlos III de Madrid technical report},
#' \code{http://EconPapers.repec.org/RePEc:cte:wsrepe:ws133228}.
#'
#'
#' @seealso \code{\link{area_under_curve}}, \code{\link{fData}}
#'
#' @examples
#'
#' P = 1e3
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' Data_1 = matrix( c( 1 * grid,
#'                     2 *  grid ),
#'                  nrow = 2, ncol = P, byrow = TRUE )
#'
#' Data_2 = matrix( 3 * ( 0.5 - abs( grid - 0.5 ) ),
#'                  nrow = 1, byrow = TRUE )
#'
#' Data_3 = rbind( Data_1, Data_1 )
#'
#'
#' fD_1 = fData( grid, Data_1 )
#' fD_2 = fData( grid, Data_2 )
#' fD_3 = fData( grid, Data_3 )
#'
#' area_ordered( fD_1, fD_2 )
#'
#' area_ordered( fD_2, fD_3 )
#'
#' @export
area_ordered = function( fData, gData )
{
  if( fData$P != gData$P ||
      fData$h != gData$h ||
      fData$t0 != gData$t0 ||
      fData$tP != gData$tP  )
  {
    stop( ' Error in area_ordered: provided fData objects have mismatching time grids')
  }

  if( fData$N != gData$N )
  {
    if( fData$N == 1 )
    {
      return( area_under_curve( gData - fData$values ) >= 0 )
    } else if( gData$N == 1 )
    {
      return( area_under_curve( fData - gData$values ) <= 0 )
    } else {
      stop( ' Error  in area_ordered: you must provide equally sized fData
            objects or at least one of them with 1 observation.')
    }
  } else {
    return( area_under_curve( fData - gData ) <= 0 )
  }
}

# cor_kendall = function( mfD, ordering = 'max' )
# {
#   if( mfD$L != 2 )
#   {
#     stop( ' Error in cor_kendall: only bivariate data are supported for now.')
#   }
#
#   N = mfD$N
#
#   if( ordering == 'area' )
#   {
#     count_concordances = function( iObs )( sum( area_ordered( mfD$fDList[[ 1 ]][ iObs, ],
#                                                               mfD$fDList[[ 1 ]][ ( iObs + 1 ) : N, ] ) ==
#                                                   area_ordered( mfD$fDList[[ 2 ]][ iObs, ],
#                                                                 mfD$fDList[[ 2 ]][ ( iObs + 1 ) : N, ] ) ) )
#   } else if ( ordering == 'max' )
#   {
#     count_concordances = function( iObs )( sum( max_ordered( mfD$fDList[[ 1 ]][ iObs, ],
#                                                              mfD$fDList[[ 1 ]][ ( iObs + 1 ) : N, ] ) ==
#                                                   max_ordered( mfD$fDList[[ 2 ]][ iObs, ],
#                                                                mfD$fDList[[ 2 ]][ ( iObs + 1 ) : N, ] ) ) )
#   }else
#   {
#     stop( ' Error in cor_kendall: unsupported ordering relation')
#   }
#
#   return( ( 2 * sum( sapply( 1 : ( N - 1 ), count_concordances ) )  - ( N * ( N - 1 ) / 2 ) ) / ( N * ( N - 1 ) / 2 ) )
# }


#'  Kendall's tau correlation coefficient for bivariate functional data
#'
#' This function computes the Kendall's tau correlation coefficient for a
#' bivariate functional dataset, with either a max or area-under-curve order
#' order relation between univariate functional elements (components).
#'
#' Given a bivariate functional dataset, with first components \eqn{X_1(t),
#' X_2(t), \ldots, X_N(t)} and second components \eqn{Y_1(t), Y_2(t), \ldots,
#' Y_N(t)}, the function exploits either the order relation based on the maxima
#' or the area-under-curve relation to compare data and produce concordances and
#' discordances, that are then used to compute the tau coefficient.
#'
#' See the references for more details.
#'
#' @param mfD a bivariate functional dataset whose Kendall's tau
#' coefficient must be computed, in form of bivariate \code{mfData} object
#' (\code{mfD$L=2}).
#' @param ordering the ordering relation to use on functional observations,
#' either \code{"max"} for the maximum relation or \code{"area"} for the
#' area under the curve relation (default is \code{"max"}).
#'
#' @return The function returns the Kendall's tau correlation coefficient for
#' the bivariate dataset provided with \code{mfData}.
#'
#' @references
#'
#' Valencia, D., Romo, J. and Lillo, R. (2015). A Kendall correlation
#' coefficient for functional dependence,
#' \emph{Universidad Carlos III de Madrid technical report},
#' \code{http://EconPapers.repec.org/RePEc:cte:wsrepe:ws133228}.
#'
#'
#' @examples
#'
#' #### TOTALLY INDEPENDENT COMPONENTS
#' N = 2e2
#' P = 1e3
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' # Creating an exponential covariance function to simulate guassian data
#' Cov = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating (independent) gaussian functional data with given center and
#' # covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#'
#' # Using the simulated data as (independent) components of a bivariate functional
#' # dataset
#' mfD = mfData( grid, list( Data_1, Data_2 ) )
#'
#' # Correlation approx. zero (components were created independently)
#' cor_kendall( mfD, ordering = 'max' )
#'
#' # Correlation approx. zero (components were created independently)
#' cor_kendall( mfD, ordering = 'area' )
#'
#' #### TOTALLY DEPENDENT COMPONENTS
#'
#' # Nonlinear transform of first component
#' Data_3 = t( apply( Data_1, 1, exp ) )
#'
#' # Creating bivariate dataset starting from nonlinearly-dependent components
#' mfD = mfData( grid, list( Data_1, Data_3 ) )
#'
#' # Correlation very high (components are nonlinearly dependent)
#' cor_kendall( mfD, ordering = 'max' )
#'
#' # Correlation very high (components are nonlinearly dependent)
#' cor_kendall( mfD, ordering = 'area' )
#'
#' @seealso \code{\link{mfData}}, \code{\link{area_ordered}},
#' \code{\link{max_ordered}}
#'
#' @export
cor_kendall = function( mfD, ordering = 'max' )
{
    if( mfD$L != 2 )
    {
      stop( 'Error in cor_kendall__var: only bivariate data are supported for now' )
    }

  N = mfD$N

  if( ordering == 'max' )
  {
    R = matrix( c( maxima( mfD$fDList[[ 1 ]] ),
                   maxima( mfD$fDList[[ 2 ]] ) ),
                ncol = 2, nrow = N, byrow = FALSE )

  } else if( ordering == 'area' )
  {
    R = matrix( c( area_under_curve( mfD$fDList[[ 1 ]] ),
                   area_under_curve( mfD$fDList[[ 2 ]] ) ),
                ncol = 2, nrow = N, byrow = FALSE )

  } else{
    stop( ' Error in cor_kendall__var: unsupported ordering relation')
  }

  aux_function = function( iRow )( sum( abs( rowSums( sign( t( t( R[ ( iRow + 1 ) : N,  ] ) - R[ iRow, ] ) ) ) ) ) )

  return( ( sum( sapply( 1 : ( N - 2 ), aux_function ) ) +
                       abs( sum( sign( R[ N,  ] - R[ N - 1, ]  ) ) ) ) / ( N * ( N - 1 ) / 2 )  - 1 )
}

#' Spearman's correlation coefficient for multivariate functional data
#'
#' This function computes the Spearman's correlation coefficient for a
#' multivariate functional dataset, with either a Modified Epigraph Index (MEI) or
#' Modified Hypograph Index (MHI) ranking of univariate elments of data
#' components.
#'
#' Given a multivariate functional dataset, with first components \eqn{X^1_1(t),
#' X^1_2(t), \ldots, X^1_N(t)}, second components \eqn{X^2_1(t), X^2_2(t), \ldots,
#' X^2_N(t)}, etc., the function exploits either the MEI or MHI to compute the matrix of
#' Spearman's correlation coefficients. Such matrix is symmetrical and has ones on the
#' diagonal. The entry (i, j) represents the Spearman correlation coefficient between
#' curves of component i and j.
#'
#' See the references for more details.
#'
#' @param mfD a multivariate functional dataset whose Spearman's correlation
#' coefficient must be computed, in form of multivariate \code{mfData} object.
#' @param ordering the ordering relation to use on functional observations,
#' either \code{"MEI"} for MEI or \code{"MHI"} for MHI (default is \code{"MEI"}).
#'
#' @return If the original dataset is bivariate, the function returns only the
#' scalar value of the correlation coefficient for the two components.
#' When the number of components is L >2, it returns the whole matrix of
#' Spearman's correlation coefficients for all the components.
#'
#' @references
#'
#' Valencia, D., Romo, J. and Lillo, R. (2015). Spearman coefficient for
#' functions, \emph{Universidad Carlos III de Madrid technical report},
#' \code{http://EconPapers.repec.org/RePEc:cte:wsrepe:ws133329}.
#'
#' @examples
#'
#' #### TOTALLY INDEPENDENT COMPONENTS
#'
#' N = 2e2
#' P = 1e3
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' # Creating an exponential covariance function to simulate guassian data
#' Cov = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating (independent) gaussian functional data with given center and
#' # covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#'
#' # Using the simulated data as (independent) components of a bivariate functional
#' # dataset
#' mfD = mfData( grid, list( Data_1, Data_2 ) )
#'
#' # Correlation approx. zero (components were created independently)
#' cor_spearman( mfD, ordering = 'MEI' )
#'
#' # Correlation approx. zero (components were created independently)
#' cor_spearman( mfD, ordering = 'MHI' )
#'
#' #### TOTALLY DEPENDENT COMPONENTS
#'
#' # Nonlinear transform of first component
#' Data_3 = t( apply( Data_1, 1, exp ) )
#'
#' # Creating bivariate dataset starting from nonlinearly-dependent components
#' mfD = mfData( grid, list( Data_1, Data_3 ) )
#'
#' # Correlation very high (components are nonlinearly dependent)
#' cor_spearman( mfD, ordering = 'MEI' )
#'
#' # Correlation very high (components are nonlinearly dependent)
#' cor_spearman( mfD, ordering = 'MHI' )
#'
#' @seealso \code{\link{mfData}}, \code{\link{MEI}}, \code{\link{MHI}}
#'
#' @export
cor_spearman = function( mfD, ordering = 'MEI' )
{
  # if( mfD$L != 2 )
  # {
  #   stop( ' Error in cor_spearman: only bivariate data are supported for now')
  # }

  if( ordering == 'MEI' )
  {
    rks = sapply(mfD$fDList, MEI)

    # rk_1 = MEI( mfD$fDList[[ 1 ]]$values )
    # rk_2 = MEI( mfD$fDList[[ 2 ]]$values )
  } else if( ordering == 'MHI' )
  {
    rks = sapply(mfD$fDList, MHI)

    # rk_1 = MHI( mfD$fDList[[ 1 ]]$values )
    # rk_2 = MHI( mfD$fDList[[ 2 ]]$values )
  }

  cor_output = cor(rks, method='pearson')

  if( mfD$L == 2 )
  {
    return(cor_output[1,2])
  } else{
    return(cor_output)
  }
  # return( cor( rk_1, rk_2, method = 'pearson' ) )
}

#' Bootstrap Spearman's correlation coefficient for multivariate functional data
#'
#' This function computes the bootstrap estimates of standard error and bias of the Spearman's
#' correlation coefficient for a multivariate functional dataset.
#'
#' Given a multivariate functional dataset \eqn{X_1^(i), \ldots, X_n^(i)}, \eqn{i=0, \ldots, L}
#' defined over the grid \eqn{I = t_0, \ldots, t_P}, having components \eqn{i=1, \ldots, L}, and a
#' chosen ordering strategy (MEI or MHI), the function computes the matrix of Speraman's correlation
#' indexes of the dataset's components, as well as their bias and standard deviation estimates
#' through a specified number of bootstrap iterations (bias and standard error are updated with
#' on-line formulas).
#'
#' @param mfD a multivariate functional dataset whose Spearman's correlation
#' coefficient must be computed, in form of multivariate \code{mfData} object.
#' @param ordering the ordering relation to use on functional observations,
#' either \code{"MEI"} for MEI or \code{"MHI"} for MHI (default is \code{"MEI"}).
#' @param bootstrap_iterations the number of bootstrap iterations to be used for estimation of bias and
#' standard error.
#' @param verbose a logical flag specifying whether to log information on the estimation progress.
#'
#'
#' @returns a list of three elements: \code{mean}, the mean of the matrix of correlation coefficients;
#' \code{bias}, a matrix containing the estimated bias (mean - point estimate of correlation coefficients);
#' \code{sd}, a matrix containing the estiated standard deviation of the coefficients' matrix. In case
#' the multivariate functional dataset has only two components, the return type is scalar and not matrix.
#'
#' @seealso \code{\link{cor_spearman}}, \code{\link{mfData}}
#'
#' @export
#'
cor_spearman_accuracy = function(mfD, ordering='MEI', bootstrap_iterations=1000,
                                 verbose=FALSE)
{
  if( verbose ){
    printout_iters = ceiling(seq(1, bootstrap_iterations, length.out=11))
  }

  if( mfD$L == 2 ){
    cor_curr = 0
  } else{
    cor_curr = matrix(0, nrow=mfD$L, ncol=mfD$L)
  }

  # Initialising data structures
  cor_mean = cor_curr
  # diag(cor_mean) = 1
  cor_sqd = cor_curr

  for( iBoot in 1:bootstrap_iterations)
  {
    if ((verbose) && (any(iBoot == printout_iters))){
      message(paste0(round(iBoot/bootstrap_iterations, 2) * 100,
                     '% bootstrap Spearman`s correlation'))
    }

    w = sample.int(mfD$N, size=mfD$N, replace=TRUE, prob=NULL)

    cor_curr = cor_spearman(mfD[w,])

    # The order here is important to use the online updates
    increment = (cor_curr - cor_mean)
    cor_mean = cor_mean + increment / iBoot
    cor_sqd = cor_sqd + increment * ( cor_curr - cor_mean )
  }

  return( list(mean = cor_mean,
               bias = cor_mean - cor_spearman(mfD, ordering = ordering),
               sd = sqrt(cor_sqd / ( bootstrap_iterations - 1 ))))
}

