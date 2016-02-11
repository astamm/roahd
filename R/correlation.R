

#'
#' \code{max.fData} maximum of an univariate functional dataset
#'
#'  It computes the maximum value of each element of the functional dataset,
#'  and optionally returns also the value of the grid where they are fulfilled
#'
#' @param fData the functional dataset containing elements whose maxima have to
#' be computed
#' @param ... additional parameters to pass to max
#' @param which logical flag specifying whether the grid values where maxima are
#' fulfilled have to be returned too
#'
max.fData = function( fData, ..., which = FALSE )
{

  if( which )
  {
    return( data.frame( grid = fData$t0 + ( apply( fData$values, 1, which.max ) - 1 ) * fData$h,
                        value = apply( fData$values, 1, max ) ) )
  } else {
    return( values = apply( fData$values, 1, max ) )
  }

}


#'
#' \code{min.fData} minimum of an univariate functional dataset
#'
#'  It computes the minimum value of each element of the functional dataset,
#'  and optionally returns also the value of the grid where they are fulfilled
#'
#' @param fData the functional dataset containing elements whose minima have to
#' be computed
#' @param ... additional parameters to pass to max
#' @param which logical flag specifying whether the grid values where minima are
#' fulfilled have to be returned too
#'
min.fData = function( fData, ..., which = FALSE )
{
  if( which )
  {
    return( data.frame( grid = fData$t0 + ( apply( fData$values, 1, which.min ) - 1 ) * fData$h,
                        value = apply( fData$values, 1, min ) ) )
  } else {
    return( values = apply( fData$values, 1, min ) )
  }

}


#'  Maximum order relation among univariate functional data
#'
#'  It implements the maximum order relation among univariate functional data,
#'  that is the pre-order relation obtained by comparing the maxima of two
#'  different functional data.
#'
#'  It accepts two fData objects, that can have either same sample size (number
#'  of elemnets), or just one element. In the former case the comparison is made
#'  element-wise between corresponding observations of the two datasets; in the
#'  latter case, the dataset with just one observation is used cyclically
#'  against the other one (if need be).
#'
#' @param fData the first univariate functional dataset containing elements to
#' be compared
#' @param gData the second univariate functional dataset containing elements to
#' be compared
#'
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

  return( max( fData ) - max( gData ) <= 0 )
}

#'
#' Area under curve of elements of a univariate functional dataset
#'
#'  It computes the area under the curve of elements of a univariate functional
#'  dataset. A plain trapezoidal rule is employed to compute the integral.
#'
#' @param fData the functional dataset containing elements whose areas under the
#' curve have to be computed
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

#'  Area under curve order relation among univariate functional data
#'
#'  It implements the area under curve order relation among univariate f
#'  unctional data, that is the pre-order relation obtained by comparing the
#'  area under curve maxima of two different functional data.
#'
#'  It accepts two fData objects, that can have either same sample size (number
#'  of elemnets), or just one element. In the former case the comparison is made
#'  element-wise between corresponding observations of the two datasets; in the
#'  latter case, the dataset with just one observation is used cyclically
#'  against the other one (if need be).
#'
#' @param fData the first univariate functional dataset containing elements to
#' be compared
#' @param gData the second univariate functional dataset containing elements to
#' be compared
#'
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

#'  Kendall's tau correlation coefficient for a bivariate functional dataset
#'
#'  It implements the computation of the Kendall's tau correlation coefficient
#'  for bivariate functional data, with the prescribed order relation (either
#'  max or area-under-curve)
#'
#' @param mfD the bivariate functional dataset whose Kendall's tau
#' coefficient must be computed
#' @param ordering the ordering relation to use on functional observations,
#' either "max" for the maximum relation or "area" for the area under the curve
#' relation
#'
cor_kendall = function( mfD, ordering = 'max' )
{
    if( mfD$L != 2 )
    {
      stop( 'Error in cor_kendall__var: only bivariate data are supported for now' )
    }

  N = mfD$N

  if( ordering == 'max' )
  {
    R = matrix( c( max( mfD$fDList[[ 1 ]] ),
                   max( mfD$fDList[[ 2 ]] ) ),
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

#'  Spearman's correlation coefficient for a bivariate functional dataset
#'
#'  It implements the computation of the Spearmans's correlation coefficient
#'  for bivariate functional data, with the prescribed order relation (either
#'  MEI or MHI)
#'
#' @param mfD the bivariate functional dataset whose Kendall's tau
#' coefficient must be computed
#' @param ordering the ordering relation to use on functional observations,
#' either "MEI" or "MHI"
#' @param ... additional parameters to be passed to rank function
#'
cor_spearman = function( mfD, ordering = 'MEI', ... )
{
  if( mfD$L != 2 )
  {
    stop( ' Error in cor_spearman: only bivariate data are supported for now')
  }

  if( ordering == 'MEI' )
  {
    rk_1 = rank( MEI( mfD$fDList[[ 1 ]]$values ), ... )
    rk_2 = rank( MEI( mfD$fDList[[ 2 ]]$values ), ... )
  } else if( ordering == 'MHI' )
  {
    rk_1 = rank( MHI( mfD$fDList[[ 1 ]]$values ), ... )
    rk_2 = rank( MHI( mfD$fDList[[ 2 ]]$values ), ... )
  }
  return( cor( rk_1, rk_2, method = 'pearson' ) )
}




