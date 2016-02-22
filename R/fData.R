
#' \code{fData} class for univariate functional data.
#'
#'  This function implements a constructor for elements of \code{S3} class
#'  \code{fData}, aimed at implementing a representation of a functional
#'  dataset.
#'
#'  The functional dataset is represented as a collection of measurement of the
#'  observations on an evenly spaced, 1D grid of discrete points (representing,
#'  e.g. time).
#'
#' @param grid the evenly spaced grid over which the functional observations are
#' measured. It must be a numeric vector of length \code{P}.
#' @param values the values of the observations in the functional dataset,
#' prodived in form of a 2D data structure (e.g. matrix or array) having as
#' rows the observations and as columns their measurements over the 1D grid of
#' length \code{P} specified in \code{grid}.
#'
#' @export
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

#' Specialised method to plot univariate functional data..
#'
#' @param x the univariate functional data object
#' @param ... additional graphical parameters to be used in plotting functions
#'
plot.fData = function( x, ... )
{
  plot_fData_default( x, ...  )
}

#' Default method to plot univariate functional data.
#'
#' @param x the univariate functional data object
#' @param lty lty graphical parameter to be used in plotting functions
#' @param col colors to be used in plotting functions
#' @param xlab the x label to add to each plot window
#' @param ylab either a string or a list of stirngs to be used as y labels for
#' each window in the plot panel, default is null
#' @param main either a string or a list of stirngs to be used as main title for
#' each window in the plot panel, default is null
#' @param ... additional graphical parameters to be used in plotting functions
plot_fData_default = function( x,
                               type = 'l', lty = 1,
                               col = fDColorPalette( min( c( x$N,
                                                             30 + x$N %% 30 ) ) ),
                               xlab = '', ylab = '', main = '',
                               ... )
{
  matplot( seq( x$t0, x$tP, length.out = x$P ),
           t( x$values ), type = type, lty = lty,
           col = col, xlab = xlab, ylab = ylab, main = main, ... )
}



#'
#' \code{mfData} generic class for (multivariate) functional data.
#'
#' @param grid the (evenly spaced) grid over which the functional data is defined
#' @param Data_list a list containing the time-by-time values of functional
#' dataset, where each node of the list (possibly named) represents a particular
#' dimension of the multivariate dataset.
#'
mfData = function( grid, Data_list )
{
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in mfData: you provided an unevenly spaced grid')

  dimMatrix = sapply( Data_list, dim )

  if( unique( apply( dimMatrix, 1,
                     function( x )( length( unique( x ) ) ) ) ) != 1 )
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

  for( nmL in nms )
  {
    fDList[[ nmL ]] = fData( grid, toRowMatrixForm( Data_list[[ nmL ]] ) )
  }

  return( structure( list( N = fDList[[ 1 ]]$N,
                           P = fDList[[ 1 ]]$P,
                           L = L,
                           t0 = fDList[[ 1 ]]$t0,
                           tP = fDList[[ 1 ]]$tP,
                           fDList = fDList ),
                     class = c( 'mfData' ) ) )
}

#' Specialised method to plot multivariate functional data.
#'
#' @param x the multivariate functional data object
#' @param ... additional graphical parameters to be used in plotting functions
#'
plot.mfData = function( x, ... )
{
  mfrow_rows = floor( sqrt( x$L ) )
  mfrow_cols = ceiling( x$L / floor( sqrt( x$L ) ) )

  par( mfrow = c( mfrow_rows, mfrow_cols ) )

  plot_mfData_default( x, ... )

}

#' Default method to plot multivariate functional data.
#'
#' @param x the multivariate functional data object
#' @param lty lty graphical parameter to be used in plotting functions
#' @param col colors to be used in plotting functions
#' @param xlab the x label to add to each plot window
#' @param ylab either a string or a list of stirngs to be used as y labels for
#' each window in the plot panel, default is null
#' @param main either a string or a list of stirngs to be used as main title for
#' each window in the plot panel, default is null
#' @param ... additional graphical parameters to be used in plotting functions
plot_mfData_default = function( x,
                                type = 'l',lty = 1,
                                col = fDColorPalette( min( c( x$N,
                                                              30 + x$N %% 30 ) ) ),
                                xlab = NULL, ylab = NULL, main = NULL,
                                ... )
{
  if( ! is.null( ylab ) )
  {
    if( length( ylab ) == 1 )
    {
      ylab = rep( ylab, x$L )
    } else if( length( ylab ) != x$L )
    {
      stop( 'Error in plot_mfData_default: you specified a wrong number of y
            labels' )
    }
  } else {
    ylab = rep( list( '' ), x$L )
  }

  if( ! is.null( main ) )
  {
    if( length( main ) == 1 )
    {
      main = rep( main, x$L )
    } else if( length( main ) != x$L )
    {
      stop( 'Error in plot_mfData_default: you specified a wrong number of
            subtitles' )
    }
  } else {
    main = rep( list( '' ), x$L )
  }

  plot_aux = function( i, ... )( plot.fData( x$fDList[[ i ]],
                                             ylab = ylab[[ i ]],
                                             main = main[[ i ]],
                                             col = col,
                                             lty = lty,
                                             ... ) )

  invisible( sapply( 1 : x$L, plot_aux, ... ) )

}


#' Operator \code{+} and \code{-} for \code{fData} objects
#'
#' These methods provide operators \code{+} and \code{-} to perform sums
#' or differences between an \code{fData} object and either another
#' \code{fData} object or other compliant data structures, like matrices or
#' vectors or arrays, representing the pointwise measurements of the second
#' term of the  sum.
#'
#' If the second term of the operation is an \code{fData} object, it must be
#' defined over the same grid as the first.
#'
#' @param fD the univariate \code{fData} object.
#' @param A either an \code{fData} object, defined on the very same grid of
#' \code{fD}, or a 1D data structure (such as 1D array or raw
#' numeric vector), or a 2D data structure (such as 2D array or raw numeric
#' matrix ), that specifies the second term of the sum.
#' In case of a 1D data structure, the sum is performed element-wise between
#' each element of  \code{fD} and \code{A}, and \code{A} must have length
#' \code{P}, size of \code{fD}'s grid.
#' In case of a 2D data structure, the sum is performed element-wise between
#' corresponding elements of \code{fD} and \code{A}'s rows. In this case,
#' \code{A} must have \code{P} columns, as the size of \code{fD}'s grid.
#'
#' @name plus-.fData
#'
#' @return The function returns an \code{fData} object, whose function values
#' have undergone the sum/difference.
#'
#'
NULL


#' @rdname plus-.fData
#'
#' @examples
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = 21, ncol = 10, byrow = TRUE ) )
#' fD + 1 : 10
#'
#' fD + array( 1, dim = c( 1, 10 ) )
#'
#' fD + fD
#'
#' @export
"+.fData" = function( fD, A )
{
  if( class( A ) == 'fData' )
  {
    if( fD$t0 != A$t0 || fD$tP != A$tP || fD$h != A$h || fD$P != A$P )
    {
      stop( 'Error in +.fData: functional data defined over
            mismatching intervals' )
    }

    if( A$N == 1 )
    {
      fD$values = t( t( fD$values ) + as.vector( A$values ) )
    } else {
      fD$values = fD$values + A$values
    }

  } else if( is.null( dim( A ) ) || nrow( A ) == 1 ){

    A = as.vector( A )

    if( length( A ) != fD$P )
    {
      stop( 'Error in +.fData: mismatching arguments')
    }

    fD$values = t( t( fD$values ) +
                     as.vector( rep( A, fD$P / length( A ) ) ) )

  } else if( is.matrix( A ) ) {

    if( ncol( A ) != fD$P || nrow( A ) != fD$N )
    {
      stop('Error in +.fData: mismatching arguments' )
    }

    fD$values = fD$values + A
  }

  return( fD )
}

#' @rdname plus-.fData
#'
#' @examples
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = 21, ncol = 10, byrow = TRUE ) )
#' fD - 2 : 11
#'
#' fD - array( 1, dim = c( 1, 10 ) )
#'
#' fD - fD
#'
#' @export
#'
"-.fData" = function( fD, A )
{
  if( class( A ) == 'fData' )
  {
    if( fD$t0 != A$t0 || fD$tP != A$tP || fD$h != A$h || fD$P != A$P )
    {
      stop( 'Error in -.fData: functional data defined over
            mismatching intervals' )
    }

    if( A$N == 1 )
    {
      fD$values = t( t( fD$values ) - as.vector( A$values ) )
    } else {
      fD$values = fD$values - A$values
    }

  } else if( is.null( dim( A ) ) || nrow( A ) == 1 ){

    A = as.vector( A )

    if( length( A ) != fD$P )
    {
      stop( 'Error in -.fData: mismatching arguments')
    }

    fD$values = t( t( fD$values ) -
                     as.vector( rep( A, fD$P / length( A ) ) ) )

  } else if( is.matrix( A ) ) {

    if( ncol( A ) != fD$P || nrow( A ) != fD$N )
    {
      stop('Error in -.fData: mismatching arguments' )
    }

    fD$values = fD$values - A
  }

  return( fD )
}

#' Operator \code{*} and \code{/} for \code{fData} objects
#'
#' These methods provide operators \code{*} and \code{/} to perform products
#' or divisions between an \code{fData} object and either a number or a
#' compliant 1D data structure, like numeric vector, array or
#' matrix. The operation is computed by performing the element-wise product
#' or division between \code{fD}'s observations and the provided value(s).
#'
#' If the second argument is a 1D data structure, it must have length \code{N}
#' equal to the number of observations in \code{fD}.
#'
#'
#' @param fD the univariate \code{fData} object.
#' @param a either a single number or a 1D data structure (such as numeric
#' raw vector, matrix or array) specifying the factor(s) to use in the
#' multiplication/division of \code{fD} elements' values.
#' In the latter case, each factor is used with the corresponding element in
#' \code{fD}, hence a must have length \code{N}, number of observations in
#' \code{fD}.
#'
#' @name times-.fData
#'
#' @return The function returns an \code{fData} object, whose function values
#' have undergone the product/division.
#'
NULL

#' @rdname times-.fData
#'
#' @examples
#'
#' N = 11
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = N, ncol = 10, byrow = TRUE ) )
#' fD * 2
#'
#' fD * seq( 1, N )
#'
#' @export
#'
"*.fData" = function( fD, a )
{
  if( ! 1 %in% dim( as.matrix( a ) ) )
  {
      stop( 'Error in *.fData: dimensions are not compliant' )
  }

  fD$values = fD$values * as.numeric( a )

  return( fD )
}

#' @rdname times-.fData
#'
#' @examples
#'
#' N = 11
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = N, ncol = 10, byrow = TRUE ) )
#' fD / 2
#'
#' fD / rep( 10, N )
#'
#' @export
"/.fData" = function( fD, a )
{
  if( ! 1 %in% dim( as.matrix( a ) ) )
  {
    stop( 'Error in *.fData: dimensions are not compliant' )
  }

  fD$values = fD$values / as.numeric( a )


  return( fD )
}

#'
#' \code{mean_fData} method to compute the sample mean of a fData object.
#'
#' It computes the \bold{cross-sectional} mean of a univariate functional
#' dataset, i.e., its time-by-time sample mean.
#'
#' @param fData the functional data object representing the dataset
#'
mean_fData = function( fData )
{
  return( fData( seq( fData$t0, fData$tP, length.out = fData$P ),
                colMeans( fData$values ) ) )

}

#'
#' \code{median_fData} method to compute the sample median of a fData object.
#'
#' It computes the depth-based median of a univariate functional dataset, i.e.
#' it finds the deepest element of the functional dataset according to a
#' specified definition of depth.
#'
#' @param fData the functional data object containing the dataset
#' @param type depth definition to use in order to find the sample median
#' (default is MBD)
#'
median_fData = function( fData, type = 'MBD' )
{
  Depths = eval( parse( text = paste( type, '( fData$values )', sep = '' ) ) )

  return( fData( seq( fData$t0, fData$tP, length.out = fData$P ),
                 fData$values[ which.max( Depths ), ] ) )
}

#'
#' Overload of [] access operator for univariate functional dataset.
#'
#' It provides an easy and natural access to subsets of a functional dataset
#' without having to deal with the inner representation of functional datasets
#' of fData
#'
#' @param fD the univariate functional dataset
#' @param i a valid expression to subset rows (observations) of the univariate
#' functional dataset
#' @param j a valid expression to subset columns (measurements) of the univariate
#' functional dataset
#' @param as_fData logical flag to specify if output should be returned as a fData
#' object containing the request subset, default is TRUE.
#'
"[.fData" = function( fD, i, j, as_fData = TRUE )
{
  if( as_fData == TRUE )
  {
    if( missing( j ) )
    {
      return( structure( list( t0 = fD$t0,
                               tP = fD$tP,
                               h = fD$h,
                               P = fD$P,
                               N = ifelse( missing( i ), fD$N, length( i ) ),
                               values = toRowMatrixForm( fD$values[ i, ] ) ),
                         class = c( 'fData' ) ) )
    } else {
      return( structure( list( t0 = fD$t0 + ( min( j ) - 1 ) * fD$h,
                               tP = fD$t0 + ( max( j ) - 1 ) * fD$h,
                               h = fD$h,
                               P = length( j ),
                               N = ifelse( missing( i ), fD$N, length( i ) ),
                               values = toRowMatrixForm( fD$values[ i, j ] ) ),
                         class = c( 'fData' ) ) )
    }
  } else {
    return( fD$values[ i, j ] )
  }
}


#'
#' Extracting list of components values from multivariate functional dataset
#'
#' It provides an easy way to extract from a multivariate functional dataset
#' with arbitrary dimensions a list of its time-by-time values stored each into
#' the corresponding matrix.
#'
#' @param mfData the multivariate functional dataset
#'
toListOfValues = function( mfData )
{
  eval( parse( text = paste( 'list(', paste( 'mfData$fDList[[ ',
                                             1 : mfData$L, ' ]]$values',
                         sep = '', collapse = ', ' ), ')' ) ) )
}


unfold = function( fData )
{
  return( fData( seq( fData$t0, fData$tP, length.out = fData$P ),
                 t( apply( fData$values,
                           1,
                           function( x )( c( x[ 1 ],
                                             x[ 1 ] +
                                               cumsum( abs( diff( x ) ) ) ) ) )
                 ) ) )
}


warp = function( fData, warpings )
{
  if( fData$t0 != warpings$t0 |
      fData$tP != warpings$tP |
      fData$P != warpings$P |
      fData$h != fData$h |
      fData$N != fData$N )
  {
    stop( ' Error in warp: you have to provide a warping for each functional
observations, and they need to be defined over the same grid' )
  }

  if( any( diff( warpings$values[ , 1 ] ) > .Machine$double.eps ) |
      any( diff( warpings$values[ , fData$P ] ) > .Machine$double.eps ) )
  {
    stop( ' Error in warp: you have to prescribe warpings with same starting
          point and ending point (at least up to .Machine$double.eps = ',
          .Machine$double.eps, ')' )
  }

  time_grid = seq( fData$t0,
                   fData$tP,
                   length.out = fData$P )

  return( fData( seq( warpings$t0,
                      warpings$tP,
                      length.out = warpings$P ),
                 t( sapply( 1 : fData$N,
                            function( i )(
                              approx( time_grid,
                                      fData$values[ i, ],
                                      xout = warpings$values[ i, ],
                                      yright = fData$values[ i, fData$P ],
                                      yleft = fData$values[ i, 1 ] )$y ) ) ) ) )
}


