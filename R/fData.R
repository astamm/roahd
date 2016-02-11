
#'
#' \code{fData} generic class for univariate functional data.
#'
#' @param grid the evenly spaced grid over which the functional data is defined
#' @param values the time-by-time values of the functional dataset, provided
#' either in a vector or matrix-like representation, with N rows (observations)
#' and P columns (time points)
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

#'
#' \code{plot.fData} specialised method to plot univariate functional data..
#'
#' @param fData the univariate functional data object
#' @param lty lty graphical parameter to be used in plotting functions
#' @param col colors to be used in plotting functions
#' @param ... additional graphical parameters to be used in plotting functions
#'
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

  if( unique( apply( dimMatrix, 1, function( x )( length( unique( x ) ) ) ) ) != 1 )
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

  return( structure( list( L = L,
                           N = fDList[[ 1 ]]$N,
                           P = fDList[[ 1 ]]$P,
                           fDList = fDList ),
                     class = c( 'mfData' ) ) )
}

#'
#' \code{plot.mfData} specialised method to plot multivariate functional data.
#'
#' @param mfData the multivariate functional data object
#' @param lty lty graphical parameter to be used in plotting functions
#' @param col colors to be used in plotting functions
#' @param ylab either a string or a list of stirngs to be used as y labels for
#' each window in the plot panel, default is null
#' @param main either a string or a list of stirngs to be used as main title for
#' each window in the plot panel, default is null
#' @param ... additional graphical parameters to be used in plotting functions
plot.mfData = function( mfData, lty = 1, col = NULL, ylab = NULL, main = NULL, ... )
{
  if( is.null( col ) )
  {
    col = set_alpha( scales::hue_pal( )( mfData$N ),
                     0.8 )
  }

  if( ! is.null( ylab ) )
  {
    if( length( ylab ) == 1 )
    {
      ylab = rep( ylab, mfData$L )
    } else if( length( ylab ) != mfData$L )
    {
      stop( 'Error in plot.mfData: you specified a wrong number of y labels' )
    }
  } else {
    ylab = rep( list( '' ), mfData$L )
  }

  if( ! is.null( main ) )
  {
    if( length( main ) == 1 )
    {
      main = rep( main, mfData$L )
    } else if( length( main ) != mfData$L )
    {
      stop( 'Error in plot.mfData: you specified a wrong number of subtitles' )
    }
  } else {
    main = rep( list( '' ), mfData$L )
  }

  mfrow_rows = ceiling( mfData$L / 2 )
  mfrow_cols = 2

  par( mfrow = c( mfrow_rows, mfrow_cols ) )

  plot_aux = function( i, ... )( plot( mfData$fDList[[ i ]], lty = lty, col = col,
                                  ylab = ylab[[ i ]],
                                  main = main[[ i ]],
                                  ... ))

  invisible( sapply( 1 : mfData$L, plot_aux, ... ) )
}


#'
#' Overload of + operator for fData objects.
#'
#'  It allows to perform sum between functional data objects or between a
#'  functional data and other compliant containers, like matrices or vectors,
#'  representing the pointwise measurements of the second term of the sum.
#'
#' @param fD the univariate functional data object
#' @param A either a functional data object, whose time_grid must be the same as
#' fD's, or a one-dimensional data structure (like 1D-array or raw numeric
#' vector), or a two-dimensional data structure (like 2D-array or raw numeric
#' matrix ), that specifies the second term of the sum.
#' In case of a one-dimensional data structure, the sum is performed element-wise
#' between each element of the functional dataset fD and A.
#' In case of a two-dimensional data structure, the sum is performet element-wise
#' between corresponding elements of fD and A's rows.
#'
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

#' Overload of - operator for fData objects.
#'
#'  It allows to perform subtractions between functional data objects or between
#'  a functional data and other compliant containers, like matrices or vectors,
#'  representing the pointwise measurements of the second term of the sum.
#'
#' @param fD the univariate functional data object
#' @param A either a functional data object, whose time_grid must be the same as
#' fD's, or a one-dimensional data structure (like 1D-array or raw numeric
#' vector), or a two-dimensional data structure (like 2D-array or raw numeric
#' matrix ), that specifies the second term of the subtraction.
#' In case of a one-dimensional data structure, the sum is performed element-wise
#' between each element of the functional dataset fD and A.
#' In case of a two-dimensional data structure, the sum is performet element-wise
#' between corresponding elements of fD and A's rows.
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
      fD$values = fD$values + A$values
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

#' Overload of * operator for fData objects.
#'
#'  It allows to perform multiplications between a functional data object and a
#'  either a numeric variable or numeric one-dimensional data structure.
#'
#' @param fD the univariate functional data object
#' @param a either a number or a numeric one-dimensional data strcture (array,
#' matrix or raw vector)
#'
"*.fData" = function( fD, a )
{
  if( ! 1 %in% dim( as.matrix( a ) ) )
  {
      stop( 'Error in *.fData: dimensions are not compliant' )
  }

  fD$values = fD$values * a

  return( fD )
}

#' Overload of / operator for fData objects.
#'
#'  It allows to perform divisions between a functional data object and a
#'  either a numeric variable or numeric one-dimensional data structure.
#'
#' @param fD the univariate functional data object
#' @param a either a number or a numeric one-dimensional data strcture ( array,
#' matrix or raw vector)
#'
"/.fData" = function( fD, a )
{
  if( ! 1 %in% dim( as.matrix( a ) ) )
  {
    stop( 'Error in /.fData: dimensions are not compliant' )
  }

  fD$values = fD$values / a

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
