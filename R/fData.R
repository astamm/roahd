
#' \code{S3} class for univariate functional datasets.
#'
#'  This function implements a constructor for elements of \code{S3} class
#'  \code{fData}, aimed at implementing a representation of a functional
#'  dataset.
#'
#'  The functional dataset is represented as a collection of measurement of the
#'  observations on an evenly spaced, 1D grid of discrete points (representing,
#'  e.g. time), namely, for functional data defined over a grid \eqn{[t_0,
#'  t_1, \ldots, t_P]}:
#'
#'  \deqn{ f_{i,j} = f_i( t_0 + j h ), \quad h =  \frac{t_P - t_0}{N},
#'  \quad \forall j = 1, \ldots, P, \quad \forall i = 1, \ldots
#'  N.}
#'
#' @param grid the evenly spaced grid over which the functional observations are
#' measured. It must be a numeric vector of length \code{P}.
#' @param values the values of the observations in the functional dataset,
#' prodived in form of a 2D data structure (e.g. matrix or array) having as
#' rows the observations and as columns their measurements over the 1D grid of
#' length \code{P} specified in \code{grid}.
#'
#' @return The function returns a \code{S3} object of class \code{fData}, containing
#' the following elements:
#' \itemize{
#'  \item{"\code{N}"}{: the number of elements in the dataset;}
#'  \item{"\code{P}"}{: the number of points in the 1D grid over which elements
#'  are measured;}
#'  \item{"\code{t0}"}{: the starting point of the 1D grid;}
#'  \item{"\code{tP}"}{: the ending point of the 1D grid;}
#'  \item{"\code{values}"}{: the matrix of measurements of the functional
#'  observations on the 1D grid provided with \code{grid}.}
#' }
#'
#' @seealso \code{\link{generate_gauss_fdata}}
#'
#' @examples
#' # Defining parameters
#' N = 20
#' P = 1e2
#'
#' # One dimensional grid
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating an exponential covariance function (see related help for more
#' # information )
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Generating a synthetic dataset with a gaussian distribution and
#' # required mean and covariance function:
#' values = generate_gauss_fdata( N,
#'                                centerline = sin( 2 * pi * grid ),
#'                                Cov = C )
#'
#' fD = fData( grid, values )
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

#' Specialised method to plot \code{fData} objects
#'
#' This function performs the plot of a functional univariate dataset stored in
#' an object of class \code{fData}. It is able to accept all the usual
#' customisable graphical parameters, otherwise it will use the default ones.
#'
#' @param x the univariate functional dataset of \code{fData} class.
#' @param ... additional graphical parameters to be used in plotting functions
#'
#' @seealso \code{\link{fData}}
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' # One dimensional grid
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating an exponential covariance function (see related help for more
#' # information )
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Generating a synthetic dataset with a gaussian distribution and
#' # required mean and covariance function:
#' values = generate_gauss_fdata( N,
#'                                centerline = sin( 2 * pi * grid ),
#'                                Cov = C )
#'
#' fD = fData( grid, values )
#'
#' plot( fD )
#'
#' @export
#'
plot.fData = function( x, ... )
{
  .plot_fData( x, ...  )
}

.plot_fData = function( x,
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



#' \code{S3} class for multivariate functional datasets
#'
#' This function implements a constructor for elements of \code{S3} class
#' \code{mfData}, aimed at implementing a representation of a multivariate
#' functional dataset.
#'
#' The functional dataset is represented as a collection of \code{L} components,
#' each one an object of class \code{fData}. Each component must contain elements
#' defined on the same grid as the others, and must contain the same number of
#' elements (\code{N}).
#'
#' @param grid the (evenly spaced) grid over which the functional dataset is
#' defined.
#' @param Data_list a \code{list} containing the \code{L} components of the
#' multivariate functional dataset, defined as 2D data structures (e.g. matrix
#' or array) having as rows the \code{N} observations and as columns the
#' \code{P} measurements on the grid provided by \code{grid}.
#'
#' @return
#' The function returns a \code{S3} object of class \code{mfData}, containing
#' the following elements:
#' \itemize{
#'  \item{"\code{N}"}{: the number of elements in the dataset;}
#'  \item{"\code{L}"}{: the number of components of the functional dataset;}
#'  \item{"\code{P}"}{: the number of points in the 1D grid over which elements
#'  are measured;}
#'  \item{"\code{t0}"}{: the starting point of the 1D grid;}
#'  \item{"\code{tP}"}{: the ending point of the 1D grid;}
#'  \item{"\code{fDList}"}{: the list of \code{fData} objects representing the
#'  \code{L} components as corresponding unviariate functional datasets.}
#' }
#'
#' @seealso \code{\link{fData}}, \code{\link{generate_gauss_fdata}},
#' \code{\link{generate_gauss_mfdata}}
#'
#' @examples
#' # Defining parameters
#' N = 1e2
#'
#' P = 1e3
#'
#' t0 = 0
#' t1 = 1
#'
#' # Defining the measurement grid
#' grid = seq( t0, t1, length.out = P )
#'
#' # Generating an exponential covariance matrix to be used in the simulation of
#' # the functional datasets (see the related help for details)
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating the measurements of two univariate functional datasets with
#' # required center and covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#'
#' # Building the mfData object
#' mfData( grid, list( Data_1, Data_2 ) )
#'
#' @export
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
                           L = L,
                           P = fDList[[ 1 ]]$P,
                           t0 = fDList[[ 1 ]]$t0,
                           tP = fDList[[ 1 ]]$tP,
                           fDList = fDList ),
                     class = c( 'mfData' ) ) )
}

#' Specialised method to plot \code{mfData} objects
#'
#' This function performs the plot of a functional multivariate dataset stored
#' in an object of class \code{mfData}. It is able to accept all the usual
#' customisable graphical parameters, otherwise it will use the default ones.
#'
#' The current active graphical device is split into a number of sub-figures,
#' each one meant to contain the plot of the corresponding dimension of the
#' \code{mfData} object. In particular, they are arranged in a rectangular
#' lattice with a number of rows equal to \eqn{ \lfloor \sqrt{ L } \rfloor }
#' and a number of columns equal to \eqn{ \lceil L / \lfloor \sqrt{L} \rfloor
#' \rceil }.
#'
#' A special use of the graphical parameters allows to set up y-labels and
#' titles for all the sub-figures in the graphical window. In particular,
#' parameters \code{ylab} and \code{main} can take as argument either a single
#' string, that are repeatedly used for all the sub-graphics, or a list of
#' different strings (one for each of the \code{L} dimensions) that have to be
#' used in the corresponding graphic.
#'
#' @param x the multivariate functional dataset of \code{mfData} class.
#' @param ... additional graphical parameters to be used in plotting functions
#' (see \code{Details} for the use of \code{ylab} and \code{main}).
#'
#' @seealso \code{\link{mfData}}, \code{\link{fData}}, \code{\link{plot.fData}}
#'
#' @examples
# Defining parameters
#' N = 1e2
#'
#' P = 1e3
#'
#' t0 = 0
#' t1 = 1
#'
#' # Defining the measurement grid
#' grid = seq( t0, t1, length.out = P )
#'
#' # Generating an exponential covariance matrix to be used in the simulation of
#' # the functional datasets (see the related help for details)
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating the measurements of two univariate functional datasets with
#' # required center and covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#'
#' # Building the mfData object and plotting tt
#' plot( mfData( grid, list( Data_1, Data_2 ) ),
#'       xlab = 'time', ylab = list( '1st dim.', '2nd dim.' ),
#'       main = list( 'An important plot here', 'And another one here' ) )
#'
#' @export
#'
plot.mfData = function( x, ... )
{
  mfrow_rows = floor( sqrt( x$L ) )
  mfrow_cols = ceiling( x$L / floor( sqrt( x$L ) ) )

  par( mfrow = c( mfrow_rows, mfrow_cols ) )

  .plot_mfData( x, ... )

}

.plot_mfData = function( x,
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


