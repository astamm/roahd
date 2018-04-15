#' Bootstrap Confidence Interval on Spearman's Correlation Coefficient between Univariate Functional Datasets
#'
#' This function computes the bootstrap confidence interval of coverage probability
#' \eqn{1 - \alpha} for the Spearman correlation coefficient between two univariate functional samples.
#'
#' The function takes two samples of compatible functional data (i.e., they must be defined over the
#' same grid and have same number of observations) and computes a bootstrap confidence interval for
#' their Spearman correlation coefficient.
#'
#' @param fD1 is the first univariate functional sample in form of an \code{fData} object.
#' @param fD2 is the first univariate functional sample in form of an \code{fData} object.
#' @param ordering is either \code{MEI} (default) or \code{MHI}, and indicates the ordering relation
#' to be used during in the Spearman's coefficient computation.
#' @param bootstrap_iterations is the number of bootstrap iterations to use in order to estimate the
#' confidence interval (default is 1000).
#' @param alpha controls the coverage probability (1-\code{alpha}).
#' @param verbose whether to log information on the progression of bootstrap iterations.
#'
#' @return The function returns a list of two elements, \code{lower} and \code{upper}, representing
#' the lower and upper end of the bootstrap confidence interval.
#'
#' @seealso \code{\link{cor_spearman}}, \code{\link{cor_spearman_accuracy}}, \code{\link{fData}},
#' \code{\link{mfData}}, \code{\link{BCIntervalSpearmanMultivariate}}
#'
#' @examples
#'
#' set.seed(1)
#'
#' N = 2e2
#' P = 1e2
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
#' BCIntervalSpearman(mfD$fDList[[1]], mfD$fDList[[2]], ordering = 'MEI')
#'
#' BCIntervalSpearman(mfD$fDList[[1]], mfD$fDList[[2]], ordering = 'MHI')
#'
#' # BC intervals contain zero since the functional samples are uncorrelated.
#'
#' @export
#'
BCIntervalSpearman = function( fD1, fD2, ordering='MEI', bootstrap_iterations=1000, alpha=0.05,
                       verbose=FALSE ){

  if( verbose ){
    printout_iters = ceiling(seq(1, bootstrap_iterations, length.out=11))
  }

  stopifnot((fD1$N == fD2$N) & (fD1$P == fD2$P) & (fD1$t0 == fD2$t0) & (fD1$tP == fD2$tP)
            & (fD1$h == fD2$h))

  N = fD1$N
  v = matrix(0, nrow=bootstrap_iterations)

  for ( iBoot in 1:bootstrap_iterations )
  {
    if ((verbose) && (any(iBoot == printout_iters))){
      message(paste0(round(iBoot/bootstrap_iterations, 2) * 100,
                     '% bootstrap Spearman`s correlation'))
    }

    w = sample.int( N, size=N, replace=TRUE, prob=NULL )
    v[iBoot] = cor_spearman( as.mfData(list(fD1[w,], fD2[w,])), ordering=ordering )
  }

  # bias-correction
  pz0 = sum( v < cor_spearman( as.mfData(list(fD1, fD2)), ordering=ordering ) ) / bootstrap_iterations
  z0 = qnorm( pz0, mean=0, sd=1 )

  # vector of jackknife cor_spearman values
  theta_i = matrix(0, nrow=1, ncol=N)

  for ( i in 1:N )
  {
    theta_i[i] = cor_spearman( as.mfData(list(fD1[-i,], fD2[-i,])), ordering=ordering )
  }
  theta_hat = mean( theta_i )

  # acceleration
  a = sum( ( theta_hat - theta_i )^3 )/( 6 * sum( ( theta_hat - theta_i )^2 )^( 1.5 ) )

  # first percentile
  alpha1 = pnorm( z0 +  (z0 + qnorm( alpha/2, mean=0, sd = 1 ) )/
                    ( 1 - a*( z0 + qnorm( alpha/2, mean = 0, sd = 1 ) ) ), mean = 0, sd = 1)
  # second percentile
  alpha2 = pnorm( z0 + (z0 + qnorm( 1-alpha/2, mean = 0,sd = 1 ) )/
                    ( 1 - a*( z0 + qnorm( 1 - alpha/2, mean = 0, sd = 1 ) ) ), mean = 0, sd = 1)
  v = sort( v )

  # it's the case in which alpha1 * bootstrap_iterations is an integer
  if( alpha1 * bootstrap_iterations - floor( alpha1 * bootstrap_iterations ) == 0 )
    int1 = v[ alpha1 * bootstrap_iterations ]
  else
    # when alpha1 * bootstrap_iterations is not an integer,
    # a modification is used
    id = min(c(floor( ( bootstrap_iterations + 1 ) * alpha1 ), bootstrap_iterations))
    id = max(c(id, 1))
    int1 = v[ id ]

  # it's the case in which alpha2 * bootstrap_iterations is an integer
  if( alpha2 * bootstrap_iterations - floor( alpha2 * bootstrap_iterations ) == 0 )
    int2 = v[ alpha2 * bootstrap_iterations ]
  else
    # when alpha2*bootstrap_iterations is not an integer,
    # a modification is used
    id = min( c( floor( ( bootstrap_iterations + 1 ) * alpha2 ), bootstrap_iterations ) )
    id = max( c( id, 1 ) )
    int2 = v[ floor( ( bootstrap_iterations + 1 ) * alpha2 ) ]

  return( list( lower=int1, upper=int2 ) )
}

#' Bootstrap Confidence Interval on Spearman's Correlation Coefficient of a Multivariate Functional Dataset
#'
#' This function computes the bootstrap confidence intervals of coverage probability
#' \eqn{1 - \alpha} for the Spearman correlation coefficients within a multivariate functional dataset.
#'
#' The function takes a multivariate functional dataset and computes a matrix of bootstrap confidence
#' intervals for its Spearman correlation coefficients.
#'
#' @param mfD is the multivariate functional sample in form of \code{mfData} object.
#' @param ordering is either \code{MEI} (default) or \code{MHI}, and indicates the ordering relation
#' to be used during in the Spearman's coefficient computation.
#' @param bootstrap_iterations is the number of bootstrap iterations to use in order to estimate the
#' confidence intervals (default is 1000).
#' @param alpha controls the coverage probability (1-\code{alpha}).
#' @param verbose whether to log information on the progression of bootstrap iterations.
#'
#' @return The function returns a list of two elements, \code{lower} and \code{upper}, representing
#' the matrices of lower and upper ends of the bootstrap confidence intervals for each pair of
#' components. The elements on the main diagonal are set to 1.
#'
#' @seealso \code{\link{cor_spearman}}, \code{\link{cor_spearman_accuracy}}, \code{\link{fData}},
#' \code{\link{mfData}}, \code{\link{BCIntervalSpearman}}
#'
#' @examples
#'
#' set.seed(1)
#'
#' N = 2e2
#' P = 1e2
#' grid = seq( 0, 1, length.out = P )
#'
#' # Creating an exponential covariance function to simulate guassian data
#' Cov = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating (independent) gaussian functional data with given center and
#' # covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 4 * pi * grid ), Cov = Cov )
#' Data_3 = generate_gauss_fdata( N, centerline = sin( 6 * pi * grid ), Cov = Cov )
#'
#' # Using the simulated data as (independent) components of a multivariate functional
#' # dataset
#' mfD = mfData( grid, list( Data_1, Data_2, Data_3 ) )
#'
#' BCIntervalSpearmanMultivariate(mfD, ordering = 'MEI')
#'
#' # BC intervals contain zero since the functional samples are uncorrelated.
#'
#' @export
#'
BCIntervalSpearmanMultivariate = function(mfD,
                                          ordering='MEI',
                                          bootstrap_iterations=1000,
                                          alpha=0.05,
                                          verbose=FALSE)
{

  lower = matrix(0, nrow=mfD$L, ncol=mfD$L)
  upper = matrix(0, nrow=mfD$L, ncol=mfD$L)

  diag(lower) = 1
  diag(upper) = 1

  for ( iL in 1:( mfD$L - 1 ) )
  {
    for( jL in (iL + 1): mfD$L )
    {
      if( verbose )
      {
        message(paste0('Bootstrap of components (', iL, ', ', jL, ')'))
      }

      interval = BCIntervalSpearman( mfD$fDList[[ iL ]],
                                     mfD$fDList[[ jL ]],
                                     ordering=ordering,
                                     bootstrap_iterations=bootstrap_iterations,
                                     alpha=alpha,
                                     verbose=verbose )

      lower[ iL, jL ] = interval$lower
      lower[ jL, iL ] = lower[iL, jL]

      upper[ iL, jL ] = interval$upper
      upper[ jL, iL ] = upper[ iL, jL ]
    }
  }

  return(list( lower=lower, upper=upper))
}

#' Bootstrap Hypothesis Test on Spearman Correlation Coefficients for Multivariate Functional Data
#'
#' This function performs a bootstrap test that checks whether the Spearman correlation structures
#' (e.g. matrices) of two populations of compatible multivariate functional data are equal or not.
#'
#' Given a first multivariate functional population, \eqn{X_1^(i), \ldots, X_n^(i)} with \eqn{i=1, \ldots, L},
#' defined on the grid \eqn{I = t_1, \ldots, t_P}, and a second multivariate functional population,
#' \eqn{Y_1^(i), \ldots, Y_m^(i)} with \eqn{i=1, \ldots, L}, defined on the same grid \eqn{I}, the
#' function performs a bootstrap test to check the hypothesis:
#'
#' \deqn{H_0: R_X = R_Y}
#' \deqn{H_1: R_X \neq R_Y,}
#'
#' where R_X and R_Y denote the L x L matrices of Spearman correlation coefficients of the two populations.
#'
#' \itemize{
#' \item{The two functional samples must have the same number of components and must be defined over the same
#' discrete interval \eqn{t_1, \ldots, t_P}.}
#' \item{The test is performed through a bootstrap argument, so
#' a number of bootstrap iterations must be specified as well. A high value for this parameter may result
#' in slow performances of the test (you may consider setting \code{verbose} to \code{TRUE} to get
#' hints on the process).}}
#'
#' @param mfD1 is the first functional dataset, specified in form of \code{mfData} object; it must
#' be compatible with \code{mfD2}.
#' @param mfD2 is the second functional dataset, specified in form of \code{mfData} object; it must
#' be compatible with \code{mfD1}.
#' @param bootstrap_iterations is the number of bootstrap iterations to be performed.
#' @param ordering is the kind of ordering to be used in the computation of Spearman's correlation
#' coefficeint (default is \code{MEI}).
#' @param normtype is the norm to be used when comparing the Spearman correlation matrices of the two
#' functional datasets (default is Frobenius, allowed values are the same as for parameter \code{type} in
#' the base function \code{norm}).
#' @param verbose a boolean flag specifying whether to print the progress of bootstrap iterations or not (default is FALSE).
#'
#' @return
#' The function returns the estimates of the test's p-value and statistics.
#'
#' @seealso \code{\link{BCIntervalSpearman}}, \code{\link{BCIntervalSpearmanMultivariate}}, \code{\link{mfData}}
#' @examples
#'
#' set.seed(1)
#' N = 2e2
#' P = 1e2
#' L = 2
#' grid = seq( 0, 1, length.out = P )
#' # Creating an exponential covariance function to simulate guassian data
#' Cov = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#'
#'
#' # Simulating two populations of bivariate functional data
#' #
#' # The first population has very high correlation between first and second component
#' centerline_1 = matrix(rep(sin(2 * pi * grid)), nrow = 2, ncol=P, byrow=TRUE)
#' values1 = generate_gauss_mfdata( N, L, correlations = 0.9,
#'                                  centerline = centerline_1, listCov = list(Cov, Cov)  )
#' mfD1 = mfData(grid, values1)
#'
#' # Pointwise estimate
#' cor_spearman(mfD1)
#'
#' # The second population has zero correlation between first and second component
#' centerline_2 = matrix(rep(cos(2 * pi * grid)), nrow = 2, ncol=P, byrow=TRUE)
#' values2 = generate_gauss_mfdata( N, L, correlations = 0,
#'                                  centerline = centerline_1, listCov = list(Cov, Cov)  )
#' mfD2 = mfData(grid, values2)
#'
#' # Pointwise estimate
#' cor_spearman(mfD2)
#'
#' # Applying the test
#' BTestSpearman(mfD1, mfD2)
#'
#' @export
#'
BTestSpearman = function( mfD1, mfD2, bootstrap_iterations=1000, ordering='MEI', normtype='f',
                          verbose=FALSE)
{
  stopifnot((mfD1$P == mfD2$P) & (mfD1$t0 == mfD2$t0) & (mfD1$tP == mfD2$tP)
            & (mfD1$fDList[[1]]$h == mfD2$fDList[[1]]$h) & (mfD1$L == mfD2$L))

  N1 = mfD1$N
  N2 = mfD2$N

  if( verbose ){
    printout_iters = ceiling(seq(1, bootstrap_iterations, length.out=11))
  }

  # Test statistic
  diff = as.matrix(cor_spearman( mfD1, ordering=ordering) - cor_spearman( mfD2, ordering=ordering ))
  phi = norm( diff, type=normtype )

  # Bootstrapping
  v = rbind( rep( 0, bootstrap_iterations ) )

  for( iBoot in 1 : bootstrap_iterations )
  {
    if ((verbose) && (any(iBoot == printout_iters))){
      message(paste0(round(iBoot/bootstrap_iterations, 2) * 100,
                     '% bootstrap Spearman`s correlation'))
    }

    idx1 = sample.int( N1 + N2, size=N1, replace=TRUE, prob=NULL )
    idx2 = sample.int( N1 + N2, size=N2, replace=TRUE, prob=NULL )

    # Building two multivariate functional datasets under H0 through resampling from mfD1 and mfD2

    idx_1_1 = idx1[idx1 <= N1]
    idx_1_2 = idx1[idx1 > N1]  - N1
    idx_2_1 = idx2[idx2 <= N1]
    idx_2_2 = idx2[idx2 > N1] - N1

    cor_H0_1 = cor_spearman( append_mfData(mfD1[idx_1_1], mfD2[idx_1_2]), ordering=ordering )
    cor_H0_2 = cor_spearman( append_mfData(mfD1[idx_2_1], mfD2[idx_2_2]), ordering=ordering )
    diff_0 =  as.matrix(cor_H0_1 - cor_H0_2)

    v[iBoot] = norm( diff_0 , type=normtype) # bootstrap replication of the test statistic under H0
  }

  pvalue =  sum( v >= phi ) / bootstrap_iterations

  return(list( pvalue=pvalue, phi=phi ))
}
