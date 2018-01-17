BCInterval = function( fD1, fD2, ordering='MEI', bootstrap_iterations=1000, alpha=0.05,
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
    int1 = v[ floor( ( bootstrap_iterations + 1 ) * alpha1 ) ]

  # it's the case in which alpha2 * bootstrap_iterations is an integer
  if( alpha2 * bootstrap_iterations - floor( alpha2 * bootstrap_iterations ) == 0 )
    int2 = v[ alpha2 * bootstrap_iterations ]
  else
    # when alpha2*bootstrap_iterations is not an integer,
    # a modification is used
    int2 = v[ floor( ( bootstrap_iterations + 1 ) * alpha2 ) ]

  return( list( lower=int1, upper=int2 ) )
}


BCIntervalMultivariate = function(mfD,
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
      interval = BCInterval( mfD$fDList[[ iL ]],
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


bootstrapTest = function( mfD1, mfD2, bootstrap_iterations=1000, ordering='MEI', normtype='f',
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
  diff = cor_spearman( mfD1, ordering=ordering) - cor_spearman( mfD2, ordering=ordering )
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
    diff_0 =  cor_H0_1 - cor_H0_2

    v[iBoot] = norm( diff_0 , type=normtype) # bootstrap replication of the test statistic under H0
  }

  pvalue =  sum( v >= phi ) / bootstrap_iterations

  return(list( pvalue=pvalue, phi=phi ))
}


append_fData = function(fD1, fD2)
{
  stopifnot((fD1$P == fD2$P) & (fD1$t0 == fD2$t0) & (fD1$tP == fD2$tP)
            & (fD1$h == fD2$h))

  grid = seq(fD1$t0, fD1$tP, length.out=fD1$P)

  return(fData(grid, values = rbind(fD1$values, fD2$values)))
}

append_mfData = function(mfD1, mfD2)
{
  stopifnot((mfD1$P == mfD2$P) & (mfD1$t0 == mfD2$t0) & (mfD1$tP == mfD2$tP)
            & (mfD1$fDList[[1]]$h == mfD2$fDList[[1]]$h) & (mfD1$L == mfD2$L))

  return(as.mfData(lapply(1:mfD1$L,
                          function(id) append_fData(mfD1$fDList[[id]], mfD2$fDList[[id]]))))
}


