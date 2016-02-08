

# TESTING OUTLIERGRAM -----------------------------------------------------

require(mvtnorm)

# Synthetic data

set.seed( 1618 )

N = 200
P = 200
time = seq( 0, 1, length.out = P )

S = matrix( sin( 4 * pi * time), ncol = P, nrow = N, byrow = T )

eps = array( rnorm( N * P ), dim = c( N,  P) )

# Exponential covariance function
C = outer( time, time, function( s, t ){ 0.2 * exp( - 0.8 * abs( s - t ) ) } )

# Spectral decomposition
C.eigen = eigen(C)

C.sqrt = C.eigen$vectors %*% sqrt( diag( C.eigen$values ) ) %*% t( C.eigen$vectors )

# Affine transformation of i.i.d gaussian data
eps = eps %*% C.sqrt

S = S + eps;

# Number of outliers
N_outliers = 4
eps.outliers = array( rnorm( N_outliers * P  ), dim = c( N_outliers, P ) )

eps.outliers = eps.outliers %*% C.sqrt

S.outliers = array( 0, dim = c( N_outliers, P ) )

S.outliers[ 1, ] = sin( 4 * pi * time ) + 2 + eps.outliers[ 1, ]
S.outliers[ 2, ] = sin( 4 * pi * time ) - 2 + eps.outliers[ 2, ]
S.outliers[ 3, ] =  eps.outliers[3,]
S.outliers[ 4, ] = sin( 4 * pi * time + pi/3 ) + eps.outliers[ 4, ]

S = rbind( S, S.outliers)

mbd = MBD( S )
mei = MEI( S )

outliergram( time, S )

outliergram( time, S, adjust = list( N_trials = 5,
                                     TPR = 0.005,
                                     VERBOSE = TRUE ) )


#### CORRECTNESS BENCHMARKING

# out_dot = dot_outliergram( time, S )
#
# out_my = my_outliergram( time, S )
#
# par_my = par_outliergram( time, S )

# l = lineprof( outliergram( time, S, display = FALSE ) )
# Rprof( filename =  'prof.out', memory = 'both' )
# invisible( my_outliergram( time, S ) )
# Rprof( NULL )
#
# summ = summaryRprof( 'prof.out' )
# print( summ )
#
# plot( parse_rprof( 'prof.out' ) )
#
# pippo = readProfileData("prof.out")


#### EFFICIENCY BENCHMARKING
# library(doParallel)
# library(foreach)
#
# tic = proc.time()
# invisible( dot_outliergram( time, S ) )
# toc = proc.time()
#
# time_old = toc - tic
#
# tic = proc.time()
# invisible( my_outliergram( time, S, p_check = 0.05 ) )
# toc = proc.time()
#
# time_1 = toc - tic
# print( time_1 )
#
# tic = proc.time()
# invisible( par_outliergram( time, S ) )
# toc = proc.time()
#
# time_2 = toc - tic
#
# time_old
# time_1
# time_2
#
# time_old / time_1
# time_old / time_2
# time_1 / time_2

# outliergram( time, S, display = TRUE, main = 'Example data', ylab = 'value', xlab = 'time' )
# outliergram( time, S, adjust = list( N_trials = 2,
#                                      trial_size = 500,
#                                      VERBOSE = 1,
#                                      tol = 1e-1 ),
#              display = TRUE, main = 'Example data', ylab = 'value', xlab = 'time' )

# l = lineprof( outliergram( time, S, display = FALSE ) )

# shine( l )

# Rprof()
# outliergram( time, S, display = FALSE )
#
# Rprof( NULL )
#
# outliergram( time, S )
