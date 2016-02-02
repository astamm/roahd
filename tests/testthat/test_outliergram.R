

# TESTING OUTLIERGRAM -----------------------------------------------------


# OUTLIER DETECTION -------------------------------------------------------
require(mvtnorm)

# Synthetic data

set.seed( 10101 )

N = 100
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

quartz()
outliergram( time, S, display = TRUE, main = 'Example data', ylab = 'value', xlab = 'time' )

# l = lineprof( outliergram( time, S, display = FALSE ) )

# shine( l )

# Rprof()
# outliergram( time, S, display = FALSE )
#
# Rprof( NULL )
#
# outliergram( time, S )
