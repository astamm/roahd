

# TESTING FDATA -----------------------------------------------------------

N = 1e2

P = 1e3
t0 = 0
t1 = 1

time_grid = seq( t0, t1, length.out = P )

# C( s, t ) = \alpha \exp( - beta | s - t | )
# Amplitude factor
alpha = 0.3

# Correlation decay factor
beta = 0.4

Cov = outer( time_grid, time_grid, function( s, t )( alpha * exp( - beta * abs( s - t ) ) ) )

Data = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )

fD = fData( time_grid, Data )

quartz()
plot( fD, xlab = 'time', ylab = 'values', main = 'A functional dataset' )


# TESTING MFDATA ----------------------------------------------------------

N = 1e2

P = 1e3

t0 = 0
t1 = 1

time_grid = seq( t0, t1, length.out = P )

# C( s, t ) = \alpha \exp( - beta | s - t | )
# Amplitude factor
alpha = 0.3

# Correlation decay factor
beta = 0.4

Cov = outer( time_grid, time_grid, function( s, t )( alpha * exp( - beta * abs( s - t ) ) ) )

Data_1 = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )
Data_2 = generate_gauss_fdata( N, center = sin( 2 * pi * time_grid ), Cov = Cov )

mfD = mfData( time_grid, list( Data_1, Data_2 ) )

quartz()
plot( mfD, xlab = 'time', ylab = 'value'  )
