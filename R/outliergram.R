

outliergram = function( time_grid, Data, MBD_data = NULL, MEI_data = NULL, q_low = 0, q_high = 1, display = TRUE, ... )
{

  require( scales )

#   #####
#
#   N = 100
#   P = 200
#   time = seq( 0, 1, length.out = P )
#
#   S = matrix( sin( 4 * pi * time), ncol = P, nrow = N, byrow = T )
#
#   eps = array( rnorm( N * P ), dim = c( N,  P) )
#
#   # Exponential covariance function
#   C = outer( time, time, function( s, t ){ 0.2 * exp( - 0.8 * abs( s - t ) ) } )
#
#   # Spectral decomposition
#   C.eigen = eigen(C)
#
#   C.sqrt = C.eigen$vectors %*% sqrt( diag( C.eigen$values ) ) %*% t( C.eigen$vectors )
#
#   # Affine transformation of i.i.d gaussian data
#   eps = eps %*% C.sqrt
#
#   S = S + eps;
#
#   # Number of outliers
#   N_outliers = 4
#   eps.outliers = array( rnorm( N_outliers * P  ), dim = c( N_outliers, P ) )
#
#   eps.outliers = eps.outliers %*% C.sqrt
#
#   S.outliers = array( 0, dim = c( N_outliers, P ) )
#
#   S.outliers[ 1, ] = sin( 4 * pi * time ) + 2 + eps.outliers[ 1, ]
#   S.outliers[ 2, ] = sin( 4 * pi * time ) - 2 + eps.outliers[ 2, ]
#   S.outliers[ 3, ] =  eps.outliers[3,]
#   S.outliers[ 4, ] = sin( 4 * pi * time + pi/3 ) + eps.outliers[ 4, ]
#
#   S = rbind( S, S.outliers)
#   Data = S
#
#   q_low = 0
#
#   q_high = 1
#
#   ####

  N = nrow( Data )

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Computing MBD
  if( is.null( MBD_data ) ){

    MBD_data = MBD( Data )
  }

  # Computing MEI
  if( is.null( MEI_data ) )
  {
    MEI_data = MEI( Data )
  }

  P = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2

  Q_d = quantile( P - MBD_data  )

  Q_d3 = Q_d[ 4 ]

  IQR = Q_d[ 4 ] - Q_d[ 2 ]

  # Computing surely outlying curves
  ID_shape_outlier = which( MBD_data <= P - Q_d3 - 1.5 * IQR )


  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : nrow( Data ), ID_shape_outlier )

  # Low MEI curves will be checked for upward shift
  ID_non_outlying_Low_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] < 0.5 ) ]

  # High MEI curves will be checked for downward shift
  ID_non_outlying_High_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] >= 0.5 ) ]

  # Managing high MEI data
  for( iObs in ID_non_outlying_High_MEI )
  {
    diff_min_vector = Data[ iObs, ] - apply( Data[ - iObs, ], 2, ifelse( q_low == 0,
                                                                         min,
                                                                         function( x )( quantile( x, probs = q_low ) ) ) )

    min_diff_min = min( diff_min_vector )

    Data_tilde = Data[ iObs, ] - rep( any ( diff_min_vector < 0 ) * min_diff_min, ncol( Data ) )

    MBD_curr = MBD( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    MEI_curr = MEI( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    P_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2

    if( MBD_curr <= P_curr - Q_d3 - 1.5 * IQR )
    {
      ID_shape_outlier = c( ID_shape_outlier, iObs )
    }
  }

  # Managing low MEI data
  for( iObs in ID_non_outlying_Low_MEI )
  {
    diff_max_vector = Data[ iObs, ] - apply( Data[ - iObs, ], 2, ifelse( q_high == 1,
                                                                         max,
                                                                         function( x )( quantile( x, probs = q_high ) ) ) )

    max_diff_max = max( diff_max_vector )

    Data_tilde = Data[ iObs, ] - rep( any ( diff_max_vector > 0 ) * max_diff_max, ncol( Data ) )

    MBD_curr = MBD( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    MEI_curr = MEI( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    P_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2

    if( MBD_curr <= P_curr - Q_d3 - 1.5 * IQR )
    {
      ID_shape_outlier = c( ID_shape_outlier, iObs )
    }
  }

#   require(foreach)
#
#   tic = proc.time()
#   foreach( idObs = ID_non_outlying_High_MEI, .combine = 'c' ) %do%
#     check_outlier_High_MEI( Data, idObs, q_low )
#
#   foreach( idObs = ID_non_outlying_Low_MEI, .combine = 'c' ) %do%
#     check_outlier_Low_MEI( Data, idObs, q_high )
#
#   toc = proc.time()
#   cat( toc - tic )


  # Plot

  if( display )
  {
    # Setting up palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ), l = 60 )( length( ID_non_outlying ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ), c = 150 )( length( ID_shape_outlier ) )
    # col_labels_outlying = scales::hue_pal( h = c( -10, 10 ), c = 1000, l = 20 )( length( ID_shape_outlier ) )

    dev.cur()
    par( mfrow = c( 1, 2 ) )

    # Plotting functional data
    matplot( time_grid, t( Data[ - ID_shape_outlier, ] ), type = 'l', lty = 1, ylim = range( Data ),
             col = col_non_outlying, ... )
    matplot( time_grid, t( Data[   ID_shape_outlier, ] ), type = 'l', lty = 1, lwd = 2, ylim = range( Data ),
             col = col_outlying, add = TRUE )

    # Adding text labels with curve ID
    w_spacing = diff( range( time_grid ) ) / ( 2 * length( ID_shape_outlier ) )

    for( iOut in seq_along( ID_shape_outlier ) )
    {
      text( time_grid[ 1 ] + ( 2 * iOut - 1 ) * w_spacing,
            Data[ ID_shape_outlier[ iOut ],
                  which.min( abs( time_grid - time_grid[ 1 ] - ( 2 * iOut - 1 ) * w_spacing ) ) ] +
              diff( range( Data[ ID_shape_outlier[ iOut ]  ] ) ) / 30,
            ID_shape_outlier[ iOut ],
            col = col_outlying[ iOut ] )
    }

    # Plotting outliergram
    plot( sort( MEI_data ), a_0_2 + a_1 * sort( MEI_data ) + a_0_2 * N^2 * sort( MEI_data )^2,
          lty = 2, type = 'l', col = 'darkblue', lwd = 2, ylim = c( 0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 ),
          main = 'Outliergram', xlab = 'MEI', ylab = 'MBD' )
    points( MEI_data[ - ID_shape_outlier ], MBD_data[ - ID_shape_outlier ],
            pch = 16, col = col_non_outlying )
    points( MEI_data[ ID_shape_outlier ], MBD_data[ ID_shape_outlier ], pch = 16, col = col_outlying )

    for( idOut in ID_shape_outlier )
    {
      text( MEI_data[ idOut ],
            MBD_data[ idOut ] + 0.5 / 30,
            idOut,
            col = col_outlying[ match( idOut, ID_shape_outlier ) ] )
    }
  }

  return( ID_shape_outlier )
}


# check_outlier_Low_MEI = function( Data, idObs, q_high = 1  )
# {
#   diff_max_vector = Data[ idObs, ] - apply( Data[ - idObs, ],
#                                             2,
#                                             ifelse( q_high == 1,
#                                                     max,
#                                                     function( x )( quantile( x, probs = q_high ) ) ) )
#
#   max_diff_max = max( diff_max_vector )
#
#   Data_tilde = Data[ idObs, ] - any ( diff_max_vector > 0 ) * max_diff_max
#
#   MBD_curr = MBD( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   MEI_curr = MEI( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   P_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2
#
#   return( as.logical( MBD_curr <= P_curr - Q_d3 - 1.5 * IQR ) )
# }
#
# check_outlier_High_MEI = function( Data, idObs, q_low = 0 )
# {
#   diff_min_vector = Data[ idObs, ] - apply( Data[ - idObs, ],
#                                             2,
#                                             ifelse( q_low == 0,
#                                                     min,
#                                                     function( x ) ( quantile( x, probs = q_low ) ) ) )
#
#   min_diff_min = min( diff_min_vector )
#
#   Data_tilde = Data[ idObs, ] - any ( diff_min_vector < 0 ) * min_diff_min
#
#   MBD_curr = MBD( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   MEI_curr = MEI( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   P_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2
#
#   return( as.logical( MBD_curr <= P_curr - Q_d3 - 1.5 * IQR ) )
# }

set_alpha = function( col, alpha )
{
  alpha = alpha * 255

  rgb_colors = rbind( col2rgb( col ), alpha = alpha )

  return( apply( rgb_colors, 2, function( x )( rgb( x[ 1 ],
                                                    x[ 2 ],
                                                    x[ 3 ],
                                                    x[ 4 ],
                                                    maxColorValue = 255 ) ) ) )
}

