

outliergram = function( time_grid = NULL, Data, MBD_data = NULL, MEI_data = NULL, q_low = 0, q_high = 1,
                        adjust = FALSE, display = TRUE, ... )
{
  require( scales )
  require( robustbase )

  N = nrow( Data )

  if( is.null( time_grid ) )
  {
    time_grid = 1 : ncol( Data)
  }

  stopifnot( length( time_grid ) == ncol( Data ) )

  if( ! is.list( adjust ) )
  {
    # Plain outliergram with default F value: F = 1.5

    out = .outliergram( time_grid, Data, MBD_data, MEI_data,
                        q_low, q_high)

  } else {

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         5 * N,
                         adjust$trial_size )

    TPR = ifelse( is.null( adjust$TPR ),
                  2 * pnorm( 4 * qnorm( 0.25 ) ),
                  adjust$FPR )

    F_min = ifelse( is.null( adjust$F_min ),
                    0.5,
                    adjust$F_min )

    F_max= ifelse( is.null( adjust$F_max ),
                    5,
                    adjust$F_max )

    tol = ifelse( is.null( adjust$tol ),
                  1e-3,
                  adjust$tol )

    maxiter = ifelse( is.null( adjust$maxiter ),
                      100,
                      adjust$maxiter )

    VERBOSE = ifelse( is.null( adjust$VERBOSE ),
                      FALSE,
                      adjust$VERBOSE )

    Cov = covOGK( Data, sigmamu = s_Qn )$cov

    CholCov <- chol( Cov )

    if( is.null ( MBD_data ) )
    {
      MBD_data = MBD( Data )
    }

    centerline = Data[ which.max( MBD_data ), ]

    Fvalues = rep( 0, N_trials )


    for( iTrial in 1 : N_trials )
    {
      if( VERBOSE > 0 )
      {
        cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
      }

      Data_gauss = generate_gauss_fdata( trial_size, centerline, CholCov = CholCov )

      cat( ' * * * * beginning optimisation\n' )
      opt = uniroot( function( F_curr )( length( .outliergram( time_grid,
                                                                    Data_gauss,
                                                                    MBD_data = NULL,
                                                                    MEI_data = NULL,
                                                                    q_low = q_low,
                                                                    q_high = q_high,
                                                                    Fvalue = F_curr )$ID_SO ) / trial_size
                                              - TPR ),
                     interval = c( F_min,
                                   F_max ),
                     tol = tol,
                     maxiter = maxiter )

      Fvalues[ iTrial ] = opt$root
    }

    Fvalue = mean( Fvalues )

    out = .outliergram( time_grid, Data, MBD_data, MEI_data, q_low, q_high,
                        Fvalue = Fvalue  )
  }

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Plot
  if( display )
  {
    # Setting up palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ), l = 60 )( length( out$ID_NO ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ), c = 150 )( length( out$ID_SO ) )
    # col_labels_outlying = scales::hue_pal( h = c( -10, 10 ), c = 1000, l = 20 )( length( out$ID_SO ) )

    dev.cur()
    par( mfrow = c( 1, 2 ) )

    # Plotting functional data
    matplot( time_grid, t( Data[ - out$ID_SO, ] ), type = 'l', lty = 1, ylim = range( Data ),
             col = col_non_outlying, ... )
    matplot( time_grid, t( Data[   out$ID_SO, ] ), type = 'l', lty = 1, lwd = 3, ylim = range( Data ),
             col = col_outlying, add = TRUE )

    # Adding text labels with curve ID
    w_spacing = diff( range( time_grid ) ) / ( 2 * length( out$ID_SO ) )

    for( iOut in seq_along( out$ID_SO ) )
    {
      text( time_grid[ 1 ] + ( 2 * iOut - 1 ) * w_spacing,
            Data[ out$ID_SO[ iOut ],
                  which.min( abs( time_grid - time_grid[ 1 ] - ( 2 * iOut - 1 ) * w_spacing ) ) ] +
              diff( range( Data[ out$ID_SO[ iOut ]  ] ) ) / 30,
            out$ID_SO[ iOut ],
            col = col_outlying[ iOut ] )
    }

    # Plotting outliergram

    # Upper parabolic limit
    grid_1D = seq( 0, 1, length.out = 100 )
    # plot( sort( out$MEI_data ), a_0_2 + a_1 * sort( out$MEI_data ) + a_0_2 * N^2 * sort( out$MEI_data )^2,
    plot( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2,
          lty = 2, type = 'l', col = 'darkblue', lwd = 2, ylim = c( 0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 ),
          main = 'Outliergram', xlab = 'MEI', ylab = 'MBD' )

    points( out$MEI_data[ - out$ID_SO ], out$MBD_data[ - out$ID_SO ],
            pch = 16, col = col_non_outlying )
    points( out$MEI_data[ out$ID_SO ], out$MBD_data[ out$ID_SO ], pch = 16, cex = 1.5, col = col_outlying )

    for( idOut in out$ID_SO )
    {
      text( out$MEI_data[ idOut ],
            out$MBD_data[ idOut ] + 0.5 / 30,
            idOut,
            col = col_outlying[ match( idOut, out$ID_SO ) ] )
    }

    # lower parabolic limit
    if( is.null( adjust ) )
    {
      lines( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 - out$Q_d3 - 1.5 * out$IQR_d,
             lty = 2, lwd = 2, col = 'lightblue' )
    } else {
      lines( grid_1D,  a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 - Fvalue * out$Q_d1,
             lty = 2, lwd = 2, col = 'lightblue' )
    }

  }

  return( out$ID_SO )
}

.outliergram = function( time_grid, Data, MBD_data = NULL, MEI_data = NULL, q_low = 0, q_high = 1, Fvalue = NULL )
{
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

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves
  ID_shape_outlier = ifelse( is.null( Fvalue ),
                             which( d >= Q_d3 + 1.5 * IQR_d ),
                             which( d >= Fvalue * Q_d1 ) )

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : nrow( Data ), ID_shape_outlier )

  # Low MEI curves will be checked for upward shift
  ID_non_outlying_Low_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] < 0.5 ) ]

  # High MEI curves will be checked for downward shift
  ID_non_outlying_High_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] >= 0.5 ) ]

  # Managing high MEI data
  for( iObs in ID_non_outlying_High_MEI )
  {
#     diff_min_vector = Data[ iObs, ] - apply( Data[ - iObs, ], 2, ifelse( q_low == 0,
#                                                                          min,
#                                                                          function( x )( quantile( x, probs = q_low ) ) ) )
    diff_min_vector = Data[ iObs, ] - apply( Data[ - iObs, ], 2, quantile, probs = q_low )

    min_diff_min = min( diff_min_vector )

    Data_tilde = Data[ iObs, ] - rep( any ( diff_min_vector < 0 ) * min_diff_min, ncol( Data ) )

    MBD_curr = MBD( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    MEI_curr = MEI( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    if( ifelse( is.null( Fvalue ),
                d_curr >= Q_d3 + 1.5 * IQR_d,
                d_curr >= Q_d1 * Fvalue ) )
    {
      ID_shape_outlier = c( ID_shape_outlier, iObs )
      ID_non_outlying = setdiff( ID_non_outlying,
                                 ID_non_outlying_High_MEI[ iObs ] )
    }
  }

  # Managing low MEI data
  for( iObs in ID_non_outlying_Low_MEI )
  {
#     diff_max_vector = Data[ iObs, ] - apply( Data[ - iObs, ], 2, ifelse( q_high == 1,
#                                                                          max,
#                                                                          function( x )( quantile( x, probs = q_high ) ) ) )
    diff_max_vector = Data[ iObs, ] - apply( Data[ - iObs, ], 2, quantile, probs = q_high )

    max_diff_max = max( diff_max_vector )

    Data_tilde = Data[ iObs, ] - rep( any ( diff_max_vector > 0 ) * max_diff_max, ncol( Data ) )

    MBD_curr = MBD( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    MEI_curr = MEI( rbind( Data[ - iObs, ], Data_tilde ) )[ N ]

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    if(ifelse( is.null( Fvalue ),
               d_curr >= Q_d3 + 1.5 * IQR_d,
               d_curr >= Q_d1 * Fvalue ) )
    {
      ID_shape_outlier = c( ID_shape_outlier, iObs )
      ID_non_outlying = setdiff( ID_non_outlying,
                                 ID_non_outlying_Low_MEI[ iObs ] )
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

  return( list( ID_SO = ID_shape_outlier,
                ID_NO = ID_non_outlying,
                MEI_data = MEI_data,
                MBD_data = MBD_data,
                Q_d3 = Q_d3,
                Q_d1 = Q_d1,
                IQR_d = IQR_d ) )

}


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

my_outliergram = function( time_grid, Data, MBD_data = NULL, MEI_data = NULL, p_check = 0.1, q_low = 0, q_high = 1, Fvalue = NULL )
{
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

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves
  ID_shape_outlier = ifelse( is.null( Fvalue ),
                             which( d >= Q_d3 + 1.5 * IQR_d ),
                             which( d >= Fvalue * Q_d1 ) )

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : nrow( Data ), ID_shape_outlier )

  # Low MEI curves will be checked for upward shift
  # ID_non_outlying_Low_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] < 0.5 ) ]
  ID_non_outlying_Low_MEI = ID_non_outlying[ which( MEI_data[ - ID_shape_outlier ] >=
                                                      quantile( MEI_data, probs = 1 - p_check ) ) ]

  # High MEI curves will be checked for downward shift
  # ID_non_outlying_High_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] >= 0.5 ) ]
  ID_non_outlying_High_MEI = ID_non_outlying[ which( MEI_data[ - ID_shape_outlier ] <=
                                                       quantile( MEI_data, probs = p_check ) ) ]
#   if( length( ID_non_outlying_High_MEI ) == 0 )
#   {
#     stop( ' Error: you provided a ')
#   }

  aux_function = function( ID )( min( Data[ ID, ] - apply( Data[ - ID, ], 2, quantile, probs = q_low ) ) )

  # Managing High MEI data
  min_diff_min = sapply( ID_non_outlying_High_MEI, aux_function )
#                       function( ID )( min( Data[ ID, ] -
#                                               apply( Data[ - ID, ],
#                                                      2, quantile, probs = q_low ) ) ) )


  ID_to_check = ID_non_outlying_High_MEI[ min_diff_min < 0 ]

  aux_function_MBD = function( ID, fun )( MBD( rbind( Data[ - ID, ],
                                               Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )
  aux_function_MEI = function( ID, fun )( MEI( rbind( Data[ - ID, ],
                                                    Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )

  if( length( ID_to_check ) > 0 )
  {
    Data_tilde = t( t( Data[ ID_to_check, ] ) - min_diff_min[ ID_to_check ] )


    MBD_curr = sapply( ID_to_check, aux_function )


    MEI_curr = sapply( ID_to_check, aux_function )

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    ID_out_extra = ID_to_check[ which( ifelse( is.null( Fvalue ),
                                               d_curr >= Q_d3 + 1.5 * IQR_d,
                                               d_curr >= Q_d1 * Fvalue ) ) ]

    ID_shape_outlier = c( ID_shape_outlier, ID_out_extra )
    ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra )
  }

  # Managing Low MEI data

  aux_function = function( ID )( max( Data[ ID, ] - apply( Data[ - ID, ], 2, quantile, probs = q_high ) ) )

  max_diff_max = sapply( ID_non_outlying_Low_MEI, aux_function )
#                          function( ID )( max( Data[ ID, ] -
#                                                  apply( Data[ - ID, ],
#                                                         2, quantile, probs = q_high ) ) ) )

  ID_to_check = ID_non_outlying_Low_MEI[ max_diff_max < 0 ]

  if( length( ID_to_check ) > 0 )
  {
    Data_tilde = t( t( Data[ ID_to_check, ] ) - max_diff_max[ ID_to_check ] )

    MBD_curr = sapply( ID_to_check, aux_function_MBD )
                       # function( ID )( MBD( rbind( Data[ - ID, ],
                                                   # Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] ) )

    MEI_curr = sapply( ID_to_check, aux_function_MEI )
                       # function( ID )( MEI( rbind( Data[ - ID, ],
                                                   # Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] ) )

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    ID_out_extra = ID_to_check[ which( ifelse( is.null( Fvalue ),
                                               d_curr >= Q_d3 + 1.5 * IQR_d,
                                               d_curr >= Q_d1 * Fvalue ) ) ]

    ID_shape_outlier = c( ID_shape_outlier, ID_out_extra )
    ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra )
  }


  return( list( ID_SO = ID_shape_outlier,
                ID_NO = ID_non_outlying,
                MEI_data = MEI_data,
                MBD_data = MBD_data ) )

}



dot_outliergram = function( time_grid, Data, MBD_data = NULL, MEI_data = NULL, q_low = 0, q_high = 1, Fvalue = NULL )
{
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

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves
  ID_shape_outlier = ifelse( is.null( Fvalue ),
                             which( d >= Q_d3 + 1.5 * IQR_d ),
                             which( d >= Fvalue * Q_d1 ) )

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

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    if( ifelse( is.null( Fvalue ),
                d_curr >= Q_d3 + 1.5 * IQR_d,
                d_curr >= Q_d1 * Fvalue ) )
    {
      ID_shape_outlier = c( ID_shape_outlier, iObs )
      ID_non_outlying = setdiff( ID_non_outlying,
                                 ID_non_outlying_High_MEI[ iObs ] )
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

    d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

    if(ifelse( is.null( Fvalue ),
               d_curr >= Q_d3 + 1.5 * IQR_d,
               d_curr >= Q_d1 * Fvalue ) )
    {
      ID_shape_outlier = c( ID_shape_outlier, iObs )
      ID_non_outlying = setdiff( ID_non_outlying,
                                 ID_non_outlying_Low_MEI[ iObs ] )
    }
  }

  return( list( ID_SO = ID_shape_outlier,
                ID_NO = ID_non_outlying,
                MEI_data = MEI_data,
                MBD_data = MBD_data ) )

}


par_outliergram = function( time_grid, Data, MBD_data = NULL, MEI_data = NULL, q_low = 0, q_high = 1, Fvalue = NULL )
{
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

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves
  ID_shape_outlier = ifelse( is.null( Fvalue ),
                             which( d >= Q_d3 + 1.5 * IQR_d ),
                             which( d >= Fvalue * Q_d1 ) )

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : nrow( Data ), ID_shape_outlier )

  # Low MEI curves will be checked for upward shift
  ID_non_outlying_Low_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] < 0.5 ) ]

  # High MEI curves will be checked for downward shift
  ID_non_outlying_High_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] >= 0.5 ) ]

  registerDoParallel( 8 )
  min_diff_min = foreach( ID = ID_non_outlying_High_MEI, .combine = 'c' ) %dopar% {
    min( Data[ ID, ] -
           apply( Data[ - ID, ],
                  2, min ) )
  }
  max_diff_max = foreach( ID = ID_non_outlying_Low_MEI, .combine = 'c' ) %dopar% {
    max( Data[ ID, ] -
           apply( Data[ - ID, ],
                  2, max ) )
  }

  ID_to_check_high = ID_non_outlying_High_MEI[ min_diff_min < 0 ]
  ID_to_check_low = ID_non_outlying_Low_MEI[ max_diff_max < 0 ]

  Data_tilde_high = t( t( Data[ ID_to_check_high, ] ) - min_diff_min[ ID_to_check_high ] )
  Data_tilde_low = t( t( Data[ ID_to_check_low, ] ) - max_diff_max[ ID_to_check_low ] )

  MBD_curr_high = foreach( ID = ID_to_check_high, .combine = 'c' ) %dopar% {
    MBD( rbind( Data[ - grep( ID, ID_non_outlying_High_MEI ), ],
                Data_tilde_high[ grep( ID, ID_to_check_high ), ] ) )[ N ]
  }
  MBD_curr_low = foreach( ID = ID_to_check_low, .combine = 'c' ) %dopar% {
    MBD( rbind( Data[ - grep( ID, ID_non_outlying_Low_MEI ), ],
                Data_tilde_low[ grep( ID, ID_to_check_low ), ] ) )[ N ]
  }

  MEI_curr_high = foreach( ID = ID_to_check_high, .combine = 'c'  ) %dopar% {
    MEI( rbind( Data[ - grep( ID, ID_non_outlying_High_MEI ), ],
                Data_tilde_high[ grep( ID, ID_to_check_high ), ] ) )[ N ]
  }
  MEI_curr_low = foreach( ID = ID_to_check_low, .combine = 'c' ) %dopar% {
    MEI( rbind( Data[ - grep( ID, ID_non_outlying_Low_MEI ), ],
                Data_tilde_low[ grep( ID, ID_to_check_low ), ] ) )[ N ]
  }

  d_curr_high = a_0_2 + a_1 * MEI_curr_high + N^2 * a_0_2 * MEI_curr_high^2 - MBD_curr_high
  d_curr_low = a_0_2 + a_1 * MEI_curr_low + N^2 * a_0_2 * MEI_curr_low^2 - MBD_curr_low

  ID_out_extra_high = ID_to_check_high[ which( ifelse( is.null( Fvalue ),
                                             d_curr_high >= Q_d3 + 1.5 * IQR_d,
                                             d_curr_high >= Q_d1 * Fvalue ) ) ]
  ID_out_extra_low = ID_to_check_low[ which( ifelse( is.null( Fvalue ),
                                                     d_curr_low >= Q_d3 + 1.5 * IQR_d,
                                                     d_curr_low >= Q_d1 * Fvalue ) ) ]

  ID_shape_outlier = c( ID_shape_outlier, ID_out_extra_high )
  ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra_high )
  ID_shape_outlier = c( ID_shape_outlier, ID_out_extra_low )
  ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra_low )


  return( list( ID_SO = ID_shape_outlier,
                ID_NO = ID_non_outlying,
                MEI_data = MEI_data,
                MBD_data = MBD_data ) )

}

