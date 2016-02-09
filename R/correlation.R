

max.fData = function( fData, ..., which = FALSE )
{

  if( which )
  {
    return( data.frame( grid = fData$t0 + ( apply( fData$values, 1, which.max ) - 1 ) * fData$h,
                        value = apply( fData$values, 1, max ) ) )
  } else {
    return( values = apply( fData$values, 1, max ) )
  }

}


min.fData = function( fData, ..., which = FALSE )
{
  if( which )
  {
    return( data.frame( grid = fData$t0 + ( apply( fData$values, 1, which.min ) - 1 ) * fData$h,
                        value = apply( fData$values, 1, min ) ) )
  } else {
    return( values = apply( fData$values, 1, min ) )
  }

}

max_ordered = function( fData, gData )
{
  if( fData$P != gData$P ||
      fData$h != gData$h ||
      fData$t0 != gData$t0 ||
      fData$tP != gData$tP  )
  {
    stop( ' Error in max_ordered: provided fData objects have mismatching time grids')
  }

  if( fData$N != gData$N )
  {
    fData$N == 1 || gData$N == 1 ||
      stop( ' Error  in max_ordered: you must provide equally sized fData
            objects or at least one of them with 1 observation.')
  }

  return( max( fData ) - max( gData ) <= 0 )
}


area_under_curve = function( fData)
{
  return( rowSums( ( fData$values[ , - 1 ] +
                       fData$values[ , - fData$P ] ) / 2 * fData$h ) )
}

area_ordered = function( fData, gData )
{
  if( fData$P != gData$P ||
      fData$h != gData$h ||
      fData$t0 != gData$t0 ||
      fData$tP != gData$tP  )
  {
    stop( ' Error in area_ordered: provided fData objects have mismatching time grids')
  }

  if( fData$N != gData$N )
  {
    if( fData$N == 1 )
    {
      return( area_under_curve( gData - fData$values ) >= 0 )
    } else if( gData$N == 1 )
    {
      return( area_under_curve( fData - gData$values ) <= 0 )
    } else {
      stop( ' Error  in area_ordered: you must provide equally sized fData
            objects or at least one of them with 1 observation.')
    }
  } else {
    return( area_under_curve( fData - gData ) <= 0 )
  }
}


