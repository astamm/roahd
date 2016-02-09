

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
