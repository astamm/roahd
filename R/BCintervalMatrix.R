BCintervalMatrix<-function( data, ordering, B, alpha ){ 
  
  h<-length( data ) # number of components of the multivariate functional dataset
  inf<-diag( length( data ) ) # this is the matrix in which we will save the lower endpoints of the intervals
  sup<-diag( length( data ) ) # this is the matrix in which we will save the upper endpoints of the intervals
  n<-dim( data[[1]] )[1]
  
  for ( i in 2:h ){
    
    if ( dim( data[[i]] )[1] != n )
      stop(" Error in spearmanMatrix: the functional samples do not have the same number of observations")
    
  }
  
  for ( i in 1:(h-1) ){
    
    for( j in (i+1):h ){
      
      interval<-BCinterval( data[[i]], data[[j]], ordering, B, alpha ) # interval is a list containing the endpoints of the confidence interval (coverage probability 1-alpha) for the ij component of the Spearman matrix
      inf[i,j]<-interval$inf
      sup[i,j]<-interval$sup
      
    }
    
  }
  
  matrices<-list( inf=inf, sup=sup ) # the function returns a list containing the matrices of lower and upper endpoints.
  return(matrices)
  
}
