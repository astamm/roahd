bootstrapTest<-function( X, Y, B, ordering, normtype ){
  
  h<-length( X ) # number of components of X
  
  if( h!=length( Y ) ){
    
    stop( "Error in bootstrapTest: only multivariate functional datasets with the same number of components can be tested" )  # X and Y have different number of components for each unit, so they are not comparable
    
  }
  
  n_x<-dim( X[[1]] )[1] # number of functional observations of X
  n_y<-dim( Y[[1]] )[1] # number of functional observations of Y
  
  for ( i in 2:h ){
    
    if( dim( X[[i]] )[1]!= n_x )
      stop( "Error in bootstrapTest: X does not have the same number of observations in each component" )
    
  }
  
  for ( i in 2:h ){
    
    if( dim( Y[[i]] )[1]!= n_y )
      stop( "Error in bootstrapTest: Y does not have the same number of observations in each component" )
    
  }
  
  diff<-spearmanMatrix( X, ordering ) - spearmanMatrix( Y, ordering ) # difference between the sample Spearman matrices of X and Y
  phi<-norm( diff, type=normtype  ) # value of the test statistic
  v<-rbind( rep( 0, B ) ) # it's the vector in which the bootstrap replications of the test statistic under H0 are saved
  
  H_star_0<-list()
  Y_star_0<-list()
  
  for( i in 1:B ){
    
    y1<-sample.int( ( n_x+n_y ), size=n_x, replace=TRUE, prob=NULL ) # it's a sample of dimension n_x from the uniform distribution on the set {1,2,..,n_x+n_y}
    y2<-sample.int( ( n_x+n_y ), size=n_y, replace=TRUE, prob=NULL )
    
    for ( j in 1:h ){
      
      H_star_0[[j]]<- rbind( X[[j]], Y[[j]] )[y1,] 
      Y_star_0[[j]]<- rbind( X[[j]], Y[[j]] )[y2,]
      
    } # this for loop builds two multivariate functional datasets under H0 with through suitable resampling from X and Y (see Section 3.5 for more details)
    
    diff_0<-spearmanMatrix( H_star_0, ordering ) - spearmanMatrix( Y_star_0, ordering ) 
    v[i]<-norm( diff_0 , type=normtype) # bootstrap replication of the test statistic under H0
    
  }
  
  pvalue<- sum( v>=phi )/B
  
  output<-list( pvalue=pvalue, phi=phi )
  return(output)
  
}
