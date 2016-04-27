
# Major fixes 

1) Modified the check of the grid provided to build fData objects. 
	Since support is provided only for evenly spaced grids, a check is needed before building an `fData` object.
	Before it was:

	```r
	all( abs( diff( unique( diff( grid  )  )  )  ) < 1e-14  )
	```

	Now it is:

	```r
		max( diff( unique( diff( grid  )  )  )  ) / diff( range( grid  )  ) < 1e-13
	```

	which is much more robust in practical cases.




# Minor fixes 

1) Fixed typos in documentation
