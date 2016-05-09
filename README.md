
# roahd 

[![Build Status](https://travis-ci.org/ntarabelloni/roahd.svg?branch=dev)](https://travis-ci.org/ntarabelloni/roahd) [![codecov](https://codecov.io/gh/ntarabelloni/roahd/branch/dev/graph/badge.svg)](https://codecov.io/gh/ntarabelloni/roahd)

Package __roahd__ (_Robust Analysis of High-dimensional Data_) allows to use
a set of statistical tools for the _exploration_ and _robustification_ of
univariate and multivariate __functional datasets__ through the use of depth-based
statistical methods.


In the implementation of functions special attention was put to their efficiency,
so that they can be profitably used also for the analysis of high-dimensional 
datasets.

(_For a full-featured description of the package, please turn to the Vignette_)

## `fData` and `mfData` objects

A simple `S3` representation of functional data object, `fData`, 
allows to encapsulate the important features of univariate functional datasets (like the 
grid of the dependent variable, the pointwise observations etc.):

```r
# Grid representing the dependent variable
grid = seq( 0, 1, length.out = 100 )

# Pointwise-measurements of the functional dataset
Data = matrix( c( sin( 2 * pi * grid ), 
                  cos ( 2 * pi * grid ),
                  sin( 2 * pi * grid + pi / 4 ) ), ncol = 100, byrow = TRUE )

# S3 object encapsulating the univariate functional dataset            
fD = fData( grid, Data )

# S3 representation of a multivariate functional dataset
mfD = mfData( grid, list( 'comp1' = Data, 'comp2' = Data ) )
```
Also, this allows to exploit simple calls to customised functions which 
simplify the exploratory analysis:

```r
# Algebra of fData objects
fD + 1 : 100 
fD * 4

fD_1 + fD_2

# Subsetting fData objects (providing other fData objects)
fD[ 1, ]
fD[ 1, 2 : 4]

# Smaple mean and (depth-based) median(s)
mean( fD )
mean( fD[ 1, 10 : 20 ] )
median_fData( fD, type = 'MBD' )

# Plotting functions
plot( fD )
plot( mean( fD ), add = TRUE )

plot( fD[ 2:3, :] )
```


## Robust methods for functional data analysis

A part of the package is specifically devoted to the computation of depths and 
other statistical indexes for functional data:

  - Band Dephts and Modified Band Depths[^1][^2],
  - Modified band depths for multivariate functional data[^3][^4],
  - Epigraph and Hypograph indexes[^5],
  - Spearman[^6] and Kendall's[^7] correlation indexes for functional data.

These also are the core of the visualization/robustification tools like 
functional boxplot[^8] (`fbplot`) and outliergram[^9] (`outliergram`), allowing
the visualization and identification of amplitude/shape outliers.

Thanks to the functions for the simulation of synthetic functional datasets, 
both `fbplot` and `outliergram` procedures can be auto-tuned to the dataset 
at hand, in order to control the false outliers rate[^9].


## References

[^1]: Lopez-Pintado, S. and Romo, J. (2009). On the Concept of Depth for Functional Data, _Journal of the American Statistical Association_, 104, 718-734.
[^2]: Lopez-Pintado, S. and Romo. J. (2007). Depth-based inference for functional data, _Computational Statistics & Data Analysis_ 51, 4957-4968.
[^3]: Ieva, F., and A.M. Paganoni. (2013). Depth Measures for Multivariate Functional Data, _Communication in Statistics - Theory and Methods_ 42 (7): 1265–76.
[^4]: Tarabelloni, N. et al. (2015). Use of Depth Measure for Multivariate Functional Data in Disease Prediction: An Application to Electrocardiograph Signals. The _International Journal of Biostatistics_ 11 (2): 189–201
[^5]: Lopez-Pintado, S. and Romo, J. (2011). A Half-Region Depth for Functional Data.” _Computational Statistics & Data Analysis 55_ (4): 1679–95.
[^6]: Valencia, D., Romo, J. and Lillo, R. (2015). Spearman coefficient for functions. _Universidad Carlos III de Madrid technical report_ `http://EconPapers.repec.org/RePEc:cte:wsrepe:ws133329`.
[^7]: Valencia, D., Romo, J. and Lillo, R. (2015). A Kendall correlation coefficient for functional dependence. _Universidad Carlos III de Madrid technical report_ `http://EconPapers.repec.org/RePEc:cte:wsrepe:ws133228`.
[^8]: Sun, Y. and Genton, M. G. (2011). Functional Boxplots. _Journal of Computational and Graphical Statistics_ 20 (2): 316–34.
[^9]: Arribas-Gil, A., and J. Romo. (2014). Shape Outlier Detection and Visualization for Functional Data: The Outliergram. _Biostatistics_ 15 (4): 603–19.
