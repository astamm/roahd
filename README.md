
# roahd 

[![Build Status](https://travis-ci.org/ntarabelloni/roahd.svg?branch=dev)](https://travis-ci.org/ntarabelloni/roahd)

[![codecov](https://codecov.io/gh/ntarabelloni/roahd/branch/master/graph/badge.svg)](https://codecov.io/gh/ntarabelloni/roahd)


Package roahd (Robust Analysis of High-dimensional Data ) allows to use
a set of statistical tools for the exploratory analysis and robustification of
univariate and multivariate functional datasets through the use of depth-based
statistical methods.

A simple object-oriented (OO) representation of functional data object, 
allows to encapsulate the important features of functional datasets (like the 
grid of the dependent variable, the pointwise observations etc.) and to exploit
simple calls to customised functions which simplify the exploratory analysis.

A part of the package is specifically devoted to the computation of depths and 
other statistical indexes for functional data, which also are the core of the
visualization/robustification tools like functional boxplot and outliergram.

In the implementation of functions special attention was put to their efficiency,
so that they can be profitably used also for the analysis of high-dimensional 
datasets.

# Examples

See the Vignette.
