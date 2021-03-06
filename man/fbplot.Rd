% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fbplot.R
\name{fbplot}
\alias{fbplot}
\alias{fbplot.fData}
\alias{fbplot.mfData}
\title{Functional boxplot of univariate and multivariate functional data}
\usage{
fbplot(
  Data,
  Depths = "MBD",
  Fvalue = 1.5,
  adjust = FALSE,
  display = TRUE,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)

\method{fbplot}{fData}(
  Data,
  Depths = "MBD",
  Fvalue = 1.5,
  adjust = FALSE,
  display = TRUE,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)

\method{fbplot}{mfData}(
  Data,
  Depths = list(def = "MBD", weights = "uniform"),
  Fvalue = 1.5,
  adjust = FALSE,
  display = TRUE,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)
}
\arguments{
\item{Data}{the univariate or multivariate functional dataset whose
functional boxplot must be determined, in form of \code{fData} or
\code{mfData} object.}

\item{Depths}{either a vector containing the depths for each element of the
dataset, or:
\itemize{
\item \emph{univariate case}: a string containing the name of the method you
want to use to compute it. The default is \code{'MBD'}.
\item \emph{multivariate case}: a list with elements \code{def}, containing the
name of the depth notion to be used to compute depths (\code{BD} or
\code{MBD}), and \code{weights}, containing the value of parameter
\code{weights} to be passed to the depth function. Default is
\code{list(def = 'MBD', weights = 'uniform')}.
}

In both cases the name of the functions to compute depths must be available
in the caller's environment.}

\item{Fvalue}{the value of the inflation factor \eqn{F}, default is \code{F =
1.5}.}

\item{adjust}{either \code{FALSE} if you would like the default value for the
inflation factor, \eqn{F = 1.5}, to be used, or (for now \bold{only in the
univariate functional case}) a list specifying the parameters required by
the adjustment:
\itemize{
\item \code{N_trials}: the number of repetitions of the adjustment procedure
based on the simulation of a gaussian population of functional data, each
one producing an adjusted value of \eqn{F}, which will lead to the averaged
adjusted value \eqn{\bar{F}}. Default is 20.
\item \code{trial_size}: the number of elements in the gaussian population of
functional data that will be simulated at each repetition of the adjustment
procedure. Default is 8 * \code{Data$N}.
\item \code{TPR}: the True Positive Rate of outliers, i.e. the proportion of
observations in a dataset without amplitude outliers that have to be
considered outliers. Default is \code{2 * pnorm(4 * qnorm(0.25))}.
\item \code{F_min}: the minimum value of \eqn{F}, defining the left boundary
for the optimization problem aimed at finding, for a given dataset of
simulated gaussian data associated to \code{Data}, the optimal value of
\eqn{F}. Default is 0.5.
\item \code{F_max}: the maximum value of \eqn{F}, defining the right boundary
for the optimization problem aimed at finding, for a given dataset of
simulated gaussian data associated to \code{Data}, the optimal value of
\eqn{F}. Default is 5.
\item \code{tol}: the tolerance to be used in the optimization problem aimed at
finding, for a given dataset of simulated gaussian data associated to
\code{Data}, the optimal value of \eqn{F}. Default is \code{1e-3}.
\item \code{maxiter}: the maximum number of iterations to solve the
optimization problem aimed at finding, for a given dataset of simulated
gaussian data associated to \code{Data}, the optimal value of \eqn{F}.
Default is \code{100}.
\item \code{VERBOSE}: a parameter controlling the verbosity of the adjustment
process.
}}

\item{display}{either a logical value indicating whether you want the
functional boxplot to be displayed, or the number of the graphical device
where you want the functional boxplot to be displayed.}

\item{xlab}{the label to use on the x axis when displaying the functional
boxplot.}

\item{ylab}{the label (or list of labels for the multivariate functional
case) to use on the y axis when displaying the functional boxplot.}

\item{main}{the main title (or list of titles for the multivariate functional
case) to be used when displaying the functional boxplot.}

\item{...}{additional graphical parameters to be used in plotting functions.}
}
\value{
Even when used in graphical way to plot the functional boxplot, the
function returns a list of three elements:
\itemize{
\item \code{Depths}: contains the depths of each element of the functional
dataset.
\item \code{Fvalue}: is the value of F used to obtain the outliers.
\item \code{ID_out}: contains the vector of indices of dataset elements flagged
as outliers (if any).
}
}
\description{
This function can be used to perform the functional boxplot of univariate or
multivariate functional data.
}
\section{Adjustment}{


In the \bold{univariate functional case}, when the adjustment option is
selected, the value of \eqn{F} is optimized for the univariate functional
dataset provided with \code{Data}.

In practice, a number \code{adjust$N_trials} of times a synthetic population
(of size \code{adjust$tiral_size} with the same covariance (robustly
estimated from data) and centerline as \code{fData} is simulated without
outliers and each time an optimized value \eqn{F_i} is computed so that a
given proportion (\code{adjust$TPR}) of observations is flagged as outliers.
The final value of \code{F} for the functional boxplot is determined as an
average of \eqn{F_1, F_2, \dots, F_{N_{trials}}}. At each time step the
optimization problem is solved using \code{stats::uniroot} (Brent's method).
}

\examples{

# UNIVARIATE FUNCTIONAL BOXPLOT - NO ADJUSTMENT

set.seed(1)

N = 2 * 100 + 1
P = 2e2

grid = seq( 0, 1, length.out = P )

D = 10 * matrix( sin( 2 * pi * grid ), nrow = N, ncol = P, byrow = TRUE )

D = D + rexp(N, rate = 0.05)


# c( 0, 1 : (( N - 1 )/2), -( ( ( N - 1 ) / 2 ) : 1 ) )^4


fD = fData( grid, D )

dev.new()
oldpar <- par(mfrow = c(1, 1))
par(mfrow = c(1, 3))

plot( fD, lwd = 2, main = 'Functional dataset',
      xlab = 'time', ylab = 'values' )

fbplot( fD, main = 'Functional boxplot', xlab = 'time', ylab = 'values', Fvalue = 1.5 )

boxplot(fD$values[,1], ylim = range(fD$values), main = 'Boxplot of functional dataset at t_0 ' )

par(oldpar)

# UNIVARIATE FUNCTIONAL BOXPLOT - WITH ADJUSTMENT


set.seed( 161803 )

P = 2e2
grid = seq( 0, 1, length.out = P )

N = 1e2

# Generating a univariate synthetic gaussian dataset
Data = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ),
                             Cov = exp_cov_function( grid,
                                                     alpha = 0.3,
                                                     beta  = 0.4 ) )
fD = fData( grid, Data )

dev.new()
\donttest{
fbplot( fD, adjust = list( N_trials = 10,
                           trial_size = 5 * N,
                           VERBOSE = TRUE ),
                     xlab = 'time', ylab = 'Values',
                     main = 'My adjusted functional boxplot' )
}

# MULTIVARIATE FUNCTIONAL BOXPLOT - NO ADJUSTMENT

set.seed( 1618033 )

P = 1e2
N = 1e2
L = 2

grid = seq( 0, 1, length.out = 1e2 )

C1 = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
C2 = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )

# Generating a bivariate functional dataset of gaussian data with partially
# correlated components
Data = generate_gauss_mfdata( N, L,
                              centerline = matrix( sin( 2 * pi * grid ),
                                                   nrow = 2, ncol = P,
                                                   byrow = TRUE ),
                              correlations = rep( 0.5, 1 ),
                              listCov = list( C1, C2 ) )

mfD = mfData( grid, Data )

dev.new()
fbplot( mfD, Fvalue = 2.5, xlab = 'time', ylab = list( 'Values 1',
                                                       'Values 2' ),
        main = list( 'First component', 'Second component' ) )



}
\references{
\enumerate{
\item Sun, Y., & Genton, M. G. (2012). Functional boxplots. Journal of
Computational and Graphical Statistics.
\item Sun, Y., & Genton, M. G. (2012). Adjusted functional boxplots for
spatio-temporal data visualization and outlier detection. Environmetrics,
23(1), 54-64.
}
}
\seealso{
\code{\link{fData}}, \code{\link{MBD}}, \code{\link{BD}},
\code{\link{mfData}}, \code{\link{multiMBD}}, \code{\link{multiBD}}
}
