% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimuFwer.R
\name{SimuFdr}
\alias{SimuFdr}
\title{Simulates Gaussian data with a given correlation matrix and applies a FDR controlling procedure on the correlations.}
\usage{
SimuFdr(
  corr_theo,
  n = 100,
  Nsimu = 1,
  alpha = 0.05,
  stat_test = "empirical",
  method = "LCTnorm",
  Nboot = 1000,
  seed = NULL
)
}
\arguments{
\item{corr_theo}{the correlation matrix of Gaussien data simulated}

\item{n}{sample size}

\item{Nsimu}{number of simulations}

\item{alpha}{level of multiple testing}

\item{stat_test}{\describe{
  \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
  \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
  \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
  \item{'gaussian'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
}}

\item{method}{choice between 'LCTnorm' and 'LCTboot', developped by Cai & Liu (2016), 
'BH', traditional Benjamini-Hochberg (1995)'s procedure,
and 'BHboot', Benjamini-Hochberg (1995)'s procedure with bootstrap evaluation of pvalues}

\item{Nboot}{number of iterations for Monte-Carlo of bootstrap quantile evaluation}

\item{seed}{seed for the Gaussian simulations}
}
\value{
Returns a line vector containing estimated values for fwer, fdr, sensitivity, specificity and accuracy.
}
\description{
Simulates Gaussian data with a given correlation matrix and applies a FDR controlling procedure on the correlations.
}
\examples{
Nsimu <- 1000
n <- 100
p <- 10
corr_theo <- diag(1,p)
corr_theo[1,3] <- 0.5
corr_theo[3,1] <- 0.5
alpha <- 0.05
SimuFdr(corr_theo,n,Nsimu,alpha,stat_test='empirical',method='LCTnorm')
}
\references{
Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.

Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240.

Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, \url{https://tel.archives-ouvertes.fr/tel-01971574v1}.
}
\seealso{
ApplyFdrCor, SimuFwer
}
