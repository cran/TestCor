% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ApplyFwer.R
\name{ApplyFwerCor_oracle}
\alias{ApplyFwerCor_oracle}
\title{Applies an oracle version of MaxTinfty procedure described in Drton & Perlman (2007) for correlation testing.}
\usage{
ApplyFwerCor_oracle(
  data,
  corr_theo,
  alpha = c(),
  stat_test = "empirical",
  method = "MaxTinfty",
  Nboot = 1000,
  stepdown = TRUE,
  vect = FALSE,
  logical = stepdown,
  arr.ind = FALSE
)
}
\arguments{
\item{data}{matrix of observations}

\item{corr_theo}{true matrix of correlations}

\item{alpha}{level of multiple testing (used if logical=TRUE)}

\item{stat_test}{\describe{
  \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
  \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
  \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
  \item{'2nd.order'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
}}

\item{method}{only 'MaxTinfty' implemented}

\item{Nboot}{number of iterations for Monte-Carlo of bootstrap quantile evaluation}

\item{stepdown}{logical, if TRUE a stepdown procedure is applied}

\item{vect}{if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix}

\item{logical}{if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
if \code{stepdown=TRUE} and \code{logical=FALSE}, returns a list of successive p-values.}

\item{arr.ind}{if TRUE, returns the indexes of the significant correlations, with repspect to level alpha}
}
\value{
Returns either \itemize{\item{the adjusted p-values, as a vector or a matrix, depending on \code{vect} (unavailable with stepdown)} \item{logicals indicating if the corresponding correlation is significant if \code{logical=TRUE}, as a vector or a matrix depending on \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
Oracle estimation of the quantile is used, based on the true correlation matrix
}
\description{
Applies oracle MaxTinfty procedure described in Drton & Perlman (2007) which controls asymptotically the FWER
for tests on a correlation matrix. It needs the true correlation matrix.
}
\examples{
n <- 100
p <- 10
corr_theo <- diag(1,p)
corr_theo[1,3] <- 0.5
corr_theo[3,1] <- 0.5
data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
# adjusted p-values:
(res <- ApplyFwerCor_oracle(data,corr_theo,stat_test='empirical',Nboot=1000,stepdown=FALSE))
# significant correlations, level alpha:
alpha <- 0.05
whichCor(res<alpha)
}
\references{
Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.

Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
}
\seealso{
ApplyFwerCor

maxTinftyCor, maxTinftyCor_SD
}
