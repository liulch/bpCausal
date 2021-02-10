\name{coefSummary}
\alias{coefSummary}
\title{Summarize the posterior distribution for coefficients}
\description{Summarize the posterior distribution for coefficients}
\usage{coefSummary(x, burn = 0)} 
\arguments{
  \item{x}{a \code{\link{bpCausal}} object.}
  \item{burn}{an integer specifying the length of burn in iterations of (remaining) MCMC 
  run. Default is 0.}
}
\value{
  \item{est.beta}{posterior means and credible intervals for constant coefficients.}
  \item{est.alpha}{posterior means and credible intervals for unit-level random 
  coefficients.}
  \item{est.xi}{posterior means and credible intervals for time-level random 
  coefficients.}
  \item{est.phi.xi}{posterior means and credible intervals for ar1 parameters of 
  time-varying parameters.}
  \item{est.phi.f}{posterior means and credible intervals for ar1 parameters of factors.}
}
\references{
 A Bayesian Alternative to Synthetic Control for Comparative Case Studies. 
 Pang et. al (2021). 
}
\seealso{
  \code{\link{bpCausal}}
}


