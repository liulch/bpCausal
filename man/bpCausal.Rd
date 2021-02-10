\name{bpCausal}
\alias{bpCausal}
\title{Bayesian Causal inference with TSCS data}
\description{Bayesian Causal inference with TSCS data}
\usage{pblasso(data, index, Yname, Dname, Xname, Zname, Aname, 
               re = "unit", ar1 = TRUE, r, 
               niter = 10000, burn = 5000,
               xlasso, zlasso, alasso, flasso,
               a1 = 0.001, a2 = 0.001, b1 = 0.001, b2 = 0.001,
               c1 = 0.001, c2 = 0.001, p1 = 0.001, p2 = 0.001)  
}
\arguments{
\item{data}{a data frame (must be with a dichotomous treatment but balanced
  is not reqiored).}
\item{index}{a two-element string vector specifying the unit (group)
  and time indicators. Must be of length 2.}
\item{Yname}{outcome.}
\item{Dname}{treatment.}
\item{Xname}{covariates that have constant effects.}
\item{Zname}{covariates that have unit-level random effects.}
\item{Aname}{covariates that have time-level random effects.}
\item{re}{a string indicating whether additive unit or time fixed effects will be 
  imposed. Must be one of the following, 
  "none", "unit", "time", or "two-way". The default is "unit".}
\item{ar1}{a logical flag indicating whether time-varying parameters and factors are 
  AR(1) process or independent across periods.}
\item{r}{an integer specifying the (upper bound of) number of factors.}
\item{niter}{an integer specifying the iterations of MCMC run. Default value is 10000.}
\item{burn}{an integer specifying the length of burn in iterations of MCMC run. 
  Default value is 5000.}
\item{xlasso}{a logical flag indicating whether to assign constant coefficients
  hierachical shrinkage priors.}
\item{zlasso}{a logical flag indicating whether to assign unit-level random coefficients 
  hierachical shrinkage priors.}
\item{alasso}{a logical flag indicating whether to assign time-level random coefficients 
  hierachical shrinkage priors.}
\item{flasso}{a logical flag indicating whether to assign factor loadings 
  hierachical shrinkage priors for factor selection.}
\item{a1}{shape parameter for the Gamma prior for constant coefficients shrinkage.}
\item{a2}{rate parameter for the Gamma prior for constant coefficients shrinkage.}
\item{b1}{shape parameter for the Gamma prior for unit-level random coefficients 
  shrinkage.}
\item{b2}{rate parameter for the Gamma prior for unit-level random coefficients shrinkage.}
\item{c1}{shape parameter for the Gamma prior for time-level random coefficients 
  shrinkage.}
\item{c2}{rate parameter for the Gamma prior for unit-level random coefficients shrinkage.}
\item{p1}{shape parameter for the Gamma prior for factor loadings shrinkage.}
\item{p2}{rate parameter for the Gamma prior for factor loadings shrinkage.}
}
\details{ 
  \code{bpCausal} implements counterfactual estimators in TSCS data analysis. It 
  simulates counterfactual outcomes under control for observations under treatment 
  based on posterior predictive distributions.  
}
\value{
  \item{yct}{trace of simulated counterfactual outcomes under control for each 
  observation under treated.}
  \item{sigma2}{trace of simulated variance for the error term.}
  \item{beta}{trace of constant coefficients.}
  \item{alpha}{trace of unit-level random coefficients.}
  \item{xi}{trace of time-level random coefficients.}
  \item{gamma}{trace of factor loadings.}
  \item{f}{trace of factors.}
  \item{pxi}{trace of ar1 parameter for time-level random coefficients.}
  \item{pf}{trace of ar1 parameter for factors.}
  \item{wa}{trace of square roots for prior variance of unit-level random coefficients.}
  \item{wxi}{trace of square roots for prior variance of time-level random coefficients.}
  \item{wg}{trace of square roots for prior variance of factor loadings.}
  \item{sb2}{trace of prior variance of constant coefficients.}
  \item{swa2}{trace of prior variance of square roots for prior variance of unit-level 
  random coefficients.}
  \item{swxi2}{trace of prior variance of square roots for prior variance of time-level 
  random coefficients.}
  \item{swg2}{trace of prior variance of square roots for prior variance of factor 
  loadings.}
  \item{raw.id.tr}{a vector of names of treated units.}
  \item{id.tr}{a vector of integer indicating which units are under treated.}
  \item{time.tr}{a vector of periods (see \code{index}) for treated units.}
  \item{tr.unit.pos}{a vector of integer indicating which observations in the dataset 
  correspond to a treated unit.}
  \item{rela.time.tr}{a vector of period relative to the occurrence of treatment for each 
  observation that corresponds to a treated unit.}
}
\references{
  A Bayesian Alternative to Synthetic Control for Comparative Case Studies. 
  Pang et. al (2021). 
}
\examples{
library(bpCausal)
data(bpCausal)
out <- pblasso(data = simdata, index = c("id", "time"), 
               Yname = "Y", Dname = "D", 
               Xname = c("X1", "X2"), Zname = NULL, Aname = c("X1", "X2"), 
               re = "time", r = 10, niter = 20000)
}


