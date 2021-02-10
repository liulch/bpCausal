\name{bpCausal-package}
\alias{bpCausal-package}
\docType{package}
\title{
Bayesian Causal Inference with TSCS Data
}
\description{
The package provides functions for estimating and making inference of treatment 
effects from time-series cross-sectional data with the Bayesian dynamic multilevel 
latent factor models.
}
\details{
The package provides functions for estimating and making inference of treatment 
effects from time-series cross-sectional data with the Bayesian dynamic multilevel 
latent factor models. Hierarchical shrinkage priors can be assigned for variable selection 
and determining the number of factors. It can accommodate both balanced and unbalanced 
panel data. 
}
\author{
Licheng Liu; Xun Pang; Yiqing Xu

Maintainer: Licheng Liu
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



