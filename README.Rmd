---
title: 'bpCausal: Bayesian Causal Inference with Time-Series Cross-Sectional Data'
output:
  html_document: default
  pdf_document: default
---

**R** package for `A Bayesain Alternative to the Synthetic Control Method`.
 
---

**Authors:** Xun Pang (Tsinghua); Licheng Liu (MIT); Yiqing Xu (Stanford).  

**Date:** Mar. 27, 2021

**Package:** bpCausal

**Version:** 0.0.1 (GitHub version). This package is still under development. 
Please report bugs!

---

## Contents

1. Installation

2. Instructions

3. Example


---
 
## Installation

The development version of the package can be installed from GitHub by typing 
the following commands:

```{r eval=FALSE}
install.packages('devtools', repos = 'http://cran.us.r-project.org') # if not already installed
devtools::install_github('liulch/bpCausal')
```

The core part of **bpcausal** is written in C++ to accelerate the computing speed, 
which depends on the packages **Rcpp** and **RcppArmadillo**. Pleases install them 
before running the functions in **bpCausal**.  


#### Notes on installation failures

1. For Rcpp, RcppArmadillo and MacOS "-lgfortran" and "-lquadmath" error, click [here]( http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/) for details.
2. Installation failure related to OpenMP on MacOS, click [here](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/) for a solution.
3. To fix these issues, try installing gfortran 6.1 from [here]( https://gcc.gnu.org/wiki/GFortranBinaries#MacOS clang4 R Binaries from https://github.com/coatless/r-macos-clang).
4. Clang error for Mac os Big Sur, you need to update Xcode, Xcode commandLine, and Clang.  

***

##  Instructions 

### Functional form

We begin with the description of the model to illustrate the syntax of the function 
`bpCausal()`. `bpCausal()` can be implemented to estimate an linear model with the 
following (reduced) functional form:

$$ y_{it} = \delta_{it} D_{it} + X_{it}^{\prime}\beta + Z_{it}^{\prime}\alpha_i 
+ A_{it}^{\prime}\Xi_t + \Gamma_i^{\prime}f_t + \varepsilon_{it} $$ 

where $D_{it}$ is a binary treatment indicator and $\delta_{it}$ 
represents heterogeneous treatment effect. $X_{it}$, $Z_{it}$ and $A_{it}$ are vectors 
of covariates that have constant, unit-level random and time-level random effects on the 
outcome, respectively. The random effects are assumed to have zero mean. Note that 
there can be overlapping among $X_{it}$, $Z_{it}$ and $A_{it}$ if some covarites have 
both constant and random effects. 

## Example 

We use the built-in dataset, `simdata`, to illustrate the functionality of `bpCausal()`.
This is a balanced panel data that contains 50 units, and the length of time periods is 
30. Among the units, 5 are exposed to the treatment from period 21 while the remaining 25 
units are never exposed to the treatment. Suppose we have 9 observed time-varying 
covariates, `x1` to `x9`. In the true data generating process, there are two latent 
factors, and the error term is sampled from i.i.d N(0, 1). We first load the data.

```{r data}
set.seed(1234)
library(bpCausal) 
data(bpCausal)
ls()
```

We implement the function `bpCausal()` to fit an linear model above with random effects 
at unit and (or) time level and latent factors. The `data` option specifies the dataset, 
which is a data frame object, the `index` option specifies the name of the unit and time 
index, and `Yname`, `Dname`, `Xname`, `Zname` and `Aname` specifies the names of the 
outcome variable, treatment indicator, and covariates respectively. Option `re` 
specifies the inclusion of two-way random effects. `ar1` is a logical flag indicating  
whether the factors (and the time-level random effects) follow ar1 process. `r` 
specifies the (upper bound of) number of factors. `Niter` is the number of MCMC draws 
and `burn` is the length of burn-in chains for MCMC simulations. `xlasso`, `zlasso`, 
`alasso` and `flasso` specify whether to shrink constant effects, unit or time level 
random effects, and factor loadings for factor selection, respectively. The other 
parameters are hyper parameters for the Gamma priors, and the default values ensure the 
priors are diffusive. For detailed explanation of the options, please refer to the help 
manual.  

```{r fit}
## with factors
out1 <- bpCausal(data = simdata, ## simulated dataset  
                 index = c("id", "time"), ## names for unit and time index
                 Yname = "Y", ## outcome variable
                 Dname = "D", ## treatment indicator  
                 Xname = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"), # covariates that have constant (fixed) effect  
                 Zname = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"), # covariates that have unit-level random effect  
                 Aname = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"), # covariates that have time-level random effect  
                 re = "both",   # two-way random effect: choose from ("unit", "time", "none", "both") 
                 ar1 = TRUE,    # whether the time-level random effects is ar1 process or jsut multilevel (independent)
                 r = 10,        # factor numbers 
                 niter = 15000, # number of mcmc draws
                 burn = 5000,   # burn-in draws 
                 xlasso = 1,    ## whether to shrink constant coefs (1 = TRUE, 0 = FALSE)
                 zlasso = 1,    ## whether to shrink unit-level random coefs (1 = TRUE, 0 = FALSE)
                 alasso = 1,    ## whether to shrink time-level coefs (1 = TRUE, 0 = FALSE)
                 flasso = 1,    ## whether to shrink factor loadings (1 = TRUE, 0 = FALSE)
                 a1 = 0.001, a2 = 0.001, ## parameters for hyper prior shrink on beta (diffuse hyper priors)
                 b1 = 0.001, b2 = 0.001, ## parameters for hyper prior shrink on alpha_i
                 c1 = 0.001, c2 = 0.001, ## parameters for hyper prior shrink on xi_t
                 p1 = 0.001, p2 = 0.001) ## parameters for hyper prior shrink on factor terms
```

The `bpCausal` package have two functions to summarize the posteriors. `coefSummary()` 
can be used to obtain summary statistics for posteriors of relevant parameters and 
`effSummary()` summaries the semi-parametric distribution of treatment effect, which 
is the difference between observed outcome under treatment and its corresponding 
posterior predictive distribution of the counterfactual outcomes. Options for individual 
level treatment effects and cumulative effects are also provided.  
 
```{r sum}
sout1 <- coefSummary(out1)  ## summary estimated parameters
eout1 <- effSummary(out1,   ## summary treatment effects
                    usr.id = NULL, ## treatment effect for individual treated units, if input NULL, calculate average TT
                    cumu = FALSE,  ## whether to calculate culmulative treatment effects
                    rela.period = TRUE) ## whether to use time relative to the occurence of treatment (1 is the first post-treatment period) or real period (like year 1998, 1999, ...)
``` 

We show the summary statistics for the posterior distributions of constant effects, where 
the first row is the intercept, and the by-period average treatment effect on the treated
: 

```{r beta}
sout1$est.beta 
``` 

```{r att}
eout1$est.eff
```

The following figure gives a graphical comparison between the true and estimated ATTs. The 
red line is the true ATT, the solid black line is the posterior mean, and the dashed lines 
indicate the 95 % posterior credible intervals. We find that the posterior means are very 
close to the true values. 

```{r attp}
x1 = c(-19:10)
y1 <- apply(matrix(simdata[which(simdata$treat==1),"eff"], 30, 5), 1, mean)

plot(x1, y1, type = "l", col = "red", ylim = c(-2, 12), 
    xlab = "Time", ylab = "ATT", cex.lab = 1.5)
abline(v = 0, lty = 3, col = "grey")
lines(x1, eout1$est.eff$estimated_ATT)
lines(x1, eout1$est.eff$estimated_ATT_ci_l, lty = 2)
lines(x1, eout1$est.eff$estimated_ATT_ci_u, lty = 2)
```

Finally, we show the results for factor number selection. The following figures are the 
posterior distribution of the square root of prior variance for each factor loading. If 
the posterior is spiked at zero, we are confident that we can remove the corresponding 
factor. If the posterior is bimodal, we are confident that there is variation among factor 
loadings, so the corresponding factor should be included to control for unobserved 
heterogeneity.

```{r factor}
xlim = c(-10,10)
dim(out1$gamma)
wg <- out1$wg
par(mfrow = c(2, 5), mar = c(3,3,2,1))
for (i in 1:10) {
    plot(density(wg[i,]), xlim = xlim, main = paste0("Loading ",i))
}
``` 














\pagebreak


---



$\square$



