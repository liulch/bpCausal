\name{effSummary}
\alias{effSummary}
\title{Summarize the posterior distribution for treatment effects}
\description{Summarize the posterior distribution for treatment effects}
\usage{effSummary(x, 
                  usr.id = NULL,        
                  burn = 0, 
                  cumu = FALSE,          ## whether to calculate cumulative effect
                  rela.period = TRUE)   ## aggregate by time relative to treatment)
      } 
\arguments{
  \item{x}{a \code{\link{bpCausal}} object.}
  \item{usr.id}{a scalar or character indicating the treated unit for which the treatment 
  effect is summarized if user is interested in individual effect. If left blank 
  (\code{NULL}), all treated units will be used.}
  \item{burn}{an integer specifying the length of burn in iterations of (remaining) MCMC 
  run. Default is 0.}
  \item{cumu}{a logical flag indicating whether to calculate cumulative effects.}
  \item{rela.period}{a logical flag indicating whether to average treatment effect by 
  time relative to treatment.}
}
\value{
  \item{est.eff}{posterior means and credible intervals for treatment effect in each 
  period.}
  \item{est.avg}{posterior means and credible intervals for average treatment effect.}
  \item{est.cumu}{posterior means and credible intervals for cumulative treatment effect.}
}
\references{
 A Bayesian Alternative to Synthetic Control for Comparative Case Studies. 
 Pang et. al (2021). 
}
\seealso{
  \code{\link{bpCausal}}
}


