## ----------------------------------------------------------------- ##
## a panel Bayesian lasso method for causal inference with tscs data ##
## ----------------------------------------------------------------- ##

## LLC PX XYQ 
## 2019.09.07

## data: a dataframe
## index: variable name for unit and time
## Yname: outcome
## Dname: treatment indicator
## Xname:
## Zname:
## Aname:
## re: unit, time, both
## r: factor numbers
## niter: number of iterations
## xlasso: a bool value
## zlasso:
## alasso:
## flasso:

bpCausal <- function(data, 
	                index, 
	                Yname, 
	                Dname, 
	                Xname,
	                Zname, 
	                Aname, 
	                re = "unit",
	                ar1 = TRUE, ## whether use ar1 or multilevel
	                r, 
	                placebo_period = 0,
	                niter = 10000,
	                burn = 5000,
	                xlasso,
	                zlasso,
	                alasso,
	                flasso,
	                a1 = 0.001,
					a2 = 0.001,
					b1 = 0.001,
					b2 = 0.001,
					c1 = 0.001,
					c2 = 0.001,
					p1 = 0.001,
					p2 = 0.001) {

	if (burn > niter) {
		stop("\"burn\" should not be greater than \"niter\".\n")
	}
	## raw data process 
	id <- index[1]
	time <- index[2]
	time_range <- range(data$time)
	
	if (placebo_period != 0) {
	  if ((placebo_period > time_range[2] - time_range[1])|| (placebo_period < 0)) {
	    warning("The placebo_period need to be in the range of time")
	  }
	
	  x_placebo_period <- tail(unique(simdata$time), placebo_period)
	  data <- data[!(data$time %in% x_placebo_period), ]
	}
	
    
    ## only save useful variables 
	data <- data[, unique(c(index, Yname, Dname, Xname, Zname, Aname))]
	
	
    ## remove missing values
	data <- na.omit(data)
	data <- data[order(data[, id], data[, time]),]

	raw.id.all <- unique(data[, id])

	## remove units always under treat
	rm.id <- tapply(data[, Dname], data[, id], sum)
	rm.id2 <- as.numeric(table(data[, id]))

	rm.id.pos <- which(rm.id == rm.id2)

	if (length(rm.id.pos) > 0) {
		rm.rela.pos <- which(data[, id] %in% raw.id.all[rm.id.pos])
		data <- data[-rm.rela.pos, ]
		raw.id.all <- raw.id.all[-rm.id.pos]
	}

	raw.id <- data[, id]

	data[, id] <- as.numeric(as.factor(data[, id]))
	data[, time] <- as.numeric(as.factor(data[, time]))

	## generate treatment indicators
	d <- data[, Dname]

	N <- max(data[, id])       ## number of units
	TT <- max(data[, time])    ## period length

	## normalize covariates
	varname <- unique(c(Xname, Zname, Aname))
	for (i in 1:length(varname)) {
		data[, varname[i]] <- data[, varname[i]] / sd(data[, varname[i]]) 
	}


	## treated status
	unit <- unique(data[, id])              ## all units
	tr <- unique(data[which(d==1), id])     ## treated units
	tr.pos <- which(data[, id] %in% tr)
	tr.unit.pos <- which(unit %in% tr)

	## gen two-way indicators for treated 
	id.tr <- as.matrix(data[tr.pos, id])
	raw.id.tr <- as.matrix(raw.id[tr.pos])
	time.tr <- as.matrix(data[tr.pos, time])
	d.tr <- as.matrix(data[tr.pos, Dname])

	## gen treatment indicator matrix
	Ntr <- length(tr)
	D.tr <- matrix(NA, TT, Ntr)
	I.tr <- matrix(0, TT, Ntr)
	for (i in 1:Ntr) {
		sub.pos <- which(c(id.tr) == tr[i])
		D.tr[c(time.tr)[sub.pos], i] <- c(d.tr)[sub.pos]
	}
	I.tr[which(!is.na(D.tr))] <- 1
    D.tr.old <- D.tr
	D.tr[which(I.tr == 0)] <- 0
	T.tr <- matrix(NA, TT, Ntr)
	for (i in 1:Ntr) {
		T.tr[, i] <- get_term(D.tr[, i], I.tr[, i]) 
	}
	rela.time.tr <- c(T.tr)[which(c(I.tr) == 1)]
	

	#T.tr[which(I.tr == 0)] <- NA
	#rela.time.tr <- na.omit(c(T.tr))

	#indicator <- cbind.data.frame(id.tr, time.tr, rela.time.tr)


	## split data: control for training model, treated for estimating effect
	y <- as.matrix(data[which(d == 0), Yname])
	y.tr <- as.matrix(data[tr.pos, Yname])

	## covar length
	k1 <- k2 <- k3 <- 0
	## gen covar matrix
	X.tr <- X <- A.tr <- A <- Z.tr <- Z <- matrix(0, 1, 0) 
	
	## gen X
	if (!is.null(Xname)) {
		X <- as.matrix(cbind(1, data[which(d == 0), Xname]))
		X.tr <- as.matrix(cbind(1, data[tr.pos, Xname])) 
		k1 <- dim(X)[2]
	} else {
		## only contain 1
		X <- as.matrix(rep(1, dim(y)[1]))
		X.tr <- as.matrix(rep(1, dim(y.tr)[1]))
		k1 <- 1
		xlasso <- 0
	}
	
	## gen Z
	if (!is.null(Zname)) {
		Z <- as.matrix(data[which(d == 0), Zname])
		Z.tr <- as.matrix(data[tr.pos, Zname])
		k2 <- length(Zname)
		if (re == "unit" || re == "both") {
			Z <- cbind(1, Z)
			Z.tr <- cbind(1, Z.tr)
			k2 <- k2 + 1
		}
	} else {
		if (re == "unit" || re == "both") {
			Z <- as.matrix(rep(1, dim(y)[1]))
			Z.tr <- as.matrix(rep(1, dim(y.tr)[1]))
			k2 <- 1
			zlasso <- 0
		}
	}

    ## gen A
	if (!is.null(Aname)) {
		A <- as.matrix(data[which(d == 0), Aname])
		A.tr <- as.matrix(data[tr.pos, Aname])
		k3 <- length(Aname)
		if (re == "time" || re == "both") {
			A <- cbind(1, A)
			A.tr <- cbind(1, A.tr)
			k3 <- k3 + 1
		}
	} else {
		if (re == "time" || re == "both") {
			A <- as.matrix(rep(1, dim(y)[1]))
			A.tr <- as.matrix(rep(1, dim(y.tr)[1]))
			k3 <- 1
			alasso <- 0
		}
	}
	
	# eff <- as.matrix(data[tr.pos, "eff"])
	# data2 <- data[order(data[, time], data[, id]),]

	               ## ------------------------- ## 
	## ------------- generate useful indicators ----------- ##
	               ## ------------------------- ## 
	data <- data[which(d == 0),c(id, time)] 
	## generate useful indicators and initialize some parameters
	data <- cbind(data, 1:dim(data)[1])
	data2 <- data[order(data[, time], data[, id]),]
	id2time <- as.matrix(data2[,dim(data2)[2]])

	data2 <- cbind(data2, 1:dim(data2)[1])
	data2 <- data2[order(data2[, id], data2[, time]),]

	time2id <- as.matrix(data2[,dim(data2)[2]])
	data2 <- data2[order(data2[, time], data2[, id]),]

	## save old 
	id.name <- id
	time.name <- time

	## vector
	id0 <- as.matrix(data[, id])
	id <- as.matrix(data2[, id])
	#time0 <- as.matrix(data2[, time])
	time <- as.matrix(data[, time])

	## break point
	IDbreak <- breakID(as.matrix(data[, id.name]))
	TIMEbreak <- breakID(as.matrix(data2[, time.name]))

	#Ntime <- as.matrix(as.numeric(table(time)))


	## initialized parameters 
	beta0 <- as.matrix(rnorm(k1)) ## k1 >= 1
	F0 <- Gamma0 <- Alpha0 <- Xi0 <- matrix(0, 1, 0)
	if (k2 > 0) {
		Alpha0 <- matrix(rnorm(k2*N), N, k2)
	} 
	if (k3 > 0) {
		Xi0 <- matrix(rnorm(k3*TT), TT, k3)
	}
	if (r > 0) {
		Gamma0 <- matrix(rnorm(r*N), N, r)
		F0 <- matrix(rnorm(r*T), TT, r)
	}
	
	omegaa0 <- matrix(1, k2, 1)
	omegaxi0 <- matrix(1, k3, 1)
	omegag0 <- matrix(1, r, 1)
	omega0 <- rbind(omegaa0, omegaxi0, omegag0)

	#if (re == "unit") {
	#	Phi0 <- matrix(0, r, 1)
	#	sigma_phi20 <- matrix(100, r, 1)
	#} else {
		Phi0 <- matrix(0, k3 + r, 1)
		sigma_phi20 <- matrix(100, k3 + r, 1)
	#}


	## relatively flat
	sigma_beta20 <- matrix(100, k1, 1)
	sigma_alpha20 <- matrix(100, k2, 1)
	sigma_xi20 <- matrix(100, k3, 1)
	sigma_gamma20 <- matrix(100, r, 1)

	lb20 <- 10
	la20 <- 10
	lxi20 <- 10
	lg20 <- 10

	## parameters in prior
	#a1 <- 0.001
	#a2 <- 0.001
	#b1 <- 0.001
	#b2 <- 0.001
	#c1 <- 0.001
	#c2 <- 0.001
	#p1 <- 0.001
	#p2 <- 0.001
	e1 <- 0.001
	e2 <- 0.001
	
	fit <- bayesLasso(y, 
                      y.tr,      
                      X,
                      X.tr,          
                      Z,
                      Z.tr,          
                      A,
                      A.tr,          
                      id,
                      id.tr,         
                      time,
                      time.tr, 
                      id0,
                      id2time,    
                      IDbreak,    
                      TIMEbreak,
                      xlasso,
                      zlasso,
                      alasso,
                      flasso, 
                      ar1, 
                      beta0,
                      Alpha0,
                      Xi0,
                      Gamma0,
                      F0,
                      omega0,
                      omegaa0,
                      omegaxi0,
                      omegag0,
                      Phi0,
                      sigma_beta20,
                      sigma_alpha20,
                      sigma_xi20,
                      sigma_gamma20,
                      sigma_phi20,
                      k1,
                      k2,
                      k3,
                      r,                
                      niter,
                      burn,
                      lb20,
                      la20,
                      lxi20,
                      lg20,
                      a1,
                      a2,
                      b1,
                      b2,
                      c1,
                      c2,
                      p1,
                      p2,
                      e1,
                      e2)
   
	
	## obserevd y and counterfactual
    #yo_t <- tapply(c(y.tr), c(time.tr), mean)
    yo_t <- y.tr
    
    ## yct_i <- fit$yct_i

    #yct_t <- matrix(0, TT, niter)
    #for (i in 1:niter) {
    #    yct_t[, i] <- tapply(c(yct_i[, i]), c(time.tr), mean)
    #}

    ## yct_t <- sapply(1:niter, function(i){return(tapply(c(yct_i[, i]), c(time.tr), mean))})

    return(c(fit, list(raw.id.tr = raw.id.tr,
    	               id.tr = id.tr, time.tr = time.tr, tr.unit.pos = tr.unit.pos,
    	               rela.time.tr = rela.time.tr, yo_t = yo_t)))

}



#################################
## support function
#################################
get_term <- function(d,
                     ii) {
    dd <- d
    iii <- ii
    ## dd <- dd[which(iii == 1)]
    first.pos <- min(which(iii == 1))
    if (first.pos != 1) {
        dd <- dd[-(1:(first.pos-1))]
        iii <- iii[-(1:(first.pos-1))]
    } 
    T <- length(dd) 
    if (0 %in% iii) {
        if (T > 1) {
            for (i in 1:(T-1)) {
                if (iii[i+1] == 0) {
                    dd[i+1] <- dd[i]
                }
            }
        }
    }

    d1 <- dd[1:(T-1)]
    d2 <- dd[2:T]
    if (sum(d1 == d2) == (T-1)) {
        term <- rep(NA, T)
    } else {
        change.pos <- which(d1 != d2) + 1
        change.length <- length(change.pos)
        term <- NULL    
        if (dd[1] == 0) {   
            for (i in 1:(change.length)) {
                if (i == 1) {
                    part.term <- (2 - change.pos[i]):0
                } else {
                    if (i %% 2 == 0) {
                        part.term <- 1:(change.pos[i] - change.pos[i-1])
                    } else {
                        part.term <- (change.pos[i-1] - change.pos[i] + 1):0
                    }
                }
                term <- c(term, part.term)
            }
        } else if (dd[1] == 1) {
            for (i in 1:(change.length)) {
                if (i == 1) {
                    part.term <- 1:(change.pos[i] - 1) 
                } else {
                    if (i %% 2 == 0) {
                        part.term <- (change.pos[i-1] - change.pos[i] + 1):0
                    } else {
                        part.term <- 1:(change.pos[i] - change.pos[i-1])
                    }
                }
                term <- c(term, part.term)
            }
        } 
        if (dd[change.pos[change.length]] == 0) {
            term <- c(term, rep(NA, (T - change.pos[change.length] + 1)))
        } else {
            term <- c(term, 1:(T - change.pos[change.length] + 1))
        }
    }
    ## term.all <- rep(NA, length(d))
    if (first.pos != 1) {
        term <- c(rep(NA, (first.pos - 1)), term)
    }
    return(term)
}



