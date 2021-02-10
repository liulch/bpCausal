bayesLasso <- function(y, 
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
                       e2) {
    
    

    nobs <- dim(y)[1]
    tr.nobs <- dim(y.tr)[1]

    N <- max(id0)
    TT <- max(time)

    ## beta fit
    xfit <- matrix(0, nobs, 1)
    ## tilde{beta}_i(or t) fit
    zfit <- matrix(0, nobs, 1)
    afit <- matrix(0, nobs, 1)
    ## factor fit
    fefit <- matrix(0, nobs, 1)
    ## residual
    res <- matrix(0, nobs, 1)

    ## initialize
    if (k1 > 0) {
    	xfit <- X %*% beta0 
    }
    if (k2 > 0) {
        zfit <- getREfit(Z, Alpha0 * matrix(rep(c(omegaa0), each = N), N, k2), id0)
    } 
    if (k3 > 0) {
        afit <- getREfit(A, Xi0 * matrix(rep(c(omegaxi0), each = TT), TT, k3), time)
    }
    if (r > 0) {
    	fefit <- getFactorFit(Gamma0 * matrix(rep(c(omegag0), each = N), N, r), F0, id0, time) 
    }

    ## parameters
    ##  ------ 1. beta ----- ##
    beta <- beta0 
    sigma_beta2 <- sigma_beta20 
    B0 <- matrix(0, k1, k1)
    for (i in 1:k1) {
    	  B0[i, i] <- sigma_beta2[i]
    }
    lb2 <- lb20 

    
    ## ----- 2. tilde_alpha ----- ##
    Alpha <- Alpha0

    ## omega_alpha
    omegaa <- omegaa0
    sigma_alpha2 <- sigma_alpha20
    la2 <- la20

    ## ----- 3. tilde_xi ----- ##
    Xi <- Xi0

    ## omega_xi
    omegaxi <- omegaxi0
    sigma_xi2 <- sigma_xi20
    lxi2 <- lxi20

    
    ## ----- 4. loading ----- ##
    Gamma <- Gamma0 
    #if (re == "unit") {
        G0 <- matrix(0, k2 + r, k2 + r)
        for (i in 1:(k2 + r)) {
            G0[i, i] <- 1
        }
    #} else {
    #    G0 <- matrix(0, r, r)
    #    for (i in 1:r) {
    #        G0[i, i] <- 1
    #    }
    #}
    

    ## omega_gamma
    omegag <- omegag0
    sigma_gamma2 <- sigma_gamma20
    lg2 <- lg20 

    ## overall
    omega <- omega0

    A0 <- matrix(0, k2 + k3 + r, k2 + k3 + r)
    if (k2 > 0) {
        for (i in 1:k2) {
            A0[i, i] <- sigma_alpha2[i]
        }
    }
    if (k3 > 0) {
        for (i in 1:k3) {
            A0[i + k2, i + k2] <- sigma_xi2[i]
        }
    }
    if (r > 0) {
        for (i in 1:r) {
            A0[i + k2 + k3, i + k2 + k3] <- sigma_gamma2[i]
        }
    }


    ## ----- 4. factor ----- ##
    F <- F0

    Phi <- Phi0
    sigma_phi2 <- sigma_phi20 

    if (k3 == 0) {
        Tau_old <- Tau <- F0
        P0 <- matrix(0, r, r)  ## P for Phi
        for (i in 1:r) {
            P0[i, i] = sigma_phi2[i, 1]
        }
    } else {
        if (r > 0) {
            Tau_old <- Tau <- cbind(Xi0, F0)
            P0 <- matrix(0, k3 + r, k3 + r)
            for (i in 1:(k3 + r)) {
                P0[i, i] = sigma_phi2[i, 1]
            }
        } else {
            Tau_old <- Tau <- Xi0
            P0 <- matrix(0, k3 + r, k3 + r)
            for (i in 1:(k3 + r)) {
                P0[i, i] = sigma_phi2[i, 1]
            }
        }
        
    }

    ## 4 error term 
    res <- y - xfit - zfit - afit - fefit 
    sigma2 <- sampleSigmaE2(res, e1, e2) 

    ## store results 
    sb2_i <- beta_i <- matrix(0, k1, niter - burn)
    if (k2 > 0) {
        alpha_i <- array(0, dim = c(N, k2, niter - burn))
        swa2_i <- wa_i <- matrix(0, k2, niter - burn)
    } 
    if (k3 > 0) {
        xi_i <- array(0, dim = c(TT, k3, niter - burn))
        swxi2_i <- wxi_i <- matrix(0, k3, niter - burn)
    }
    if (r > 0) {
        gamma_i <- array(0, dim = c(N, r, niter - burn))
        f_i <- array(0, dim = c(TT, r, niter - burn))
        swg2_i <- wg_i <- matrix(0, r, niter - burn)
    }

    
    
    ## ar1 parameters
    ## if (re == "unit") {
    ##   p_i <- matrix(0, r, niter)
    ## } else {
        if (k3 + r > 0) {
            p_i <- matrix(0, k3 + r, niter - burn)
        }
    ## }
    
    ## error term 
    sigma2_i <- matrix(0, 1, niter - burn) 

    ## counterfactual
    yct_i <- matrix(0, tr.nobs, niter - burn)

    for (i in 1:niter) {
    	
      	## 1. update beta 
      	if (k1 > 0) {
      	    ## 1.1 sample beta
      	    res <- y - zfit - afit - fefit 
            beta <- sampleN(X, res, B0, sigma2)
            xfit <- X %*% beta
      	  
      	    ## 1.2 sample sigma_beta  	  
            if (xlasso == TRUE) {
                for (j in 1:k1) {
                  	pa1 <- sqrt(lb2/beta[j,1]^2)
                  	pa2 <- lb2
                    sigma_beta2[j, 1] <- 1 / rrinvgauss(pa1, pa2) 
            	    B0[j, j] <- sigma_beta2[j, 1] 
          	    }
                if (i > burn) {
                    sb2_i[, i - burn] <- sigma_beta2
                }
                

          	    ## 1.3 sample lambda_beta2 
          	    pa1 <- a1 + k1
          	    pa2 <- 1/(a2 + sum(sigma_beta2)/2)
          	    lb2 <- try(sampleG(pa1, pa2), silent = TRUE)
                if ('try-error' %in% class(lb2)) {
                    return(list(sb2 = sb2_i))
                }
            }
      	    
            ## 1.4 save results
            if (i > burn) {
                beta_i[, i - burn] <- beta 
            }
                   
      	}
    	
      	## 2. jointly update random effects and loadings
        #con1 <- re == "unit" && k + r > 0
        #con2 <- re == "time" && r > 0
        #if (con1 || con2) {
            ## jointly sample
        #    if (re == "unit") {
        #        res <- y - xfit 
        #        tX <- genTildeZ(X, F, omega, time)
        #        Gamma_all <- sampleAlpha(tX, res, G0, IDbreak, N, r, sigma2)
        #        Alpha <- as.matrix(Gamma_all[, 1:k])
        #        Gamma <- as.matrix(Gamma_all[, (k + 1):(k + r)])
                ## update refit 
        #        refit <- getREfit(X, Alpha * matrix(rep(c(omegaa), each = N), N, k), id0)
        #    } else {
        #        res <- y - xfit - refit
        #        tX <- genTildeZ(matrix(0, 1, 0), F, omegag, time)
        #        Gamma <- sampleAlpha(tX, res, G0, IDbreak, N, r, sigma2)
        #    }
        #}

        if (k2 + r > 0) {
            res <- y - xfit - afit
            if (k2 > 0) {
                tZ <- genTildeZ(Z, F, rbind(omegaa, omegag), time)
                Gamma_all <- sampleAlpha(tZ, res, G0, IDbreak, N, r, sigma2)
                Alpha <- as.matrix(Gamma_all[, 1:k2])
                #alpha_i[,,i] <- Alpha
                if (r > 0) {
                    Gamma <- as.matrix(Gamma_all[, (k2 + 1):(k2 + r)])
                    #if (i > burn) {
                    #    gamma_i[,,i - burn] <- Gamma
                    #}
                    
                }
                ## update zfit 
                zfit <- getREfit(Z, Alpha * matrix(rep(c(omegaa), each = N), N, k2), id0)
            } else {
                tZ <- genTildeZ(matrix(0, 1, 0), F, omegag, time)
                Gamma <- sampleAlpha(tZ, res, G0, IDbreak, N, r, sigma2)
                #if (i > burn) {
                #    gamma_i[,,i - burn] <- Gamma
                #}
                
            }
        }


        ## 3. jointly update random effects and factors
        #con1 <- re == "unit" && r > 0
        #con2 <- re == "time" && k + r > 0
        #if (con1 || con2) {
        #    if (re == "time") {
        #        res <- y - xfit
        #        res <- sortX(res, id2time)  
        #        tA <- sortX(genTildeZ(X, Gamma, omega, id0), id2time)
        #        Tau <- sampleXi(tA, res, Tau_old, Phi, TIMEbreak, T, sigma2) 
        #        Tau_old <- Tau
        #        Alpha <- as.matrix(Tau[, 1:k])
        #        F <- as.matrix(Tau[, (k + 1):(k + r)])  
        #    } else {
        #        res <- y - xfit - refit
        #        res <- sortX(res, id2time)  
        #        tA <- genTildeZ(matrix(0, 1, 0), Gamma, omegag, id)
        #        Tau <- sampleXi(tA, res, Tau_old, Phi, TIMEbreak, T, sigma2) 
        #        Tau_old <- Tau
        #        F <- Tau 
        #    }
        #}

        if (k3 + r > 0) {
            res <- y - xfit - zfit 
            res <- sortX(res, id2time)
            if (k3 > 0) {
                tA <- sortX(genTildeZ(A, Gamma, rbind(omegaxi, omegag), id0), id2time)
                Tau <- sampleXi(tA, res, Tau_old, Phi, TIMEbreak, TT, sigma2)
                Tau_old <- Tau
                Xi <- as.matrix(Tau[, 1:k3])
                #xi_i[,,i] <- Xi
                if (r > 0) {
                    F <- as.matrix(Tau[, (k3 + 1):(k3 + r)])
                    #f_i[,,i] <- F 
                }
            } else {
                tA <- genTildeZ(matrix(0, 1, 0), Gamma, omegag, id)
                Tau <- sampleXi(tA, res, Tau_old, Phi, TIMEbreak, TT, sigma2) 
                Tau_old <- Tau
                F <- Tau 
                #f_i[,,i] <- F
            } 
        }

        ## 4. update omega 
        if (k2 + k3 + r > 0) {
            ## 4.1 update omega
            res <- y - xfit
            #if (re == "unit") {
            #    tTau <- genTildeTau(X, Alpha, Gamma, F, id0, time, 1) 
            #} else {
            #    tTau <- genTildeTau(X, Alpha, Gamma, F, id0, time, 0) 
            #}
            tTau <- genTildeTau(Z, A, Alpha, Xi, Gamma, F, id0, time)
            omega <- sampleN(tTau, res, A0, sigma2)

            #omegaa <- as.matrix(omega[1:k])
            #omegag <- as.matrix(omega[(k+1):(k+r)])

            ## 4.2 sample sigma_omega      
            if (k2 > 0) {
                omegaa <- as.matrix(omega[1:k2])
                if (zlasso == TRUE) {
                    for (j in 1:k2) {
                        pa1 <- sqrt(la2/omegaa[j,1]^2)
                        pa2 <- la2
                        sigma_alpha2[j, 1] <- 1 / rrinvgauss(pa1, pa2) 
                        A0[j, j] <- sigma_alpha2[j, 1] 
                    }
                    ## sample lambda_alpha2 
                    pa1 <- b1 + k2
                    pa2 <- 1/(b2 + sum(sigma_alpha2)/2)
                    la2 <- try(sampleG(pa1, pa2), silent = TRUE)
                    
                    if ('try-error' %in% class(la2)) {
                        return(list(sa2 = sigma_alpha2))
                    }
                    
                    if (i > burn) {
                        swa2_i[, i - burn] <- sigma_alpha2
                    }
                    
                }
                ## permutation
                permu1 <- permute(omegaa, Alpha)
                omegaa <- permu1$omega
                Alpha <- permu1$Xi 
                ## store 
                if (i > burn) {
                    alpha_i[,,i - burn] <- Alpha
                    wa_i[, i - burn] <- omegaa
                }
                
            }

            if (k3 > 0) {
                omegaxi <- as.matrix(omega[(k2+1):(k2+k3)])
                if (alasso == TRUE) {
                    for (j in 1:k3) {
                        pa1 <- sqrt(lxi2/omegaxi[j,1]^2)
                        pa2 <- lxi2
                        sigma_xi2[j, 1] <- 1 / rrinvgauss(pa1, pa2) 
                        A0[j + k2, j + k2] <- sigma_xi2[j, 1] 
                    }
                    ## sample lambda_xi2 
                    pa1 <- c1 + k3
                    pa2 <- 1/(c2 + sum(sigma_xi2)/2)
                    lxi2 <- try(sampleG(pa1, pa2), silent = TRUE)
                    if ('try-error' %in% class(lxi2)) {
                        return(list(sxi2 = sigma_xi2))
                    }

                    if (i > burn) {
                        swxi2_i[, i - burn] <- sigma_xi2
                    }
                    
                }
                ## permutation
                permu2 <- permute(omegaxi, Xi)
                omegaxi <- permu2$omega
                Xi <- permu2$Xi 
                ## store
                if (i > burn) {
                    xi_i[,,i - burn] <- Xi
                    wxi_i[, i - burn] <- omegaxi
                }
                
            }
            
            if (r > 0) {
                omegag <- as.matrix(omega[(k2+k3+1):(k2+k3+r)])
                if (flasso == TRUE) {
                    for (j in 1:r) {
                        pa1 <- sqrt(lg2/omegag[j,1]^2)
                        pa2 <- lg2
                        sigma_gamma2[j, 1] <- 1 / rrinvgauss(pa1, pa2) 
                        A0[j + k2 + k3, j + k2 + k3] <- sigma_gamma2[j, 1] 
                    }
                    ## sample lambda_gamma2 
                    pa1 <- p1 + r
                    pa2 <- 1/(p2 + sum(sigma_gamma2)/2)
                    lg2 <- try(sampleG(pa1, pa2), silent = TRUE)
                    if ('try-error' %in% class(lg2)) {
                        return(list(sg2 = sigma_gamma2))
                    }
                    if (i > burn) {
                        swg2_i[, i - burn] <- sigma_gamma2
                    }
                    
                }
                ## permutation
                #permu3 <- permute(omegag, F)
                #omegag <- permu3$omega
                #F <- permu3$Xi 
                ## restriction on factor loadings 
                for (j in 1:r) {
                    con <- runif(1, 0, 1)
                    if (con > 1/4 && con <= 1/2) {
                        Gamma[, j] <- -Gamma[, j]
                        F[, j] <- -F[, j]
                    } else if (con > 1/2 && con <= 3/4) {
                        Gamma[, j] <- -Gamma[, j]
                        omegag[j, 1] <- -omegag[j, 1] 

                    } else if (con > 3/4) {
                         F[, j] <- -F[, j]
                         omegag[j, 1] <- -omegag[j, 1] 
                    }
                }
                
                #for (i in 1:r) {
                #    if (omegag[i, 1] * Gamma[i, i] <= 0) {
                #        omegag[i, 1] <- -omegag[i, 1]
                #    }
                #}
                #pos <- which(c(omegag) * diag(Gamma) <= 0)
                #omegag[pos, 1] <- - omegag[pos, 1]


                #permu3 <- permuteF(omegag, Gamma, F)
                #omegag <- permu3$omega
                #Gamma <- permu3$Gamma
                #F <- permu3$F
                ## store
                if (i > burn) {
                    f_i[,,i - burn] <- F 
                    gamma_i[,,i - burn] <- Gamma
                    wg_i[, i - burn] <- omegag
                }
                
            }
        }

        ## 5. permutation 
        #permu1 <- permute(omegaa, Alpha)
        #omegaa <- permu1$omega
        #Alpha <- permu1$Xi 

        #permu2 <- permuteF(omegag, Gamma, F)
        #omegag <- permu2$omega
        #Gamma <- permu2$Gamma
        #F <- permu2$F
        #permu2 <- permute(omegag, Gamma)
        #omegag <- permu2$omega
        #Gamma <- permu2$Xi


        ## wb_i[, i] <- omegaa
        ## wg_i[, i] <- omegag

        ## update AR(1) parameters 
      	if (k3 + r > 0) {
            ##  update Phi
            if (ar1 == TRUE) {
                Phi <- samplePhi(Tau, P0) 
                if (i > burn) {
                    p_i[, i - burn] <- Phi
                }
            }            
        }
        
      	## 8. error term 
        #if (k > 0) {
        #    if (re == "unit") {
        #        refit <- getREfit(X, Alpha * matrix(rep(c(omegaa), each = N), N, k), id0)
        #    } 
        #    else if (re == "time") {
        #        refit <- getREfit(X, Alpha * matrix(rep(c(omegaa), each = T), T, k), time)
        #    }
        #}

        if (k2 > 0) {
            zfit <- getREfit(Z, Alpha * matrix(rep(c(omegaa), each = N), N, k2), id0)
        }

        if (k3 > 0) {
            afit <- getREfit(A, Xi * matrix(rep(c(omegaxi), each = TT), TT, k3), time)
        }

        if (r > 0) {
            fefit <- getFactorFit(Gamma * matrix(rep(c(omegag), each = N), N, r), F, id0, time) 
        }

        res <- y - xfit - zfit - afit - fefit 
        sigma2 <- sampleSigmaE2(res, e1, e2) 
        if (i > burn) {
            sigma2_i[1, i - burn] <- sigma2 
        }
        

        ## 9. predict conterfactual
        y.ct <- X.tr %*% beta 
        if (k2 > 0) {
            y.ct <- y.ct + getREfit(Z.tr, Alpha * matrix(rep(c(omegaa), each = N), N, k2), id.tr)
        }
        if (k3 > 0) {
            y.ct <- y.ct + getREfit(A.tr, Xi * matrix(rep(c(omegaxi), each = TT), TT, k3), time.tr)
        }
        if (r > 0) {
            y.ct <- y.ct + getFactorFit(Gamma * matrix(rep(c(omegag), each = N), N, r), F, id.tr, time.tr)
        }
        #if (re == "unit") {
        #    y.ct <- y.ct + getREfit(X.tr, Alpha * matrix(rep(c(omegaa), each = N), N, k), id.tr)
        #} else {
        #    y.ct <- y.ct + getREfit(X.tr, Alpha * matrix(rep(c(omegaa), each = T), T, k), time.tr)
        #}
        #y.ct <- y.ct + getFactorFit(Gamma * matrix(rep(c(omegag), each = N), N, r), F, id.tr, time.tr)

        ## error term
        y.ct <- y.ct + as.matrix(rnorm(tr.nobs, sd = sqrt(sigma2)))

        if (i > burn) {
            yct_i[, i - burn] <- y.ct
        }
        

        if (i %% 100 == 0) {
            cat(paste("Simulated times: ", i, "\n", sep = ""))
        }
        
    }

    out <- list(yct = yct_i, sigma2 = sigma2_i)
    if (k1 > 0) {
        out <- c(out, list(beta = beta_i))
        if (xlasso == TRUE) {
            out <- c(out, list(sb2 = sb2_i))
        }
    }
    if (k2 > 0) {
        out <- c(out, list(alpha = alpha_i, wa = wa_i))
        if (zlasso == TRUE) {
            out <- c(out, list(swa2 = swa2_i))
        }
    }
    if (k3 > 0) {
        out <- c(out, list(xi = xi_i, wxi = wxi_i))
        if (alasso == TRUE) {
            out <- c(out, list(swxi2 = swxi2_i))
        }
    }
    if (r > 0) {
        out <- c(out, list(gamma = gamma_i, f = f_i, wg = wg_i))
        if (flasso == TRUE) {
            out <- c(out, list(swg2 = swg2_i))
        }
    }
    
    if (ar1 == TRUE) {
        if (r + k3 > 0) {
            if (k3 > 0) {
                out <- c(out, list(pxi = as.matrix(p_i[1:k3,])))
            }
            if (r > 0) {
                out <- c(out, list(pf = as.matrix(p_i[(k3+1):(k3+r),])))

            }
            out <- c(out, list(p = p_i))
        }
    }
    

    return(out)   
}

