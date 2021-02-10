## summary function for Bayesian Causal 
## 09/08/2019: Xun made a change by adding niter to the other two funcetions
##
## 1. treated and conterfactual outcomes
## 2. estimated atts
## 3. estimated beta
## 4. multi-level part
## 5. dynamic part
## 6.

### The plot does not plot the fourth parameter 
coefSummary <- function(x,                 ## estimation results
                        burn = 0) {        ## burn-in length

    niter <- dim(x$sigma2)[2]
    out <- NULL

    ## ----------------------------- ##
    ## 1. beta
    beta_i <- x$beta
    if (!is.null(beta_i)) {
        p <- dim(beta_i)[1]

        beta_i <- matrix(c(beta_i[, (burn + 1):niter]), p, niter - burn)

        est_beta_mean <- apply(beta_i, 1, mean)
        est_beta_ci <- t(apply(beta_i, 1, quantile, c(0.025, 0.975)))
        data <- cbind.data.frame(est_beta_mean, est_beta_ci)

        names(data) <- c("mean", "ci_l", "ci_u")
        est.beta <- data

        out <- c(out, list(est.beta = est.beta))

    }


    ## ------------------------------ ##
    ## 2. multi-level coefficient 

    tr.unit.pos <- x$tr.unit.pos ## label treated units

    alpha_i <- x$alpha 
    wa_i <- x$wa
    if (!is.null(alpha_i)) {

        N <- dim(alpha_i)[1]
        p <- dim(alpha_i)[2]

        for (i in (burn + 1):niter) {
            alpha_i[,, i] <- alpha_i[,,i] * matrix(rep(wa_i[, i], each = N), N, p)
        }

        est_alpha_ci_l <- est_alpha_ci_u <- est_alpha_mean <- matrix(NA, N, p)

        for (i in 1:N) {
            sub_alpha <- matrix(alpha_i[i,,], p, niter)
            if (p == 1) {
                sub_alpha <- matrix(sub_alpha[, (burn + 1):niter], 1)
            } else {
                sub_alpha <- as.matrix(sub_alpha[, (burn + 1):niter])
            }
            est_alpha_mean[i,] <- apply(sub_alpha, 1, mean)
            est_alpha_ci_l[i,] <- apply(sub_alpha, 1, quantile, 0.025)
            est_alpha_ci_u[i,] <- apply(sub_alpha, 1, quantile, 0.975)
        }

        data <- NULL

        for (i in 1:p) {
            subdata <- cbind.data.frame(est_alpha_mean[, i], est_alpha_ci_l[, i], est_alpha_ci_u[, i])
            colnames(subdata) <- c("mean", "ci_l", "ci_u")
            subdata$id <- N:1
            subdata$tr <- 0
            subdata[tr.unit.pos, "tr"] <- 1
            subdata$tr <- as.factor(subdata$tr)
            data <- c(data, list(subdata))
        }

        est.alpha <- data
        out <- c(out, list(est.alpha = est.alpha))
    }


    ## ------------------------------ ##
    ## 3. time-varying coefficient 

    xi_i <- x$xi 
    wxi_i <- x$wxi
    
    if (!is.null(xi_i)) {
        TT <- dim(xi_i)[1]
        p <- dim(xi_i)[2]

        for (i in (burn + 1):niter) {
            xi_i[,, i] <- xi_i[,,i] * matrix(rep(wxi_i[, i], each = TT), TT, p)
        }

        est_xi_ci_l <- est_xi_ci_u <- est_xi_mean <- matrix(NA, TT, p)

        for (i in 1:TT) {
            sub_xi <- matrix(xi_i[i,,], p, niter)
            if (p == 1) {
                sub_xi <- matrix(sub_xi[, (burn + 1):niter], 1)
            } else {
                sub_xi <- as.matrix(sub_xi[, (burn + 1):niter])
            }
            est_xi_mean[i,] <- apply(sub_xi, 1, mean)
            est_xi_ci_l[i,] <- apply(sub_xi, 1, quantile, 0.025)
            est_xi_ci_u[i,] <- apply(sub_xi, 1, quantile, 0.975)
        }

        data <- NULL

        for (i in 1:p) {
            subdata <- cbind.data.frame(est_xi_mean[, i], est_xi_ci_l[, i], est_xi_ci_u[, i])
            colnames(subdata) <- c("mean", "ci_l", "ci_u")
            subdata$time <- 1:TT
            data <- c(data, list(subdata))
        }

        est.xi <- data
        out <- c(out, list(est.xi = est.xi))
    }


    ## ------------------------------ ##
    ## 4. ar1 coefficients

    p_i <- x$p

    if (!is.null(p_i)) {
        k <- dim(p_i)[1]

        p_i <- matrix(c(p_i[, (burn + 1):niter]), k, niter - burn)

        est_phi_mean <- apply(p_i, 1, mean)
        est_phi_ci <- t(apply(p_i, 1, quantile, c(0.025, 0.975)))
        data <- cbind.data.frame(est_phi_mean, est_phi_ci)

        names(data) <- c("mean", "ci_l", "ci_u")

        est.phi.xi <- est.phi.f <- NULL
        r <- k1 <- 0
        if (!is.null(x$xi)) {
            k1 <- dim(x$xi)[2]
            est.phi.xi <- data[1:k1,]
            out <- c(out, list(est.phi.xi = est.phi.xi))

        }
        if (!is.null(x$f)) {
            r <- dim(x$f)[2]
            est.phi.f <- data[(k1+1):(k1+r),]
            out <- c(out, list(est.phi.f = est.phi.f))
        }

    }

    return(out)
} 



## --------------------------------------- ##
## 5. observed and counterfactual outcomes  
## 6. plot by period att
## 7. plot cumulative att

effSummary <- function(x, 
                       usr.id = NULL,         ## individual effect, if left blank, all treated units will be used
                       burn = 0, 
                       cumu = FALSE,          ## whether to calculate cumulative effect
                       rela.period = TRUE) {  ## aggregate by time relative to treatment

    niter <- dim(x$sigma2)[2] 

    if (cumu) {
        rela.period <- TRUE
    }

    id.tr <- x$raw.id.tr
    time.tr <- x$time.tr
    rela.time.tr <- x$rela.time.tr  


    ## id indicator
    id.pos <- NULL
    unique.tr <- c(unique(id.tr))
    if (is.null(usr.id)) {
        id.pos <- 1:length(c(id.tr))  
    } else {
        if (sum(usr.id %in% unique.tr) != length(usr.id)) {
            stop("Some specified ids are not in treated group, please check input.\n")
        }
        id.pos <- which(c(id.tr) %in% usr.id)
    }

    yo_t <- x$yo_t
    yo_t <- yo_t[id.pos]

    time.tr <- time.tr[id.pos]
    rela.time.tr <- rela.time.tr[id.pos]

    yct_i <- x$yct
    yct_i <- matrix(c(yct_i[id.pos, (burn + 1):niter]), length(id.pos), niter - burn)

    ## mean observed and counterfactual
    count.tr <- NULL ## num of observations at each period

    if (rela.period) { ## relative to treatment occurence
        m_yo <- tapply(yo_t, rela.time.tr, mean)
        m_yct <- sapply(1:(niter - burn), function(i){tapply(yct_i[, i], rela.time.tr, mean)})

        count.tr <- as.numeric(table(rela.time.tr))
    } else { ## real time period
        m_yo <- tapply(yo_t, time.tr, mean)
        m_yct <- sapply(1:(niter - burn), function(i){tapply(yct_i[, i], time.tr, mean)})

        count.tr <- as.numeric(table(rela.time.tr))
    }


    ## outcomes -------------------

    m_yct_mean <- apply(m_yct, 1, mean)
    m_yct_ci_l <- apply(m_yct, 1, quantile, 0.025)
    m_yct_ci_u <- apply(m_yct, 1, quantile, 0.975)


    ## effect ---------------------

    eff_i <- matrix(rep(c(m_yo), niter - burn), length(c(m_yo)), niter - burn) - m_yct

    eff_mean <- apply(eff_i, 1, mean)
    eff_ci_l <- apply(eff_i, 1, quantile, 0.025)
    eff_ci_u <- apply(eff_i, 1, quantile, 0.975)


    data <- cbind.data.frame(m_yo, m_yct_mean, m_yct_ci_u, m_yct_ci_l, eff_mean, eff_ci_l, eff_ci_u)
    names(data) <- c("observed", "estimated_counterfactual", 
                     "counterfactual_ci_l", "counterfactual_ci_u",
                     "estimated_ATT", "estimated_ATT_ci_l", "estimated_ATT_ci_u")
    if(rela.period) {
        data$time <- sort(unique(rela.time.tr))
        data$count <- count.tr
    } else {
        data$time <- sort(unique(time.tr))
    }

    est.eff <- data


    ## cumulative effects ---------
    est.cumu <- NULL
    if (cumu) {
        relatime <- sort(unique(rela.time.tr))

        st.pos <- which(relatime == 1) ## start point 

        eff_sub_i <- matrix(c(eff_i[st.pos:length(relatime), ]), length(relatime) - st.pos +1)

        eff_cumu_i <- matrix(NA, length(relatime) - st.pos +1, niter - burn)

        count.tr.sub <- count.tr[st.pos:length(relatime)]

        eff_cumu_i[1, ] <- eff_sub_i[1, ]
        if (length(relatime) - st.pos >= 2) {
            for (j in 2:(length(relatime) - st.pos +1)) {
                eff_cumu_i[j, ] <- sapply(1:(niter - burn), function(i) {sum(eff_sub_i[1:j, i] * count.tr.sub[1:j])/sum(count.tr.sub[1:j])} ) * j
            }
        }

        eff_cumu_mean <- apply(eff_cumu_i, 1, mean)
        eff_cumu_ci_l <- apply(eff_cumu_i, 1, quantile, 0.025)
        eff_cumu_ci_u <- apply(eff_cumu_i, 1, quantile, 0.975)

        data <- cbind.data.frame(eff_cumu_mean, eff_cumu_ci_l, eff_cumu_ci_u)
        names(data) <- c("mean", "ci_l", "ci_u")

        data$count <- count.tr[st.pos:length(relatime)]
        data$time <- relatime[st.pos:length(relatime)]

        est.cumu <- data
    }
    


    ## average effetcs ------------

    t.post <- which(rela.time.tr > 0)

    eff_avg_i <- sapply(1:(niter - burn), function(i) {mean(yo_t[t.post] - yct_i[t.post, i])})

    eff_avg_mean <- mean(eff_avg_i)
    eff_avg_ci_l <- quantile(eff_avg_i, 0.025)
    eff_avg_ci_u <- quantile(eff_avg_i, 0.975)


    est.avg <- cbind(eff_avg_mean, eff_avg_ci_l, eff_avg_ci_u)
    colnames(est.avg) <- c("mean", "ci_l", "ci_u")


    out <- list(est.eff = est.eff, 
                est.avg = est.avg)

    if (!is.null(est.cumu)) {
        out <- c(out, list(est.cumu = est.cumu))
    }

    return(out)

}





