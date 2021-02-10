# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;


/*  ----------------- useful sub-functions --------------  */

/* Obtain factors and loadings by QR decomposition */
// [[Rcpp::export]]
List qr_factor (arma::mat F, arma::mat L) {
  
  int T = F.n_rows ;
  int N = L.n_rows ;
  int r = F.n_cols ;

  arma::mat factor(T, r, arma::fill::zeros) ;
  arma::mat lambda(N, r, arma::fill::zeros) ;
  arma::mat FE (T, N, arma::fill::zeros) ;

  arma::mat Q1 ; 
  arma::mat R1 ;
  arma::qr(Q1, R1, F) ;
  arma::mat Qf = Q1.cols(0, r-1) ;
  arma::mat Rf = R1.rows(0, r-1) ;

  // arma::mat RL = Rf * L.t() ;
  // arma::mat RL = L * Rf.t() ;

  arma::mat Q2 ;
  arma::mat R2 ;
  arma::qr(Q2, R2, Rf * L.t()) ;
  // arma::mat Ql = Q2.cols(0, r-1) ;
  // arma::mat Rl = R2.rows(0, r-1) ;

  factor = sqrt(double(T)) * Qf * Q2 ;
  // factor = Qf * Rfl.t() ;
  lambda = R2.t() / sqrt(double(T)) ;
  // factor = sqrt(double(T)) * Qf ;
  // lambda = L * Rf.t() / sqrt(double(T)) ;

  FE = factor * lambda.t() ;
  List result ;
  result["lambda"] = lambda ;
  result["factor"] = factor ;
  result["FE"] = FE ;
  return(result) ;
  
}

// sort matrix by index
// [[Rcpp::export]]
arma::mat sortX(arma::mat X, arma::mat index) {
  
  int p = X.n_cols ;
  int n = X.n_rows ;

  arma::mat sx(n, p, arma::fill::zeros) ;

  for (int i = 0; i < n; i++) {
    sx.row(i) = X.row(index(i, 0) - 1) ;
  }

  return(sx) ;
}

// generate a sub-matrix
// [[Rcpp::export]]
arma::mat getSubX(arma::mat X, arma::mat index) {
  int p = X.n_cols ;
  int n = index.n_rows ;
  int n2 = X.n_rows ;

  arma::mat subx(n, p, arma::fill::zeros) ;

  if (n2 == n) {
    subx = X ;
  }
  else {
    for (int i = 0; i < n; i++) {
      subx.row(i) = X.row(index(i, 0) - 1) ;
    }
  }

  return(subx) ;
}

// reshape from matrix to vector
// [[Rcpp::export]]
arma::mat reShape(arma::mat fit, 
                  arma::mat index, 
                  arma::mat idBreak) {

  int nobs = index.n_rows ;
  int N = fit.n_cols ;

  arma::mat rfit(nobs, 1, arma::fill::zeros) ; // reshape fit

  int n_start = 0 ;
  int n_end = 0 ;

  for (int i = 0; i < N; i++) {
    n_start = idBreak(i, 0) ;
    if (i < (N - 1)) {
      n_end = idBreak(i+1, 0) - 1 ;
    } 
    else {
      n_end = nobs - 1 ;
    }
    
    rfit.rows(n_start, n_end) = getSubX(fit.col(i), index.rows(n_start, n_end)) ;
  }
  return(rfit) ;  
}

// generate a vector of break point among units
// [[Rcpp::export]]
arma::mat breakID (arma::mat id) {
  
  int N = int(id.max()) ;
  int nobs = id.n_rows ;

  arma::mat idBreak(N, 1, arma::fill::zeros) ;

  // sort x by period
  int j = 1 ; // first period
  // int i_old = 0; // first position
  for (int i = 0; i < nobs; i++) {
  	if (id(i, 0) > j) {
  	  idBreak(j, 0) = i ;
  	  j = j + 1 ;
  	}
  	else {
  	  if (j == N) {
  		break ;
  	  }
  	}
  }

  return(idBreak) ;

}


// generate matrix useful for iterated beta
// [[Rcpp::export]]
arma::mat genXY(arma::mat x,
                arma::mat y,
                arma::mat idBreak) {

  // arma::cube XX(k, k, T, arma::fill::zeros) ;
  // arma::cube XY(k, 1, T, arma::fill::zeros) ;
 
  int k = x.n_cols ;
  int nobs = x.n_rows ;
  int T = idBreak.n_rows ;
  int k2 = y.n_cols ;

  arma::mat XY(k*T, k2, arma::fill::zeros) ;

  int n_start = 0 ;
  int n_end = 0 ;

  for (int i = 0; i < T; i++) {
    n_start = idBreak(i, 0) ;
    if (i < (T - 1)) {
      n_end = idBreak(i+1, 0) - 1 ;
    } 
    else {
      n_end = nobs - 1 ;
    }
    XY.rows(i*k, (i+1)*k - 1) = x.rows(n_start, n_end).t() * y.rows(n_start, n_end) ;

  }
  return(XY) ;
}


/*  ----------------- sub-functions for beta --------------  */
// generate covariance matrix
// [[Rcpp::export]]
arma::mat genCov1(arma::mat X, 
	                arma::mat invCov0,
	                double sigma2) {
  
  int k = X.n_cols ;
  double invsigma2 = 1/sigma2 ;
  arma::mat invCov1(k, k, arma::fill::zeros) ;

  invCov1 = invsigma2 * X.t() * X + invCov0 ;
  return(arma::inv(invCov1)) ;
}


// generate mean value matrix 
// [[Rcpp::export]]
arma::mat genMu(arma::mat X,
	              arma::mat Y,
	              arma::mat Cov1,
	              double sigma2) {
  
  int k = X.n_cols ;
  double invsigma2 = 1/sigma2 ;
  arma::mat mu(k, 1, arma::fill::zeros) ;

  mu = Cov1 * (invsigma2 * X.t() * Y) ;
  return(mu) ;
}


// sample multivariate normal
// [[Rcpp::export]]
arma::mat sampleN(arma::mat X,
  	              arma::mat Y,
  	              arma::mat Cov0,
  	              double sigma2) {  

  // arma::mat invCov0 = arma::inv(Cov0) ;
  // invCov0(0, 0) = 0 ;

  int k = X.n_cols ;
  arma::mat invCov0(k, k, arma::fill::zeros) ;
  for (int j = 1; j < k; j++) {
    invCov0(j, j) = 1 / Cov0(j, j) ; 
  }

  arma::mat S = genCov1(X, invCov0, sigma2) ;
  arma::mat M = genMu(X, Y, S, sigma2) ;

  // int k = X.n_rows ;
  arma::mat beta(k, 1, arma::fill::zeros) ;
  if (k > 1) {
  	beta = arma::mvnrnd(M.col(0), S) ;
  }
  else {
  	beta(0, 0) = R::rnorm(M(0, 0), sqrt(S(0, 0))) ; 
  }
  return(beta) ;

}


              /* ----------------------------- */    
/*------------ unit-level re and factor loadings ---------------*/
              /* ----------------------------- */ 

// sample conditional multivariate normal
// [[Rcpp::export]]
arma::mat sampleCN(arma::mat M,   // mean vector
                   arma::mat S,   // covariance matrix
                   arma::mat m) { // conditional on outcome

  int n1 = M.n_rows ;
  int n2 = m.n_rows ;
  int n3 = n1 - n2 ;

  arma::mat s11 = S.submat(0, 0, n3 - 1, n3 - 1) ;
  arma::mat s22 = S.submat(n3, n3, n1 - 1, n1 - 1) ;
  arma::mat s12 = S.submat(0, n3, n3 - 1, n1 - 1) ;
  arma::mat s21 = S.submat(n3, 0, n1 - 1, n3 - 1) ;

  arma::mat mu1 = M.rows(0, n3 - 1) ;
  arma::mat mu2 = M.rows(n3, n1 - 1) ;

  arma::mat invs22 = arma::inv(s22) ;
  arma::mat mu = mu1 + s12 * invs22 * (m - mu2) ;
  arma::mat s = s11 - s12 * invs22 * s21 ;

  arma::mat sv(n3, 1) ;
  sv.fill(-1) ;

  // double q = 0 ;

  if (n3 == 1) {
    // while(q <= 0) {
      // q = R::rnorm(mu(0, 0), sqrt(s(0, 0))) ;
    // }
    sv(0, 0) = R::rnorm(mu(0, 0), sqrt(s(0, 0))) ;
  } 
  else {
    // while (q <= 0) {
      sv = arma::mvnrnd(mu.col(0), s) ;
      // q = sv(n3 - 1, 0) ;
    // }
  }
  return(sv) ;
}


// sampling unit-level re and factor loadings
// [[Rcpp::export]]
arma::mat sampleSubAlpha(arma::mat M,  // mean vector
                         arma::mat S, // covariance matrix
                         int r, // number of factors
                         int id // id of the unit: to ensure the restrictions on loadings
                         ) {

  int k = M.n_rows ; // unit-level random effects plus factor loadings
  arma::mat sV(k, 1, arma::fill::zeros) ;
  // double q = 0 ;

  if (r == 0) { // no factors 
    if (k == 1) {
      sV(0,0) = R::rnorm(M(0, 0), sqrt(S(0, 0))) ;  
    }
    else {
      sV = arma::mvnrnd(M.col(0), S) ;
    }
  }
  else {
    if (k == r) { // only factors 
      if (r == 1) {
        if (id != 1) {
          sV(0,0) = R::rnorm(M(0, 0), sqrt(S(0, 0))) ;
        } 
        else {
          // while(q <= 0) {
          //   q = R::rnorm(M(0, 0), sqrt(S(0, 0))) ;
          // }
          //while(sV(0,0) <= 0) {
            sV(0,0) = R::rnorm(M(0, 0), sqrt(S(0, 0))) ;
          //}
        }
      }
      else { // r > 1
        if (id > r) {
          sV = arma::mvnrnd(M.col(0), S) ;
        }
        else { // 
          //if (id == 1) {
            // while(q <= 0) {
            // q = R::rnorm(M(0, 0), sqrt(S(0, 0))) ;
            //}
            //while(sV(0,0) <= 0) {
              // sV(0,0) = R::rnorm(M(0, 0), sqrt(S(0, 0))) ;
            //  sV(0, 0) = sampleCN(M, S, sV.rows(1, r - 1)) ;
            //}
          //}
          //else {
            if (id <= r - 1) {
              //while(sV(id - 1, 0) <= 0) {
                sV.rows(0, id - 1) = sampleCN(M, S, sV.rows(id, r - 1)) ; 
              //}
            } else {
              // while(q <= 0) {
              //if (id == r) {
                //while(sV(r - 1, 0) <= 0) {  
                  sV = arma::mvnrnd(M.col(0), S) ;
                //}
              //} else {
              //  sV = arma::mvnrnd(M.col(0), S) ;
                // q = sV(r-1, 0) ;
              //}
              // }
            }
          //}
        }
      }
    }
    else { // both factors and unit-level random effects 
      if (id > r) {
        sV = arma::mvnrnd(M.col(0), S) ;
      }
      else {
        if (id == r) {
          //while(sV(k - 1, 0) <= 0) {
            sV = arma::mvnrnd(M.col(0), S) ;
          //}
        }
        // sV(k - r + id - 1, 0) = 1 ;
        else {
          //while(sV(k - r + id - 1, 0) <= 0) {
            sV.rows(0, k - r + id - 1) = sampleCN(M, S, sV.rows(k - r + id, k - 1)) ;
          //}
        }
        // sV = arma::mvnrnd(M.col(0), S) ;
      }
    }
  }
  return(sV) ;
}


// iterated alpha sampling for each unit
// [[Rcpp::export]]
arma::mat iterGenAlpha(arma::mat XX,
  	                   arma::mat XY,
  	                   arma::mat A0,
  	                   int N,
  	                   int k,
  	                   int r,
  	                   double sigma2) {
  // int k1 = k - r ;

  arma::mat invA0 = arma::inv(A0) ;

  arma::mat bara(k, 1, arma::fill::zeros) ;
  arma::mat A1(k, k, arma::fill::zeros) ;

  // store all random effects
  arma::mat Alpha(N, k, arma::fill::zeros) ;
  // store factor loadings
  arma::mat L(N, r, arma::fill::zeros) ;

  arma::mat a(k, 1, arma::fill::zeros) ;

  for (int i = 0; i < N; i++) {
  	A1 = arma::inv(XX.rows(i*k, (i+1)*k - 1) / sigma2 + invA0) ;
  	bara = A1 * XY.rows(i*k, (i+1)*k - 1) / sigma2 ;
  	// we need to add restrictions on the factor loadings 
  	a = sampleSubAlpha(bara, A1, r, i+1) ; // i+1 unit indicator
  	Alpha.row(i) = a.t() ;
  }

  // if (r > 0) {
  //	L = Alpha.cols(k1, k - 1) ;
  // 	Alpha.shed_cols(k1, k - 1) ;
  // }
  // List result ;
  // if (k1 > 0) {
  // 	result["Alpha"] = Alpha ;
  // }
  // if (r > 0) {
  // 	result["L"] = L ;
  // }
  return(Alpha) ;
}


// iterated alpha sampling for each unit: a N*k matrix 
// [[Rcpp::export]]
arma::mat sampleAlpha(arma::mat X,
  	                  arma::mat y,
  	                  arma::mat A0,
  	                  arma::mat idBreak, // indicator vector
      			          int N,
      			          int r,
      			          double sigma2) {
  
  int k = X.n_cols ;

  arma::mat XY = genXY(X, y, idBreak) ;
  arma::mat XX = genXY(X, X, idBreak) ;

  arma::mat Alpha = iterGenAlpha(XX, XY, A0, N, k, r, sigma2) ;

  return(Alpha) ;
}

                   /* ------------------------------------- */    
      /*------------ state space modeling for xi and factors ---------------*/
                   /* ------------------------------------- */


/* ------------ sub-functions ------------ */
// permutation: generate binary value
// [[Rcpp::export]]
arma::mat genBi(int k) {
  arma::mat Bi(k, 1, arma::fill::ones) ;
  double rn = 0 ;
  for (int i = 0; i < k; i++) {
  	rn = R::rnorm(0.0, 1.0) ;
  	if (rn < 0) {
  		Bi(i, 0) = -1 ;
  	}
  }
  return(Bi) ;
}


// permutation result
// [[Rcpp::export]]
List permute(arma::mat omega,
	           arma::mat xi) {
  
  int k = xi.n_cols ;
  int T = xi.n_rows ;
  arma::mat coef = genBi(k) ; // k*1

  arma::mat omega_p = omega % coef ;
  arma::mat xi_p = xi % repmat(coef.t(), T, 1) ;

  List result ;
  result["omega"] = omega_p ;
  result["Xi"] = xi_p ;
  return(result) ;
}


// permutation result2
// [[Rcpp::export]]
List permuteF(arma::mat omega,
              arma::mat gamma,
              arma::mat F) {
  
  int r = gamma.n_cols ;

  arma::mat omega2 = omega ;
  arma::mat gamma2 = gamma ;
  arma::mat F2 = F ;

  double a1 = 0 ;
  double a2 = 0 ;

  for (int j = 0; j < r; j++) {
    a1 = R::runif(0.0, 1.0) ;
    a2 = R::runif(0.0, 1.0) ;
    if (a1 < 1/3) {
      if (a2 > 1/2) {
        gamma2.col(j) = gamma2.col(j) * (-1) ;
        F2.col(j) = F2.col(j) * (-1) ;
      }
    } 
    else if (a1 > 2/3) {
      if (a2 > 1/2) {
        omega2.row(j) = omega2.row(j) * (-1) ;
        F2.col(j) = F2.col(j) * (-1) ;
      }
    }
    else {
      if (a2 > 1/2) {
        omega2.row(j) = omega2.row(j) * (-1) ;
        gamma2.col(j) = gamma2.col(j) * (-1) ;
      }
    }
  }

  List out ;
  out["omega"] = omega2 ;
  out["Gamma"] = gamma2 ;
  out["F"] = F2 ;
  return(out) ;
  
}


// generate tilde A
// [[Rcpp::export]]
arma::mat genTildeA(arma::mat A,
                    arma::mat gamma,
                    arma::mat omega,
                    arma::mat id, // unit index
                    int r) {
  
  int nobs = id.n_rows ;
  int k = A.n_cols ;

  arma::mat newgamma(nobs, r, arma::fill::zeros) ;
  arma::mat newA = A ;

  if (k > 0) {
    if (r > 0) {
      for (int i = 0; i < nobs; i++) {
        newgamma.row(i) = gamma.row(id(i, 0) - 1) ;
      }
      // newA.insert_cols(k-1, newgamma) ;
    }
    return(arma::join_rows(newA, newgamma) % repmat(omega.t(), nobs, 1)) ;
  }
  else {
    for (int i = 0; i < nobs; i++) {
      newgamma.row(i) = gamma.row(id(i, 0) - 1) ;
    }
    return(newgamma % repmat(omega.t(), nobs, 1)) ;
  }
}


// generate tilde Z
// [[Rcpp::export]]
arma::mat genTildeZ(arma::mat Z,
                    arma::mat f,
                    arma::mat omega,
                    arma::mat time) { // time index
  
  int nobs = time.n_rows ;
  int k = Z.n_cols ;
  int r = f.n_cols ;

  arma::mat newf(nobs, r, arma::fill::zeros) ;
  arma::mat newZ = Z ;

  if (k > 0) {
    if (r > 0) {
      for (int i = 0; i < nobs; i++) {
        newf.row(i) = f.row(time(i, 0) - 1) ;
      }
      // newZ.insert_cols(k-1, newf % repmat(omegaf.t(), nobs, 1)) ;
    }
    return(arma::join_rows(newZ, newf) % repmat(omega.t(), nobs, 1)) ;
  }
  else {
    for (int i = 0; i < nobs; i++) {
      newf.row(i) = f.row(time(i, 0) - 1) ;
    }
    return(newf % repmat(omega.t(), nobs, 1)) ;
  }
}

// generate tilde Tau
// [[Rcpp::export]]
arma::mat genTildeTau(arma::mat Z,
                      arma::mat A,
                      arma::mat Alpha,
                      arma::mat Xi,
                      arma::mat Gamma,
                      arma::mat F,
                      arma::mat id0,
                      arma::mat time) {
  
  int nobs = id0.n_rows ;
  int k1 = Z.n_cols ;
  int k2 = A.n_cols ;
  int r = F.n_cols ;

  arma::mat ze(nobs, k1, arma::fill::zeros) ;
  arma::mat ae(nobs, k2, arma::fill::zeros) ;
  arma::mat fe(nobs, r, arma::fill::zeros) ;

  for (int i = 0; i < nobs; i++) {
    if (k1 > 0) {
      ze.row(i) = Z.row(i) % Alpha.row(id0(i, 0) - 1) ;
    }
    if (k2 > 0) {
      ae.row(i) = A.row(i) % Xi.row(time(i, 0) - 1) ;
    }
    if (r > 0) {
      fe.row(i) = Gamma.row(id0(i, 0) - 1) % F.row(time(i, 0) - 1) ;
    }
  }
  return(arma::join_rows(arma::join_rows(ze, ae), fe)) ;
}

/* -------- sample xi and factor --------- */

// iterated xi sampling for each time period: a T*k matrix 
// [[Rcpp::export]]
arma::mat iterGenXi(arma::mat Xi_old, // a T*k matrix
                    arma::mat Phi, // AR(1) coefficients
                    arma::mat XX,
                    arma::mat XY,
                    int T,
                    int k,
                    double sigma2) {

  arma::mat I(k, k, arma::fill::eye) ;

  // store results
  arma::mat Xi(T, k, arma::fill::zeros) ;

  // covariance matrix
  arma::mat omegat(k, k, arma::fill::zeros) ;
  arma::mat invomegat(k, k, arma::fill::zeros) ;
  
  // mean value matrix
  arma::mat mut(k, 1, arma::fill::zeros) ;
  arma::mat iomu(k, 1, arma::fill::zeros) ;

  // we first set xi0 as a zero vector 
  arma::mat xi0(k, 1, arma::fill::zeros) ;
  arma::mat xi_t1(1, k, arma::fill::zeros) ;
  arma::mat xi_t2(1, k, arma::fill::zeros) ;
  arma::mat xi(k, 1, arma::fill::zeros) ;
  // iomu = Phi * Xi_old.rows(0, k - 1) ;
  // xi0 = arma::mvnrnd(iomu.col(0), arma::inv(Phi*Phi)) ;
  arma::mat mPhi(k, k, arma::fill::zeros) ; // in matrix form 
  for (int j = 0; j < k; j++) {
  	mPhi(j, j) = Phi(j, 0) ;
  }
  

  arma::mat pp = mPhi * mPhi ;

  for (int i = 0; i < T; i++) {

    if (i < (T - 1)) {
      omegat = I + XX.rows(i*k, (i+1)*k - 1) / sigma2 + pp ;
      mut = XY.rows(i*k, (i+1)*k - 1) / sigma2 ;
      
      if (i == 0) {
        // for the first period: only conditional on beta0
        xi_t2 = Xi_old.row(1) ;
        mut = mut + mPhi * xi_t2.t() ;
      }
      else {
        // then: conditional on updated beta_t-1 and beta_t+1 from old
        xi_t1 = Xi_old.row(i-1) ;
        xi_t2 = Xi_old.row(i+1) ;
        mut = mut + mPhi * (xi_t1.t() + xi_t2.t()) ;
      }

    }
    else {
      omegat = I + XX.rows(i*k, (i+1)*k - 1) / sigma2 ;
      xi_t1 = Xi_old.row(i-1) ; 
      mut = XY.rows(i*k, (i+1)*k - 1) / sigma2 + mPhi * xi_t1.t() ;
    }
    
    // sample xi_t
    invomegat = arma::inv(omegat) ;
    iomu = invomegat * mut ;

    xi = arma::mvnrnd(iomu.col(0), invomegat) ;
    Xi.row(i) = xi.t() ; // mean should be a colvec object 

  }
  return(Xi) ;
}


// sample time-varying effects: a T*k matrix 
// [[Rcpp::export]]
arma::mat sampleXi(arma::mat X,
  	               arma::mat y,
  	               arma::mat Xi_old,
  	               arma::mat Phi,
  	               arma::mat timeBreak, // indicator vector
  	               int T,
  	               double sigma2) {
  
  int k = X.n_cols ;

  arma::mat XY = genXY(X, y, timeBreak) ;
  arma::mat XX = genXY(X, X, timeBreak) ;

  arma::mat Xi = iterGenXi(Xi_old, Phi, XX, XY, T, k, sigma2) ;

  return(Xi) ;
}


/* ------------- AR(1) parameter ------------ */

// sample phi for each ar(1) process
// [[Rcpp::export]]
arma::mat samplePhi(arma::mat Xi,  // a T*k matrix
                    arma::mat P0) { // a k*k diagnoal matrix
  
  int T = Xi.n_rows ;
  int k = Xi.n_cols ;
  arma::mat lagXi = Xi ;
  lagXi.shed_row(T-1) ;

  arma::mat nu(1, k, arma::fill::zeros) ;
  lagXi.insert_rows(0, nu) ;

  arma::mat Phi(k, 1, arma::fill::zeros) ;

  // arma::mat muphi(1, 1, arma::fill::zeros) ;
  // arma::mat varphi(1, 1, arma::fill::zerps) ;
  double muphi = 0 ;
  double varphi = 0 ;

  arma::mat xx = lagXi.t() * lagXi ;
  arma::mat xy = lagXi.t() * Xi ;

  for (int i = 0; i < k; i++) {
  	varphi = 1/(xx(i, i) + 1/P0(i, i)) ;
  	muphi = varphi * xy(i, i) ;
  	Phi(i, 0) = R::rnorm(muphi, sqrt(varphi)) ;
  }
  return(Phi) ;
}


/* ------------------- partial fit -------------------- */


// generate random effect
// [[Rcpp::export]]
arma::mat getREfit(arma::mat X,
                   arma::mat eff,
                   arma::mat index) {
  
  int nobs = X.n_rows ;
  arma::mat REfit(nobs, 1, arma::fill::zeros) ;

  for (int i = 0; i < nobs; i++) {
    REfit(i, 0) = arma::accu(X.row(i) % eff.row(index(i, 0) - 1)) ;
  }
  return(REfit) ;
}


// generate factor structure effect
// [[Rcpp::export]]
arma::mat getFactorFit(arma::mat gamma,
                       arma::mat f,
                       arma::mat index1,
                       arma::mat index2) {
  
  int nobs = index1.n_rows ;
  arma::mat Ffit(nobs, 1, arma::fill::zeros) ;

  for (int i = 0; i < nobs; i++) {
    Ffit(i, 0) = arma::accu(gamma.row(index1(i, 0) - 1) % f.row(index2(i, 0) - 1)) ;
  }
  return(Ffit) ;
}

/* -------------- error term ------------- */

// [[Rcpp::export]]
double sampleSigmaE2(arma::mat res, double c0, double d0) {
  
  int nobs = res.n_rows ;
  double c1 = c0 + nobs ;
  double d1 = d0 + arma::accu(pow(res, 2)) ; 
  double sigmaE_2 = 0 ;
  // sampling
  sigmaE_2 = arma::randg<double>( arma::distr_param(c1/2, 2/d1) ) ;
  double sigmaE2 = 1/sigmaE_2 ;

  return(sigmaE2) ; 
}

// [[Rcpp::export]]
double sampleG(double a, double b) {
  

  double c = arma::randg<double>( arma::distr_param(a, b) ) ;


  return(c) ; 
}

// the project
// [[Rcpp::export]]
arma::mat insertB(arma::mat B) {
  arma::mat BB = B;
  int k = B.n_cols ;
  arma::mat nu(1, k, arma::fill::zeros) ;
  BB.insert_rows(0, nu) ;
  return(BB) ;
}

// sample inv-gaussian
// [[Rcpp::export]]
double rrinvgauss(double mu, double lambda) {
 
  double z, y, x, u, m ;
 
  z = R::rnorm(0,1) ;
  y = z * z ;
  x = mu + 0.5 * mu * mu * y / lambda - 0.5 * (mu / lambda) * sqrt(4 * mu * lambda * y + mu * mu * y * y) ;
  u = R::runif(0, 1) ;
  if(u <= mu / (mu + x)) {
    m = x ;
  } else {
    m = mu * mu / x ;
  } 

  return(m) ;
}


/* --------------------- Simulation function ---------------------------- */














