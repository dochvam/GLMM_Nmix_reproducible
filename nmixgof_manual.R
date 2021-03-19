# This code is modified from Knape et al. 2018 (DOI: 10.1111/2041-210X.13062)
# Published with permission from the authors

library(Rcpp)
#u <- runif(F(z-1), F(z))
#rqres = pnorm(u)


# Site-sum residuals
# F(y_si) = Sum{ F_BinSum(y_si|N, p_i1,...,p_iT) * P_i(N) } 

# F_BinSum is the CDF of a sum of independent binomial variables with same N, diff P
# Computed by brute force: 
# F_BinSum(y_si|N, p) = sum{ Binom(k1|N,p1) * ... * Binom(kT|N, pT) }
#   for all combos of ks where k1 + ... + kT <= y_si


Rcpp::sourceCpp("nmixgof_manual.cpp")

# pbinsumRow(c(10,2,2), 10, c(.5,.4,.3))

rqRes = function(y, pFun, size, prob, s = NULL) {
  if (!is.null(s)) {
    # We're doing beta-binomial, which we need to do row-by-row
    b <- y * 0
    a <- y * 0
    
    suppressMessages(
    for (i in 1:nrow(y)) {
      len <- sum(!is.na(y[i,]))
      
      if (any(prob[i,] == 1, na.rm = T)) {
        for (j in 1:len) prob[i,j] <- min(1 - 1e-5, prob[i,j])
      }
      
      b[i, 1:len] <- pFun(q = y[i, 1:len], size = size[i, 1:len], 
                          m = prob[i,1:len], s)
      a[i, 1:len] <- pFun(y[i, 1:len] - 1, size = size[i, 1:len], 
                          m = prob[i,1:len], s = s)
    })
  } else {
    b = pFun(y, size, prob)
    a = pFun(y - 1, size, prob)
  }
  
  if (is.null(dim(y))) {
    res = numeric(length(y)) + NA
  } else {
    res = array(data = NA, dim = dim(y))
  }
  isnna = which(!is.na(y) & !is.na(a) & !is.na(b))
  res[isnna] = rqRes0(a[isnna], b[isnna])
  # if (any(is.infinite(res))) {
  #   warning("Some residuals infinite.")
  # }
  res
}

rqRes0 = function(a, b) {
  stopifnot(length(a) == length(b))
  stats::qnorm(runif(length(a), a, b))
}



# y: mtx
# lam: vector
# p: mtx
rqResS = function(y, lam, p, mixture, K, theta = NULL, s = NULL) {

  cumProb = matrix(0, nrow = nrow(y), ncol = 2)
  dfun = switch(mixture,
                `B-P` = function(N) {stats::dpois(N, lam)},
                `BB-P` = function(N) {stats::dpois(N, lam)},
                `B-NB` = function(N) {stats::dnbinom(N, mu = lam, size = theta)},
                `BB-NB` = function(N) {stats::dnbinom(N, mu = lam, size = theta)},
  )
  
  if (mixture %in% c("B-P", "B-NB")) {
    for (N in 0:K) {
      cumProb = cumProb + kronecker(dfun(N), t(c(1,1))) * 
                pbinsum(y, rep(N, nrow(y)), p)[,2:3]
    }
  } else {
    for (N in 0:K) {
      cumProb = cumProb + kronecker(dfun(N), t(c(1,1))) * 
        pbbinsum(y, rep(N, nrow(y)), p, 1/s)[,2:3] ## TODO: Look at parameterization in the pbbinsuum to see if s = theta (probably doesn't)
    }
  }
  
  res = rqRes0(cumProb[,1], cumProb[,2])
  # if (any(is.infinite(res))) {
  #   warning(paste(sum(is.infinite(res)), " residuals infinite."))
  # }
  res
}

manual_ranef <- function(y, lam, p, K, theta, mixture) {
  R <- length(lam)
  N <- 0:K
  
  post <- array(NA_real_, c(R, length(N), 1))
  colnames(post) <- N
  
  for(i in 1:R) {
    if (grepl(mixture, "NB")) {
      f <- dnbinom(N, mu=lam[i], size=alpha, log=TRUE)
    } else {
      f <- dpois(N, lam[i], log=TRUE)
    }
    
    g <- rep(0, K+1)
    for(j in 1:sum(!is.na(y[i,]))) {
      if(is.na(y[i,j]) | is.na(p[i,j]))
        next
      g <- g + dbinom(y[i,j], N, p[i,j], log=TRUE)
    }
    fudge <- exp(f+g)
    # if (sum(fudge) == 0) {
    #   post[i,,1] <- -Inf
    # } else {
      post[i,,1] <- fudge / sum(fudge)
    # }
  }
  post
}

rqResObs = function(y, lam, p, mixture, K, theta = NULL, s = NULL) {
  rN = integer(nrow(y)) + NA
  
  sample_block_na <- function(x, size, replace, prob) {
    if (any(is.na(prob))) {
      return(rep(NaN, size))
    }
    sample(x, size, replace, prob)
  }
  
  rN = apply(manual_ranef(y, lam, p, K, theta, mixture)[,,1], 1, 
             sample_block_na, x = K + 1, size = 1, replace = FALSE) - 1
  
  sites_to_flag <- which(is.nan(rN))
  
  if (grepl("BB", mixture)) {
    res = rqRes(y, pFun = rmutil::pbetabinom, s = s,
                size = kronecker(rN, t(rep(1, ncol(p)))), prob=p)
  } else {
    res = rqRes(y, pFun = pbinom,
                size = kronecker(rN, t(rep(1, ncol(p)))), prob=p)
  }
  
  
  for (i in sites_to_flag) {
    len <- sum(!is.na(y[i,]))
    res[i, 1:len] <- NaN
  }
  
  if (any(is.infinite(res))) {
    #browser()
    warning(paste(sum(is.infinite(res)), " residuals infinite."))
  }
  res
}

