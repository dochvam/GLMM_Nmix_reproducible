# model_code.R
# Author: Benjamin R. Goldstein

# This script contains NIMBLE code specifying four N-mixture models,
# representing the four considered (B-P, B-NB, BB-P, BB-NB). 
# Only the first model is commented thoroughly.

BP_nmix_code_ebd <- nimbleCode({
  
  # Logit-linked observation probabilities for each observation. 
  # Observations are provided as a vector, with a second vector indicating
  # where the obs for each site begin and end
  for (o in 1:nobs) {
    logit(p[o]) <- p_coeffs[1] +
      sum(p_coeffs[2:npcov] * p_dat[o, 1:(npcov-1)])
  }
  
  # Loop over sites with at least one observation
  for (s in 1:(nsite - nSite1Obs)) {
    # Log-linked mean abundance
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    # Data are distributed according to the custom distribution
    y[index_start[s]:index_end[s]] ~ 
      dNmixture_v(lambda = lambda[s],
                     prob = p[index_start[s]:index_end[s]],
                     Nmin = 0,
                     Nmax = K,
                     len = length(index_start[s]:index_end[s]))
  }
  
  # Sites with only one observation are processed separately
  for (s in (nsite - nSite1Obs + 1):nsite) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]] ~ dpois(lambda = lambda[s] * p[index_start[s]])
    
  }
})

BNB_nmix_code_ebd <- nimbleCode({
  
  for (o in 1:nobs) {
    logit(p[o]) <- p_coeffs[1] +
      sum(p_coeffs[2:npcov] * p_dat[o, 1:(npcov-1)])
  }
  
  theta <- exp(ltheta)
  
  for (s in 1:(nsite - nSite1Obs)) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]:index_end[s]] ~
      dNmixture_BNB_v(lambda = lambda[s],
                      theta = theta,
                      prob = p[index_start[s]:index_end[s]],
                      Nmin = 0,
                      Nmax = K,
                      len = length(index_start[s]:index_end[s]))
  }
  
  for (s in (nsite - nSite1Obs + 1):nsite) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]] ~ dNmixture_BNB_oneObs(lambda = lambda[s],
                                             prob = p[index_start[s]],
                                             theta = theta,
                                             Nmin = 0,
                                             Nmax = K,
                                             len = 1)
  }
})

BBP_nmix_code_ebd <- nimbleCode({
  
  for (o in 1:nobs) {
    logit(p[o]) <- p_coeffs[1] +
      sum(p_coeffs[2:npcov] * p_dat[o, 1:(npcov-1)])
  }
  
  s_bb <- exp(log_s)
  
  for (s in 1:(nsite - nSite1Obs)) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]:index_end[s]] ~
      dNmixture_BBP_v(lambda = lambda[s],
                      prob = p[index_start[s]:index_end[s]],
                      s = s_bb,
                      Nmin = 0,
                      Nmax = K,
                      len = length(index_start[s]:index_end[s]))
  }
  
  for (s in (nsite - nSite1Obs + 1):nsite) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]] ~ dNmixture_BBP_oneObs(lambda = lambda[s],
                                             prob = p[index_start[s]],
                                             s = s_bb,
                                             Nmin = 0,
                                             Nmax = K,
                                             len = 1)
  }
})


BBNB_nmix_code_ebd <- nimbleCode({
  
  for (o in 1:nobs) {
    logit(p[o]) <- p_coeffs[1] +
      sum(p_coeffs[2:npcov] * p_dat[o, 1:(npcov-1)])  
  }
  
  s_bb <- exp(log_s)
  theta <- exp(ltheta)  
  
  for (s in 1:(nsite - nSite1Obs)) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]:index_end[s]] ~
      dNmixture_BBNB_v(lambda = lambda[s],
                       theta = theta,
                       prob = p[index_start[s]:index_end[s]],
                       s = s_bb,
                       Nmin = 0,
                       Nmax = K,
                       len = length(index_start[s]:index_end[s]))
  }
  
  for (s in (nsite - nSite1Obs + 1):nsite) {
    log(lambda[s]) <- 
      abund_coeffs[1] + # intercept
      sum(abund_coeffs[2:nlcov] * # data columns
            abund_dat[index_start[s], 1:(nlcov-1)])
    
    y[index_start[s]] ~ dNmixture_BBNB_oneObs(lambda = lambda[s],
                                              theta = theta,
                                              prob = p[index_start[s]],
                                              s = s_bb,
                                              Nmin = 0,
                                              Nmax = K,
                                              len = 1)
  }
})
