# 06_post_tests.R
# Author: Benjamin R. Goldstein
# Date: 2/1/2021

### Script 3: Main and helper functions for model fitting.
# This file contains functions to batch-run post tests on the model,
# including goodness-of-fit checking, checking residuals for autocorrelation,
# stability checks, etc.


library(DHARMa)
library(unmarked)
library(nmixgof)
library(tidyverse)
library(parallel)
library(ape)
source("03_fit_one_ssr_final.R")
source("nmixgof_manual.R")


#### Helper functions ####
refit_as_nimble <- function(nmix_res_list, K = NULL, dat = NULL) {
  
  if (is.null(K)) {
    K <- nmix_res_list$K
  }
  mixture <- nmix_res_list$mixture
  best_row <- nmix_res_list$best_row
  
  if (is.null(dat)) {
    dat_df <- make_dat_df(nmix_res_list$species, nmix_res_list$subregion)
  } else {
    dat_df <- dat
  }
  
  accept_det_names <- clean_vec_names(c("duration", "num_observers"), dat_df)
  accept_abd_names <- c("elevation")#clean_vec_names(c("elevation", "precip", "tmax"), dat_df)
  
  modellist <- make_nmix_model(dat_df, names(best_row)[-length(best_row)], 
                               mixture, K[1])
  nmix_model <- modellist$nmix_model
  Cnmix_model <- modellist$Cnmix_model
  
  all_zero_nodes <- c(
    nmix_model$expandNodeNames("abund_coeffs"), 
    nmix_model$expandNodeNames("p_coeffs"),
    if (mixture %in% c("B-NB", "BB-NB")) "ltheta",
    if (mixture %in% c("BB-P", "BB-NB")) "log_s"
  )
  
  zero_assignments <- nimble::values(Cnmix_model, all_zero_nodes)
  
  target_nodes <- row_to_nodes(this_row = nmix_res_list$best_row, 
                               mixture = nmix_res_list$mixture, 
                               dat_df = dat_df, 
                               accept_det_names = accept_det_names, 
                               accept_abd_names = accept_abd_names)
  refit <- list()
  for (i in 1:length(K)) {
    nmix_model$K <- K[i]
    Cnmix_model$K <- K[i]
    refit[[i]] <- Cnmix_model_optim(Cnmix_model, nmix_model,
                                    target_nodes, all_zero_nodes,
                                    zero_assignments)
  }
  
  return(refit)
}

refit_as_unmarked <- function(nmix_res) {
  
  dat_df <- make_dat_df(nmix_res$species, nmix_res$subregion)
  mixture <- nmix_res$mixture
  
  if (!mixture %in% c("B-P", "B-NB")) {
    stop("unmarked only supports B-P and B-NB mixtures.")
  }
  umMix <- substr(mixture, 3, nchar(mixture))
  
  accept_det_names <- clean_vec_names(c("duration", "num_observers"), dat_df)
  accept_abd_names <- c("elevation") #clean_vec_names(c("elevation", "precip", "tmax"), dat_df)
  
  f <- unmarked_formula_from_row(nmix_res$best_row, 
                                 accept_det_names = accept_det_names,
                                 accept_abd_names = accept_abd_names)
  
  umfpc <- unmarkedFramePCount(
    y = get_var_wide(dat_df, "total_count"), 
    obsCovs = list(
      duration = get_var_wide(dat_df, "duration"), 
      num_observers = get_var_wide(dat_df, "num_observers"), 
      tod = get_var_wide(dat_df, "tod"), 
      tod_sq = get_var_wide(dat_df, "tod_sq"), 
      yday = get_var_wide(dat_df, "yday"), 
      yday_sq = get_var_wide(dat_df, "yday_sq") 
      # distance = get_var_wide(dat_df, "distance") 
    ),
    siteCovs = data.frame(
      elevation = get_var_wide(dat_df, "elevation")[,1],
      # precip = get_var_wide(dat_df, "precip")[,1], 
      # tmax = get_var_wide(dat_df, "tmax")[,1], 
      pct_water = get_var_wide(dat_df, "pct_water")[,1], 
      pct_ag = get_var_wide(dat_df, "pct_ag")[,1], 
      pct_tree = get_var_wide(dat_df, "pct_tree")[,1], 
      pct_veg = get_var_wide(dat_df, "pct_veg")[,1]
    )
  )
  
  umfit <- pcount(f, umfpc, K = nmix_res$K, mixture = umMix)
  
  return(list(
    refit_AIC = umfit@AIC,
    original_AIC = nmix_res$AIC,
    coeffs = umfit@estimates
  ))
}



rqresid_from_nmixfit <- function(nmix_res, type = "Site-Sum") {
  
  if (!type %in% c("Site-Sum", "Observation", "Marginal")) {
    stop("'type' must be one of: 'Site-Sum', 'Observation', 'Marginal'")
  }
  
  mixture <- nmix_res$mixture
  
  dat_df <- make_dat_df(nmix_res$species, nmix_res$subregion)
  
  coeffs <- nmix_res$coefficients
  
  accept_det_names <- clean_vec_names(c("duration", "num_observers"), dat_df)
  accept_abd_names <- c("elevation") #clean_vec_names(c("elevation", "precip", "tmax"), dat_df)
  
  det_pars <- coeffs %>% 
    filter(param %in% accept_det_names |
             grepl("det", param)) %>% 
    mutate(param = gsub("-det", "", param))
  abd_pars <- coeffs %>% 
    filter(param %in% accept_abd_names |
             grepl("abd", param)) %>% 
    mutate(param = gsub("-abd", "", param))
  
  det_pars$param[det_pars$param == "tod*tod"] <- "tod_sq"
  det_pars$param[det_pars$param == "yday*yday"] <- "yday_sq"
  
  det_formula <- as.formula(paste0("~", paste(det_pars$param, collapse = "+")))
  if (nrow(abd_pars) > 0) {
    abd_formula <- as.formula(paste0("~", paste(abd_pars$param, collapse = "+")))
  } else {
    abd_formula <- ~1
  }
  
  det_mm <- model.matrix(det_formula, dat_df)
  abd_mm <- model.matrix(abd_formula, dat_df)
  
  logit_p_intercept <- logit(
    exp((coeffs$est[coeffs$param == "log(p * lambda)"] +
           coeffs$est[coeffs$param == "log(p / lambda)"]) / 2)
  )
  log_lam_intercept <- log(
    exp((coeffs$est[coeffs$param == "log(p * lambda)"] -
           coeffs$est[coeffs$param == "log(p / lambda)"]) / 2)
  )
  
  dat_df$p <- det_mm %*% 
    c(logit_p_intercept,
      det_pars$est) %>% 
    expit()
  dat_df$p <- unlist(lapply(dat_df$p, function(x) min(x, 1 - 1e-5)))
  
  dat_df$lambda <- abd_mm %*% 
    c(log_lam_intercept,
      abd_pars$est) %>% 
    exp() 
  
  y_mtx <- get_var_wide(dat_df, "total_count")
  p_mtx <- get_var_wide(dat_df, "p")
  lambda_mtx <- get_var_wide(dat_df, "lambda")
  
  if (mixture %in% c("B-NB", "BB-NB")) {
    theta <- exp(coeffs$est[coeffs$param == "ltheta"])
  } else {
    theta <- NULL
  }
  
  if (mixture %in% c("BB-P", "BB-NB")) {
    s <- exp(coeffs$est[coeffs$param == "log_s"])
  } else {
    s <- NULL
  }
  
  
  if (type == "Site-Sum") {
    return(
      rqResS(y = y_mtx, lam = lambda_mtx[,1], p = p_mtx, 
             mixture = nmix_res$mixture, K = nmix_res$K,
             theta = theta, s = s)
    )
  } else if (type == "Observation") {
    return(
      rqResObs(y = y_mtx, lam = lambda_mtx[,1], p = p_mtx, 
               mixture = nmix_res$mixture, K = nmix_res$K,
               theta = theta, s = s)
    )
  } else if (type == "Marginal") {
    return(
      rqResMarginal(y = y_mtx, lam_mtx = lambda_mtx, p = p_mtx, 
                    mixture = nmix_res$mixture, K = nmix_res$K,
                    theta = theta, s = s)
    )
  }
}


modifiedSimResiduals <- 
  function (fittedModel, n = 250, refit = F, integerResponse = NULL, 
            plot = F, seed = 123, method = c("PIT", "traditional"), 
            ...) {
    if (n < 2) 
      stop("error in DHARMa::simulateResiduals: n > 1 is required to calculate scaled residuals")
    DHARMa:::checkModel(fittedModel)
    match.arg(method)
    randomState <- getRandomState(seed)
    on.exit({
      randomState$restoreCurrent()
    })
    ptm <- proc.time()
    out = list()
    family = family(fittedModel)
    out$fittedModel = fittedModel
    out$modelClass = class(fittedModel)[1]
    out$nObs = nobs(fittedModel)
    out$nSim = n
    out$refit = refit
    out$observedResponse = getObservedResponse(fittedModel)
    if (is.null(integerResponse)) {
      if (family$family %in% c("binomial", "poisson", "quasibinomial", 
                               "quasipoisson", "Negative Binom", "nbinom2", "nbinom1", 
                               "genpois", "compois", "truncated_poisson", "truncated_nbinom2", 
                               "truncated_nbinom1", "betabinomial", "Poisson", 
                               "Tpoisson", "COMPoisson", "negbin", "Tnegbin") | 
          grepl("Negative Binomial", family$family)) 
        integerResponse = TRUE
      else integerResponse = FALSE
    }
    out$integerResponse = integerResponse
    out$problems = list()
    if (out$modelClass %in% c("HLfit")) {
      out$fittedPredictedResponse = predict(fittedModel, type = "response", 
                                            re.form = ~0)[, 1L]
    }
    else {
      out$fittedPredictedResponse = predict(fittedModel, type = "response", 
                                            re.form = ~0)
    }
    out$fittedFixedEffects = getFixedEffects(fittedModel)
    out$fittedResiduals = residuals(fittedModel, type = "response")
    if (refit == FALSE) {
      out$simulatedResponse = getSimulations(fittedModel, 
                                             nsim = n, type = "normal", ...)
      out$simulatedResponse[is.nan(out$simulatedResponse)] <- 0
      DHARMa:::checkSimulations(out$simulatedResponse, out$nObs, out$nSim)
      out$scaledResiduals = getQuantile(simulations = out$simulatedResponse, 
                                        observed = out$observedResponse, integerResponse = integerResponse, 
                                        method = method)
    }
    else {
      out$refittedPredictedResponse <- matrix(nrow = out$nObs, 
                                              ncol = n)
      out$refittedFixedEffects <- matrix(nrow = length(out$fittedFixedEffects), 
                                         ncol = n)
      out$refittedResiduals = matrix(nrow = out$nObs, ncol = n)
      out$refittedPearsonResiduals = matrix(nrow = out$nObs, 
                                            ncol = n)
      out$simulatedResponse = getSimulations(fittedModel, 
                                             nsim = n, type = "refit", ...)
      for (i in 1:n) {
        simObserved = out$simulatedResponse[[i]]
        try({
          refittedModel = getRefit(fittedModel, simObserved)
          out$refittedPredictedResponse[, i] = predict(refittedModel, 
                                                       type = "response")
          out$refittedFixedEffects[, i] = getFixedEffects(refittedModel)
          out$refittedResiduals[, i] = residuals(refittedModel, 
                                                 type = "response")
          out$refittedPearsonResiduals[, i] = residuals(refittedModel, 
                                                        type = "pearson")
        }, silent = TRUE)
      }
      if (anyNA(out$refittedResiduals)) 
        warning("DHARMa::simulateResiduals warning: on refit = TRUE, at least one of the refitted models produced an error. Inspect the refitted model values. Results may not be reliable.")
      dup = sum(duplicated(out$refittedFixedEffects, MARGIN = 2))
      if (dup > 0) {
        if (dup < n/3) {
          warning(paste("There were", dup, "of", n, "duplicate parameter estimates in the refitted models. This may hint towards a problem with optimizer convergence in the fitted models. Results may not be reliable. The suggested action is to not use the refitting procedure, and diagnose with tools available for the normal (not refitted) simulated residuals. If you absolutely require the refitting procedure, try changing tolerance / iterations in the optimizer settings."))
        }
        else {
          warning(paste("There were", dup, "of", n, "duplicate parameter estimates in the refitted models. This may hint towards a problem with optimizer convergence in the fitted models. Results are likely not reliable. The suggested action is to not use the refitting procedure, and diagnose with tools available for the normal (not refitted) simulated residuals. If you absolutely require the refitting procedure, try changing tolerance / iterations in the optimizer settings."))
          out$problems[[length(out$problems) + 1]] = "error in refit"
        }
      }
      out$scaledResiduals = getQuantile(simulations = out$refittedResiduals, 
                                        observed = out$fittedResiduals, integerResponse = integerResponse, 
                                        method = method)
    }
    out$time = proc.time() - ptm
    out$randomState = randomState
    class(out) = "DHARMa"
    if (plot == TRUE) 
      plot(out)
    return(out)
  }


#### Function for doing gof given list of 4 model results ####
# four_res_list should have a result_list object for each of the four models in
# alphabetical order by filename (GLMM_Nbin, GLMM_Pois, Nmix_Nbin, Nmix_Pois)
gof_by_ssr <- function(ssr_str, onemodel_path) {
  
  glmm_nbin_res <- readRDS(paste0(onemodel_path, "GLMM_Nbin", ssr_str, ".RDS"))
  glmm_pois_res <- readRDS(paste0(onemodel_path, "GLMM_Pois", ssr_str, ".RDS"))
  
  nmix_BBNB_res <- readRDS(paste0(onemodel_path, "Nmix_BBNB", ssr_str, ".RDS"))
  nmix_BBP_res  <- readRDS(paste0(onemodel_path, "Nmix_BBP", ssr_str, ".RDS"))
  nmix_BNB_res  <- readRDS(paste0(onemodel_path, "Nmix_BNB", ssr_str, ".RDS"))
  nmix_BP_res   <- readRDS(paste0(onemodel_path, "Nmix_BP", ssr_str, ".RDS"))
  
  six_res_list <- list(glmm_nbin_res,
                       glmm_pois_res,
                       nmix_BBNB_res,
                       nmix_BBP_res ,
                       nmix_BNB_res ,
                       nmix_BP_res)
  simResidP <- NULL
  simResidNB <- NULL
  
  tryCatch({
    simResidP  <- simulateResiduals(fittedModel = glmm_pois_res$fit, n = 1000)
    simResidNB <- simulateResiduals(fittedModel = glmm_nbin_res$fit, n = 1000)
  }, error = function(e) {})
  
  tryCatch({
    if (is.null(simResidP)) {
      simResidP  <- simulateResiduals(fittedModel = glmm_pois_res$fit, n = 10000)
    }
    if (is.null(simResidNB)) {
      simResidNB <- simulateResiduals(fittedModel = glmm_nbin_res$fit, n = 10000)
    }
  }, error = function(e) {})
  
  tryCatch({
    if (is.null(simResidP)) {
      simResidP  <- modifiedSimResiduals(fittedModel = glmm_pois_res$fit, n = 10000)
    }
    if (is.null(simResidNB)) {
      simResidNB <- modifiedSimResiduals(fittedModel = glmm_nbin_res$fit, n = 10000)
    }
  }, error = function(e) {})
  
  gof <- data.frame(
    test = c(rep(c("Uniformity", "Dispersion"), 2),
             rep("Uniformity", 8)),
    resid_type = rep(c("DHARMa", "Site-Sum RQ", "Marginal RQ"), each = 4),
    fit_model = c(rep(c("GLMM_P", "GLMM_NB"), each = 2),
                  rep(c("Nmix_BBNB", "Nmix_BBP", "Nmix_BNB", "Nmix_BP"), 2)),
    species = glmm_nbin_res$species,
    subregion = glmm_nbin_res$subregion,
    pvalue = NA, stat = NA,
    chosenAIConly = c("GLMM_NB", "GLMM_P", "Nmix_BBNB",
                      "Nmix_BBP", "Nmix_BNB", "Nmix_BP")[
                        which.min(unlist(lapply(six_res_list, function(x) x$AIC)))
                      ]
  )
  
  # GLMM gof
  result_unif_P <- testUniformity(simResidP, plot = F)
  result_disp_P <- testDispersion(simResidP, plot = F)
  result_unif_NB <- testUniformity(simResidNB, plot = F)
  result_disp_NB <- testDispersion(simResidNB, plot = F)
  
  gof$pvalue[1] <- result_unif_P$p.value
  gof$pvalue[2] <- result_disp_P$p.value
  gof$pvalue[3] <- result_unif_NB$p.value
  gof$pvalue[4] <- result_disp_NB$p.value
  
  gof$stat[1] <- result_unif_P$statistic
  gof$stat[2] <- result_disp_P$statistic
  gof$stat[3] <- result_unif_NB$statistic
  gof$stat[4] <- result_disp_NB$statistic
  
  # get RQ resids for each N-mixture fit
  if (!file.exists(paste0("output/residuals/rqres_ss", ssr_str, ".RDS"))) {
    nmix_BP_resids   <- rqresid_from_nmixfit(nmix_res = nmix_BP_res, type = "Site-Sum")
    nmix_BNB_resids  <- rqresid_from_nmixfit(nmix_res = nmix_BNB_res, type = "Site-Sum")
    nmix_BBP_resids  <- rqresid_from_nmixfit(nmix_res = nmix_BBP_res, type = "Site-Sum")
    nmix_BBNB_resids <- rqresid_from_nmixfit(nmix_res = nmix_BBNB_res, type = "Site-Sum")
    all_nmix_resids_ss <- list(nmix_BP_resids,
                               nmix_BNB_resids,
                               nmix_BBP_resids,
                               nmix_BBNB_resids)
    saveRDS(all_nmix_resids_ss, paste0("output/residuals/rqres_ss", ssr_str, ".RDS"))
  } else {
    all_nmix_resids_ss <- readRDS(paste0("output/residuals/rqres_ss", ssr_str, ".RDS"))
  }
  if (!file.exists(paste0("output/residuals/rqres_mar", ssr_str, ".RDS"))) {
    nmix_BP_resids   <- rqresid_from_nmixfit(nmix_res = nmix_BP_res, type = "Marginal")
    nmix_BNB_resids  <- rqresid_from_nmixfit(nmix_res = nmix_BNB_res, type = "Marginal")
    nmix_BBP_resids  <- rqresid_from_nmixfit(nmix_res = nmix_BBP_res, type = "Marginal")
    nmix_BBNB_resids <- rqresid_from_nmixfit(nmix_res = nmix_BBNB_res, type = "Marginal")
    all_nmix_resids_mar <- list(nmix_BP_resids,
                                nmix_BNB_resids,
                                nmix_BBP_resids,
                                nmix_BBNB_resids)
    saveRDS(all_nmix_resids_mar, paste0("output/residuals/rqres_mar", ssr_str, ".RDS"))
  } else {
    all_nmix_resids_mar <- readRDS(paste0("output/residuals/rqres_mar", ssr_str, ".RDS"))
  }
  
  
  ### Site-sum first
  nmix_BP_resids   <- all_nmix_resids_ss[[1]]
  nmix_BNB_resids  <- all_nmix_resids_ss[[2]]
  nmix_BBP_resids  <- all_nmix_resids_ss[[3]]
  nmix_BBNB_resids <- all_nmix_resids_ss[[4]]
  
  # N-mixture gof
  gof$pvalue[8] <- ks.test(nmix_BP_resids  , pnorm)$p.value
  gof$pvalue[7] <- ks.test(nmix_BNB_resids , pnorm)$p.value
  gof$pvalue[6] <- ks.test(nmix_BBP_resids , pnorm)$p.value
  gof$pvalue[5] <- ks.test(nmix_BBNB_resids, pnorm)$p.value  
  
  gof$stat[8] <- ks.test(nmix_BP_resids  , pnorm)$statistic
  gof$stat[7] <- ks.test(nmix_BNB_resids , pnorm)$statistic
  gof$stat[6] <- ks.test(nmix_BBP_resids , pnorm)$statistic
  gof$stat[5] <- ks.test(nmix_BBNB_resids, pnorm)$statistic
  
  
  
  ### Now do marginal
  nmix_BP_resids   <- numeric(nrow(all_nmix_resids_mar[[1]]))
  nmix_BNB_resids  <- numeric(nrow(all_nmix_resids_mar[[1]]))
  nmix_BBP_resids  <- numeric(nrow(all_nmix_resids_mar[[1]]))
  nmix_BBNB_resids <- numeric(nrow(all_nmix_resids_mar[[1]]))
  
  for (i in 1:length(nmix_BP_resids)) {
    ind <- sample(1:sum(!is.na(all_nmix_resids_mar[[1]][i,])), 1)
    nmix_BP_resids[i]   <- all_nmix_resids_mar[[1]][i,ind]
    nmix_BNB_resids[i]  <- all_nmix_resids_mar[[2]][i,ind]
    nmix_BBP_resids[i]  <- all_nmix_resids_mar[[3]][i,ind]
    nmix_BBNB_resids[i] <- all_nmix_resids_mar[[4]][i,ind]
  }
  
  # N-mixture gof
  gof$pvalue[12] <- ks.test(nmix_BP_resids  , pnorm)$p.value
  gof$pvalue[11] <- ks.test(nmix_BNB_resids , pnorm)$p.value
  gof$pvalue[10] <- ks.test(nmix_BBP_resids , pnorm)$p.value
  gof$pvalue[9]  <- ks.test(nmix_BBNB_resids, pnorm)$p.value  
  
  gof$stat[12] <- ks.test(nmix_BP_resids  , pnorm)$statistic
  gof$stat[11] <- ks.test(nmix_BNB_resids , pnorm)$statistic
  gof$stat[10] <- ks.test(nmix_BBP_resids , pnorm)$statistic
  gof$stat[9]  <- ks.test(nmix_BBNB_resids, pnorm)$statistic
  
  return(gof)
}

stability_check <- function(ssr_str, moddist, onemodel_path) {
  file <- list.files(onemodel_path, 
                      pattern = paste0(moddist, "_*", ssr_str), full.names = TRUE)

  nmix_res <- readRDS(file)
  
  test_refits <- refit_as_nimble(nmix_res, K = c(
    nmix_res$K + 2000, nmix_res$K + 4000
  ))

  AICs <- lapply(test_refits, function(x) AIC_optim(x)) %>% unlist()
  
  coeffs_df <- nmix_res$coefficients %>% 
    select(-se) %>% 
    bind_rows(data.frame(
      param = "AIC", est = nmix_res$AIC
    )) %>% 
    mutate(K = nmix_res$K)
  
  refit_coeffs_df_list <- lapply(test_refits, function(x) {
    data.frame(
      param = coeffs_df$param,
      est = c(log(expit(x$par[1])) + x$par[2], 
              log(expit(x$par[1])) - x$par[2],
              x$par[3:length(x$par)],
              AIC_optim(x))
    )
  })
  refit_coeffs_df_list[[1]]$K <- nmix_res$K + 2000
  refit_coeffs_df_list[[2]]$K <- nmix_res$K + 4000
  refit_coeffs_df_list[[3]] <- coeffs_df
  
  result_df <- do.call(rbind, refit_coeffs_df_list)
  
  return(result_df)
}


check_autocorr <- function(species, subregion, resid_path, onemodel_path,
                           glmm_only = FALSE) {
  
  ssr_str <- paste0("_", subregion, "_", species)
  dat_df <- make_dat_df(species, subregion, include_ll = TRUE) %>% 
    mutate(
      lon = scale(gx),
      lat = scale(gy)
    )
  
  glmm_mods <- lapply(list.files(
    path = onemodel_path, full.names = TRUE, pattern = paste0("GLMM_.*", ssr_str)
  ), readRDS)
  
  if (!glmm_only) {
    nmix_ss_resids_file <- list.files(resid_path, pattern = paste0("rqres_ss", ssr_str),
                                   full.names = TRUE)
    nmix_ss_resids_list <- readRDS(nmix_ss_resids_file)
    nmix_mar_resids_file <- list.files(resid_path, pattern = paste0("rqres_mar", ssr_str),
                                      full.names = TRUE)
    nmix_mar_resids_list <- readRDS(nmix_mar_resids_file)
    
    # Randomly subsample the marginal resids
    nmix_BP_resids   <- numeric(nrow(nmix_mar_resids_list[[1]]))
    nmix_BNB_resids  <- numeric(nrow(nmix_mar_resids_list[[1]]))
    nmix_BBP_resids  <- numeric(nrow(nmix_mar_resids_list[[1]]))
    nmix_BBNB_resids <- numeric(nrow(nmix_mar_resids_list[[1]]))
    
    for (i in 1:length(nmix_BP_resids)) {
      ind <- sample(1:sum(!is.na(nmix_mar_resids_list[[1]][i,])), 1)
      nmix_BP_resids[i]   <- nmix_mar_resids_list[[1]][i,ind]
      nmix_BNB_resids[i]  <- nmix_mar_resids_list[[2]][i,ind]
      nmix_BBP_resids[i]  <- nmix_mar_resids_list[[3]][i,ind]
      nmix_BBNB_resids[i] <- nmix_mar_resids_list[[4]][i,ind]
    }
    
    sites <- dat_df %>% 
      distinct(locID, gx, gy) %>% 
      mutate(
        ss_BP_resids    = nmix_ss_resids_list[[1]],
        ss_BNB_resids   = nmix_ss_resids_list[[2]],
        ss_BBP_resids   = nmix_ss_resids_list[[3]],
        ss_BBNB_resids  = nmix_ss_resids_list[[4]],
        mar_BP_resids   = nmix_BP_resids,
        mar_BNB_resids  = nmix_BNB_resids,
        mar_BBP_resids  = nmix_BBP_resids,
        mar_BBNB_resids = nmix_BBNB_resids
      )
  } else {
    sites <- dat_df %>% 
      distinct(locID, gx, gy)
  }
  
  pois_moran <- list(p.value = NA)
  nbin_moran <- list(p.value = NA)
  tryCatch({
    glmm_pois_resids <- simulateResiduals(glmm_mods[[2]]$fit, n = 500)
    glmm_nbin_resids <- simulateResiduals(glmm_mods[[1]]$fit, n = 500)
    
    glmm_p_agg_resids <- recalculateResiduals(glmm_pois_resids, group = dat_df$locID)
    glmm_nb_agg_resids <- recalculateResiduals(glmm_nbin_resids, group = dat_df$locID)
    
    pois_moran <- testSpatialAutocorrelation(glmm_p_agg_resids, 
                                             x = sites$gx, y = sites$gy, plot = F)
    nbin_moran <- testSpatialAutocorrelation(glmm_nb_agg_resids, 
                                             x = sites$gx, y = sites$gy, plot = F)
  }, error = function(err) {})
  # glmm_dists <- as.matrix(dist(dat_df[, c("lon", "lat")]))
  # glmm_weights <- 1 / (glmm_dists)
  # glmm_weights[is.infinite(glmm_weights)] <- 0
  # diag(glmm_weights) <- 0
  # # if (exclude_samesite) glmm_weights[glmm_dists == 0] <- 0
  # pois_moran <- NA
  # nbin_moran <- NA
  # tryCatch({
  #   pois_moran <- Moran.I(dat_df$pois_resids, weight = glmm_weights)
  #   nbin_moran <- Moran.I(dat_df$nbin_resids, weight = glmm_weights)
  # }, error = function(err) {})
  # if (is.na(pois_moran[[1]])) pois_moran <- list(p.value = NA)
  # if (is.na(nbin_moran[[1]])) nbin_moran <- list(p.value = NA)
  
  if (!glmm_only) {
    site_dists <- as.matrix(dist(sites[, c("gx", "gy")], ))
    site_weights <- 1 / (site_dists)
    site_weights[is.infinite(site_weights)] <- 0
    diag(site_weights) <- 0
    # if (exclude_samesite) site_weights[site_dists == 0] <- 0
    
    bp_sites <- which(is.finite(sites$ss_BP_resids))
    bp_moran <- Moran.I(sites$ss_BP_resids[bp_sites], 
                        weight = site_weights[bp_sites, bp_sites], 
                        na.rm = TRUE)
  
    bnb_sites <- which(is.finite(sites$ss_BNB_resids))
    bnb_moran <- Moran.I(sites$ss_BNB_resids[bnb_sites], 
                        weight = site_weights[bnb_sites, bnb_sites], 
                        na.rm = TRUE)
    
    bbp_sites <- which(is.finite(sites$ss_BBP_resids))
    bbp_moran <- Moran.I(sites$ss_BBP_resids[bbp_sites], 
                         weight = site_weights[bbp_sites, bbp_sites], 
                         na.rm = TRUE)
    bbnb_sites <- which(is.finite(sites$ss_BBNB_resids))
    bbnb_moran <- Moran.I(sites$ss_BBNB_resids[bbnb_sites], 
                         weight = site_weights[bbnb_sites, bbnb_sites], 
                         na.rm = TRUE)
    
    
    # Same for marginal resids
    mar_bp_sites <- which(is.finite(sites$mar_BP_resids))
    mar_bp_moran <- Moran.I(sites$mar_BP_resids[mar_bp_sites], 
                        weight = site_weights[mar_bp_sites, mar_bp_sites], 
                        na.rm = TRUE)
    
    mar_bnb_sites <- which(is.finite(sites$mar_BNB_resids))
    mar_bnb_moran <- Moran.I(sites$mar_BNB_resids[mar_bnb_sites], 
                         weight = site_weights[mar_bnb_sites, mar_bnb_sites], 
                         na.rm = TRUE)
    
    mar_bbp_sites <- which(is.finite(sites$mar_BBP_resids))
    mar_bbp_moran <- Moran.I(sites$mar_BBP_resids[mar_bbp_sites], 
                         weight = site_weights[mar_bbp_sites, mar_bbp_sites], 
                         na.rm = TRUE)
    mar_bbnb_sites <- which(is.finite(sites$mar_BBNB_resids))
    mar_bbnb_moran <- Moran.I(sites$mar_BBNB_resids[mar_bbnb_sites], 
                          weight = site_weights[mar_bbnb_sites, mar_bbnb_sites], 
                          na.rm = TRUE)
  }
  
  # sp::coordinates(sites) = ~lon + lat
  # variog_resid <- variogram(BBNB_resids~1, data = sites)
  # saveRDS(variog_resid, file = "variogram_residuals.Rds")
  # var_resid_fit <- fit.variogram(variog_resid, vgm("Exp"), fit.kappa = TRUE)
  # plot(variogramLine(var_resid_fit, maxdist = 5), type = "l")
  if (glmm_only) {
    findings <- data.frame(
      species = species,
      subregion = subregion,
      moddist = c("GLMM_Pois", "GLMM_Nbin"),
      resid_type = "DHARMa",
      p.value = c(pois_moran$p.value, nbin_moran$p.value)
    )
  } else {
    findings <- data.frame(
      species = species,
      subregion = subregion,
      moddist = c("Nmix_BP", "Nmix_BNB", "Nmix_BBP", "Nmix_BBNB",
                  "Nmix_BP", "Nmix_BNB", "Nmix_BBP", "Nmix_BBNB",
                  "GLMM_Pois", "GLMM_Nbin"),
      resid_type = c(rep(c("Site-Sum RQ", "Marginal RQ"), each = 4),
                     "DHARMa", "DHARMa"),
      p.value = c(bp_moran$p.value, bnb_moran$p.value, 
                  bbp_moran$p.value, bbnb_moran$p.value,
                  mar_bp_moran$p.value, mar_bnb_moran$p.value, 
                  mar_bbp_moran$p.value, mar_bbnb_moran$p.value,
                  pois_moran$p.value, nbin_moran$p.value)
    )
  }
  
  return(findings)
}



ca_obs <- read_csv("intermediate/CA_obs_processed.csv")
ca_checklists <- read_csv("intermediate/CA_checklists_w_covariates.csv")
