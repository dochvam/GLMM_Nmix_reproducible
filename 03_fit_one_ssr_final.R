### Script 3: Main and helper functions for model fitting.

library(furrr)
library(glmmTMB)
library(lme4)
library(tidyverse)
library(nimble)
library(unmarked)
library(parallel)
library(numDeriv)

source("dNmixture_dists.R")
source("model_code.R")


# Create a data frame giving the data for a species and subregion.
make_dat_df <- function(species, subregion, include_ll = FALSE) {
  # Gather data
  this_obs <- ca_obs %>% filter(center == subregion, 
                                name_clean == species)
  this_cl <- ca_checklists %>% 
    filter(center == subregion) %>% 
    left_join(this_obs[, c("SAMPLING.EVENT.IDENTIFIER", "total_count")], 
              by = "SAMPLING.EVENT.IDENTIFIER") %>% 
    mutate(total_count = ifelse(!is.na(total_count), total_count, 0))
  
  loc_ids <- this_cl %>% 
    count(gx, gy) %>% 
    arrange(-n) %>% 
    mutate(locID = row_number())
  
  this_cl <- this_cl %>% left_join(
    loc_ids[, c("gx", "gy", "locID")], 
    by = c("gx", "gy")
  ) %>% 
    arrange(locID)
  
  if (include_ll) {
    dat_df <- this_cl %>% 
      mutate(#lon = scale(LONGITUDE),
        #lat = scale(LATITUDE),
        yday = scale(lubridate::yday(OBSERVATION.DATE)),
        tod = scale(as.numeric(TIME.OBSERVATIONS.STARTED)),
        duration = scale(DURATION.MINUTES),
        distance = scale(EFFORT.DISTANCE.KM),
        protocol = as.factor(PROTOCOL.CODE),
        num_observers = scale(NUMBER.OBSERVERS),
        year = as.factor(year),
        elevation = scale(elevation),
        precip = scale(precip),
        tmax = scale(tmax),
        pct_veg_comb = scale(pct_veg_comb),
        pct_veg = scale(pct_veg),
        pct_water = scale(pct_water),
        pct_urban = scale(pct_urban),
        pct_ag = scale(pct_ag),
        pct_tree = scale(pct_tree)
      ) %>% 
      mutate(tod_sq = tod^2,
             yday_sq = yday^2) %>% 
      select(total_count, locID, yday, yday_sq, tod, tod_sq, duration, 
             distance, protocol, num_observers, year, elevation, precip, tmax,
             pct_veg, pct_water, pct_tree, pct_ag, gx, gy)
  } else {
    dat_df <- this_cl %>% 
      mutate(#lon = scale(LONGITUDE),
        #lat = scale(LATITUDE),
        yday = scale(lubridate::yday(OBSERVATION.DATE)),
        tod = scale(as.numeric(TIME.OBSERVATIONS.STARTED)),
        duration = scale(DURATION.MINUTES),
        distance = scale(EFFORT.DISTANCE.KM),
        protocol = as.factor(PROTOCOL.CODE),
        num_observers = scale(NUMBER.OBSERVERS),
        year = as.factor(year),
        elevation = scale(elevation),
        precip = scale(precip),
        tmax = scale(tmax),
        pct_veg_comb = scale(pct_veg_comb),
        pct_veg = scale(pct_veg),
        pct_water = scale(pct_water),
        pct_urban = scale(pct_urban),
        pct_ag = scale(pct_ag),
        pct_tree = scale(pct_tree)
      ) %>% 
      mutate(tod_sq = tod^2,
             yday_sq = yday^2) %>% 
      select(total_count, locID, yday, yday_sq, tod, tod_sq, duration, 
             distance, protocol, num_observers, elevation, precip, tmax,
             pct_veg, pct_water, pct_tree, pct_ag)
  }
  
  dat_df
}

# Convert a column in a dataframe to wide form (each row = 1 location). Used
# for refitting in unmarked
get_var_wide <- function(df, col) {
  df <- df[, c("locID", col)]
  colnames(df)[2] <- "target_col"
  
  df <- df %>% 
    group_by(locID) %>% 
    mutate(id = paste0("X", str_pad(1:n(), 3, pad = "0"))) %>%
    spread(key = id, value = target_col) %>% 
    arrange(locID) %>% 
    ungroup() %>% 
    dplyr::select(-locID) %>% 
    as.matrix()
} 



# Get AIC from nimble optim object
AIC_optim <- function(optim_fit) {
  if (!is.list(optim_fit)) return(Inf)
  return(2 * length(optim_fit$par) - 2 * optim_fit$value)
}

# Get AIC from a fit unmarked model
AIC_unmarked <- function(umfit) {
  SEs <- c(summary(best_fit)$det$`Std. Error`,
           summary(best_fit)$state$`Std. Error`)
  if (any(is.na(SEs)) || any(is.nan(SEs))) return(Inf)
  return(umfit@AIC)
}

# Get AIC from a (potentially) fit glmmTMB object
AIC_tmb <- function(tmb_fit) {
  if (!is.list(tmb_fit)) return(Inf) # catch NAs
  
  if (is.na(AIC(tmb_fit))) return(Inf)
  
  AIC(tmb_fit)
}

# Make a GLMM formula from a row of a model specification matrix
glmm_formula_from_row <- function(mrow, accept_det_covs, accept_abd_covs) {
  
  names(mrow)[names(mrow) == "tod*tod"] <- "tod_sq"
  names(mrow)[names(mrow) == "yday*yday"] <- "yday_sq"
  
  if (sum(mrow) > 0) {
    f <- paste0("total_count ~ ", paste(names(mrow)[which(mrow == 1)], collapse = " + "),
                    paste0(" + ", paste0(accept_det_covs, collapse = "+"), " + ",
                           paste0(accept_abd_covs, collapse = "+"), " + (1|locID)"))
  } else f <- paste0("total_count ~ ", paste0(accept_det_covs, collapse = "+"), " + ",
                     paste0(accept_abd_covs, collapse = "+"), " + (1|locID)")

  as.formula(f)
}

# Make an unmarked formula from a row of a model specification matrix
unmarked_formula_from_row <- function(mrow, accept_det_names, accept_abd_names) {

  thisnames <- names(mrow)[mrow == 1]
  
  terms <- lapply(str_split(thisnames, pattern = "-"), function(x) x[[1]]) %>% 
    unlist()
  type <- lapply(str_split(thisnames, pattern = "-"), function(x) x[[2]]) %>% 
    unlist()
  
  
  terms[terms == "tod*tod"] <- "tod_sq"
  terms[terms == "yday*yday"] <- "yday_sq"
  
  det_f <- paste0("~", paste(c(accept_det_names, terms[type == "det"]), collapse = " + "))
  abd_f <- paste0("~", paste(c(accept_abd_names, terms[type == "abd"]), collapse = " + "))
  
  return(as.formula(paste(det_f, abd_f)))
}


# Turn a row from a model specification matrix into a set of NIMBLE model nodes
# for optimizing a specific set of covariates
row_to_nodes <- function(this_row, mixture, dat_df, 
                         accept_det_names, accept_abd_names) {
  terms <- names(this_row)
  
  int_inds <- which(grepl("\\*", terms))
  for (i in int_inds) {
    if (this_row[i] == 1) {
      sp <- unlist(str_split(terms[i], "-"))
      type <- sp[2]
      first_orders <- unlist(str_split(sp[1], "\\*"))
      this_row[grepl(first_orders[1], terms) &
               grepl(type, terms) & !grepl("\\*", terms)] <- 1
      this_row[grepl(first_orders[2], terms) &
               grepl(type, terms) & !grepl("\\*", terms)] <- 1
    }
  }
  
  
  abd_terms <- terms[grepl("abd", terms)]
  det_terms <- terms[grepl("det", terms)]
  
  intercepts <- c("p_coeffs[1]", "abund_coeffs[1]")
  
  this_target_nodes <-
    c(intercepts, # Always intercepts
      if (mixture %in% c("B-NB", "BB-NB")) "ltheta",
      if (mixture %in% c("BB-P", "BB-NB")) "log_s",
      paste0("p_coeffs[", 2:(1 + length(accept_det_names)), "]"), # accept det covariates
      paste0("abund_coeffs[", 2:(1 + length(accept_abd_names)), "]"), # accept abd covariates
      paste0("abund_coeffs[", (2 + length(accept_abd_names)):(length(abd_terms) + length(accept_abd_names) + 1), "]")[
           which(this_row[names(this_row) %in% abd_terms] == 1)
         ],
      paste0("p_coeffs[", (2 + length(accept_det_names)):(length(det_terms) + 1 + length(accept_det_names)), "]")[
           which(this_row[names(this_row) %in% det_terms] == 1)
         ])
  return(this_target_nodes)
}

# A version of nimble::setAndCalculate that returns large negative number if NaN
# is produced so that optimization doesn't break
noNanSetAndCalculate <- nimbleFunction(
  name = 'noNanSetAndCalculate',
  setup = function(model, targetNodes) {
    targetNodesAsScalar <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(targetNodes)
  },
  run = function(targetValues = double(1)) {
    values(model, targetNodesAsScalar) <<- targetValues
    lp <- calculate(model, calcNodes)
    returnType(double())
    if (!is.nan(lp)) return(lp)
    else return(-1e6)
  }
)


# Function to optimize a compiled N-mixture model. Takes arguments to set
# all parameters to 0, not just the ones being optimized, to "reset" the model
Cnmix_model_optim <- function(Cnmix_model, nmix_model, this_target_nodes,
                              all_zero_nodes, zero_assignments) {
  nimble::values(nmix_model, all_zero_nodes) <- zero_assignments
  nimble::values(Cnmix_model, all_zero_nodes) <- zero_assignments
  sc <- noNanSetAndCalculate(model = nmix_model, 
                             targetNodes = this_target_nodes)
  Csc <- suppressMessages(compileNimble(sc))
  fitct <- 0
  
  lower <- rep(-Inf, length(this_target_nodes))
  lower[this_target_nodes == "p_coeffs[1]"] <- logit(0.001)
  lower[this_target_nodes %in% c("ltheta", "log_s")] <- -10
  
  upper <- rep(Inf, length(this_target_nodes))
  # upper[this_target_nodes == "p_coeffs[1]"] <- logit(0.999)
  upper[this_target_nodes %in% c("ltheta", "log_s")] <- 10
  

  while (fitct < 5) {
    fitct <- fitct + 1
    this_fit <- tryCatch(
      optim(nimble::values(Cnmix_model, this_target_nodes),
            Csc$run, hessian = TRUE, method = "L-BFGS-B",
            lower = lower, upper = upper,
            control = list(fnscale = -1, maxit = 100000, trace = 0)),
      error = function(err) {
        return(NA)
      }
    )
    if (is.list(this_fit)) {
      fitct <- Inf
    } else {
      nimble::values(Cnmix_model, all_zero_nodes) <- zero_assignments
      nimble::values(nmix_model, all_zero_nodes) <- zero_assignments
      
      nimble::values(Cnmix_model, this_target_nodes) <- 
        nimble::values(Cnmix_model, this_target_nodes) + 
        rnorm(n = length(this_target_nodes), sd = 0.5)
    }
  }
  
  
  if (!is.na(this_fit) &&
      !any(this_fit$par[this_target_nodes %in% c("ltheta", "log_s")] == -10)) {
    tryCatch(
      this_fit <- optim(this_fit$par,
            Csc$run, hessian = TRUE, method = "L-BFGS-B",
            lower = lower, upper = upper,
            control = list(fnscale = -1, maxit = 100000, trace = 0)),
      error = function(e) {NA}
    )
  }
  
  return(this_fit)
}

# setAndCalculate reparameterized for use in uncertainty quantification
scReparam <- function(x, Csc, two_ints = TRUE, logit_p = NA,
                      ltheta = NULL, log_s = NULL,
                      this_target_nodes) {
  
  if (two_ints) {
    par <- x
    par[1] <- logit(exp( (x[1] + x[2]) / 2 ))
    par[2] <- (x[1] - x[2]) / 2
  } else {
    par <- c(NA, x)
    par[1] <- logit_p # logit(p)
    par[2] <- x[1] - log(expit(logit_p)) # log(p*lam) - log(p)
  }
  
  # Add fixed ltheta and log_s back in to par, if provided
  if (!is.null(ltheta)) {
    par <- c(par[1:2], ltheta, par[3:length(par)])
  }
  if (!is.null(log_s)) {
    if ("ltheta" %in% this_target_nodes) {
      par <- c(par[1:3], log_s, par[4:length(par)])
    } else {
      par <- c(par[1:2], log_s, par[3:length(par)])
    }
  }
  
  Csc$run(par)
}

# Get standard errors for parameters from fitted Cnmix_model object
Cnmix_model_SEs <- function(Cnmix_model, nmix_model, 
                            this_target_nodes, fit_values,
                            all_zero_nodes, zero_assignments) {
  nimble::values(Cnmix_model, all_zero_nodes) <- zero_assignments
  nimble::values(nmix_model, all_zero_nodes) <- zero_assignments
  sc <- noNanSetAndCalculate(model = nmix_model, 
                             targetNodes = this_target_nodes)
  Csc <- suppressMessages(compileNimble(sc))
  
  logit_p <- NULL
  two_ints <- TRUE
  ltheta <- NULL
  log_s <- NULL
  
  if (expit(fit_values[1]) >= 0.00125) {
    x <- fit_values
    x[1] <- log(expit(fit_values[1])) + fit_values[2] # log(p) + log(lam)
    x[2] <- log(expit(fit_values[1])) - fit_values[2] # log(p) - log(lam)
  } else {
    x <- fit_values
    x[1] <- log(expit(fit_values[1])) + fit_values[2] # log(p) + log(lam)
    x <- x[-2]
    logit_p <- fit_values[1] # store logit p
    two_ints <- F
  }
  
  # If ltheta was estimated at boundary, give it a fixed value
  if ("ltheta" %in% this_target_nodes &&
      fit_values[this_target_nodes == "ltheta"] %in% c(10, -10)) {
    ltheta <- fit_values[this_target_nodes == "ltheta"]
    x <- x[-which(this_target_nodes == "ltheta")]
  }
  
  # same for log_s
  if ("log_s" %in% this_target_nodes &&
      fit_values[this_target_nodes == "log_s"] %in% c(10, -10)) {
    log_s <- fit_values[this_target_nodes == "log_s"]
    x <- x[-which(this_target_nodes == "log_s")]
  }
  
  
  hess <- -hessian(scReparam, x, Csc = Csc, 
                   two_ints = two_ints, 
                   logit_p = logit_p,
                   ltheta = ltheta,
                   log_s = log_s,
                   this_target_nodes = this_target_nodes)
  
  SEs <- sqrt(diag(solve(hess)))
  
  if (!is.null(logit_p)) {
    SEs <- c(SEs[1], NA, SEs[2:length(SEs)])
  }
  if (!is.null(ltheta)) {
    SEs <- c(SEs[1:2], NA, SEs[3:length(SEs)])
  }
  if (!is.null(log_s)) {
    if ("ltheta" %in% this_target_nodes) {
      SEs <- c(SEs[1:3], NA, SEs[4:length(SEs)])
    } else {
      SEs <- c(SEs[1:2], NA, SEs[3:length(SEs)])
    }
  }
  
  return(SEs)
}


# Run forward AIC for GLMMs.
# dat_df: a data frame as returned by make_dat_df
# covs_tbl: describes which parameters in dat_df are being considered for fwd AIC
# fam: distribution for GLMM (nbinom2 or poisson)
GLMM_fwd_AIC <- function(dat_df, covs_tbl, fam, 
                         species = NA, subregion = NA,
                         verbose = F, parallelize = TRUE) {
  start_time <- Sys.time()
  
  accept_det_names <- clean_vec_names(c("duration", "num_observers"), dat_df)
  accept_abd_names <- clean_vec_names(c("elevation", "precip", "tmax"), dat_df)
  
  models_tried <- as.data.frame(matrix(
    0, ncol = nrow(covs_tbl), nrow = 1
  ))
  colnames(models_tried) <- covs_tbl$param
  null_formula <- glmm_formula_from_row(models_tried[1,1:nrow(covs_tbl)],
                                        accept_det_names, accept_abd_names)
  models_tried$AIC <- AIC_tmb(glmmTMB(null_formula, data = dat_df, family = fam))
  best_row <- models_tried[1,1:nrow(covs_tbl)]
  best_index <- 1
  
  done <- FALSE
  nmodel <- 0
  # First order terms loop
  while (!done) {
    last_best <- best_index
    covs_to_test <- which(!covs_tbl$accepted & !covs_tbl$interaction)
    
    if (length(covs_to_test) > 0) {
      if (parallelize) {
        plan(multicore)

        parAICs <- future_map_dbl(covs_to_test, function(x) {
          this_row <- best_row
          this_row[x] <- 1
          
          this_formula <- glmm_formula_from_row(this_row, accept_det_names, 
                                                accept_abd_names)
          
          this_fit <- glmmTMB(this_formula, data = dat_df, family = fam)
          return(AIC_tmb(this_fit))
        }) %>% unlist()
      } else {
        parAICs <- lapply(covs_to_test, function(x) {
          nmodel <<- nmodel + 1
          if (verbose) cat("\nModel ", nmodel)
          
          this_row <- best_row
          this_row[x] <- 1
          
          this_formula <- glmm_formula_from_row(this_row, accept_det_names,
                                                accept_abd_names)
          
          this_fit <- glmmTMB(this_formula, data = dat_df, family = fam)
          return(AIC_tmb(this_fit))
        }) %>% unlist()
      }
      
      for (i in 1:length(parAICs)) {
        this_row <- best_row
        this_row[covs_to_test[i]] <- 1
        this_row$AIC <- parAICs[i]
        models_tried <- rbind(models_tried, this_row)
      }
      
      best_index <- which.min(models_tried$AIC)
      if (best_index == last_best) {
        done <- TRUE
      } else {
        best_row <- models_tried[best_index, 1:nrow(covs_tbl)]
        covs_tbl$accepted <- as.logical(best_row)
      }
    } else {
      done <- TRUE
    }
  }
  
  done <- FALSE
  # Second order terms loop
  while (!done) {
    last_best <- best_index
    covs_to_test <- which(!covs_tbl$accepted & covs_tbl$interaction)
    
    if (length(covs_to_test) > 0) {
      if (parallelize) {
        plan(multicore)
        
        parAICs <- future_map_dbl(covs_to_test, function(x) {
          this_row <- best_row
          this_row[x] <- 1
          
          this_formula <- glmm_formula_from_row(this_row, accept_det_names,
                                                accept_abd_names)
          
          this_fit <- glmmTMB(this_formula, data = dat_df, family = fam)
          return(AIC_tmb(this_fit))
        }) %>% unlist()
      } else {
        parAICs <- lapply(covs_to_test, function(x) {
          nmodel <<- nmodel + 1
          if (verbose) cat("\nModel ", nmodel)
          
          this_row <- best_row
          this_row[x] <- 1
          
          this_formula <- glmm_formula_from_row(this_row, accept_det_names,
                                                accept_abd_names)
          
          this_fit <- glmmTMB(this_formula, data = dat_df, family = fam)
          return(AIC_tmb(this_fit))
        }) %>% unlist()
      }
      
      for (i in 1:length(parAICs)) {
        this_row <- best_row
        this_row[covs_to_test[i]] <- 1
        this_row$AIC <- parAICs[i]
        models_tried <- rbind(models_tried, this_row)
      }
      
      best_index <- which.min(models_tried$AIC)
      if (best_index == last_best) {
        done <- TRUE
      } else {
        best_row <- models_tried[best_index, 1:nrow(covs_tbl)]
        covs_tbl$accepted <- as.logical(best_row)
      }
    } else {
      done <- TRUE
    }
  }
  
  best_fit <- glmmTMB(glmm_formula_from_row(best_row, accept_det_names, accept_abd_names), 
                      data = dat_df, family = fam)
  
  coeffs_df <- summary(best_fit)$coefficients$cond %>% as.data.frame() %>% 
    mutate(param = rownames(summary(best_fit)$coefficients$cond)) %>% 
    mutate(est = Estimate, se = `Std. Error`) %>% 
    select(param, est, se)
  rownames(coeffs_df) <- NULL
  
  if (identical(fam, nbinom1) || identical(fam, nbinom2)) {
    sigma <- summary(best_fit)$sigma
  } else {
    sigma <- NA
  }

  end_time <- Sys.time()
  
  notes <- list()
  
  best_model_rtn_info <- list(
    species = species,
    subregion = subregion,
    fam = fam,
    AIC = AIC_tmb(best_fit),
    models_tried = models_tried,
    fit = best_fit,
    best_row = best_row,
    coefficients = coeffs_df,
    dispersion_param = sigma,
    var_locID = summary(best_fit)$varcor$cond$locID[[1]],
    time_taken = start_time - end_time,
    notes = notes
  )
  best_model_rtn_info
}

# Function to produce a NIMBLE model for the appropriate N-mixture, relevant data
make_nmix_model <- function(dat_df, vars, mixture, K) {
  accept_det_names <- clean_vec_names(c("duration", "num_observers"), dat_df)
  accept_abd_names <- clean_vec_names(c("elevation", "precip", "tmax"), dat_df)
  
  det_var_final_order <- c(accept_det_names, gsub("-det", "", vars[grepl("-det", vars)]))
  abd_var_final_order <- c(accept_abd_names, gsub("-abd", "", vars[grepl("-abd", vars)]))
  
  index_start <- numeric(max(dat_df$locID))
  index_end <- numeric(max(dat_df$locID))
  index_start[1] <- 1
  index_end[max(dat_df$locID)] <- nrow(dat_df)
  current_site <- 1
  for (i in 1:nrow(dat_df)) {
    this_site <- dat_df$locID[i]
    if (this_site != current_site) {
      index_end[current_site] <- i - 1
      index_start[this_site] <- i
      current_site <- this_site
    }
  }
  
  nSite1Obs <- sum(index_start == index_end)
  
  # Make the detection data matrix incl. interaction columns
  det_1st_mtx <- dat_df[, c(accept_det_names,
                            gsub("-det", "", 
                                 vars[grepl("det", vars) & !grepl("\\*", vars)]))]
  abd_1st_mtx <- dat_df[, c(accept_abd_names,
                            gsub("-abd", "", 
                                 vars[grepl("abd", vars) & !grepl("\\*", vars)]))]
  
  int_names <- vars[grepl("det", vars) & grepl("\\*", vars)] %>% 
    gsub(pattern="-det", replacement="")
  int_col_list <- list()
  if (length(int_names) > 0 ) {
    for (i in 1:length(int_names)) {
      depends <- strsplit(int_names[i], "\\*")[[1]]
      int_col_list[[i]] <- 
        dat_df[, depends[1]] * dat_df[, depends[2]]
    }
    int_mtx <- do.call(cbind, int_col_list)
    colnames(int_mtx) <- int_names
    p_dat <- cbind(det_1st_mtx, int_mtx)
  } else  {
    p_dat <- det_1st_mtx
  }
  abund_dat <- abd_1st_mtx
  
  p_dat     <- p_dat[, det_var_final_order]
  abund_dat <- abund_dat[, abd_var_final_order]
  
  use_code <- switch(mixture,
                     `B-P` = BP_nmix_code_ebd,
                     `B-NB` = BNB_nmix_code_ebd,
                     `BB-P` = BBP_nmix_code_ebd,
                     `BB-NB` = BBNB_nmix_code_ebd)
  
  nmix_model <- nimbleModel(
    code = use_code,
    data = list(abund_dat = abund_dat,
                p_dat = p_dat,
                K = K,
                y = as.numeric(dat_df$total_count)),
    constants = list(
      nobs = nrow(dat_df),
      npcov = ncol(p_dat) + 1,
      nlcov = ncol(abund_dat) + 1,
      nsite = length(index_start),
      index_start = index_start,
      index_end = index_end,
      nSite1Obs = nSite1Obs),
    inits = list(
      p_coeffs = c(logit(0.1), rep(0, ncol(p_dat))),
      abund_coeffs = 
        c(log(max(dat_df$total_count)), 
          rep(0, ncol(abund_dat))),
      ltheta = 0, log_s = 0
    ), calculate = FALSE)
  Cnmix_model <- compileNimble(nmix_model)
  
  return(list(
    nmix_model = nmix_model,
    Cnmix_model = Cnmix_model
  ))
}

# This function takes the standard parameters that are a priori accepted
# and makes sure they're not too correlated. If two are, the latter one is
# dropped. This is more of an issue for the automatically accepted ones because
# otherwise it'd largely be taken care of by stepwise AIC.
clean_vec_names <- function(accept_names, dat_df) {
  
  combos <- as.data.frame(t(combn(accept_names, 2)))
  dropped <- data.frame(param = accept_names, dropped = F)
  
  for (i in 1:nrow(combos)) {
    if (!dropped$dropped[dropped$param == combos$V1[i]] &&
        !dropped$dropped[dropped$param == combos$V2[i]]) {
      if (abs(cor(dat_df[, combos$V1[i]], dat_df[, combos$V2[i]])) > 0.9) {
        dropped$dropped[dropped$param == combos$V2[i]] <- TRUE
      }
    }
  }
  
  return(as.character(dropped$param[!dropped$dropped]))
}

# Run forward AIC for N-mixture models.
Nmix_fwd_AIC <- function(dat_df, covs_tbl, mixture, K, species = NA, 
                         subregion = NA, verbose = F, nCores, clearComp = T) {
  on.exit({
    if (exists("in.cl")) stopCluster(in.cl)
  })
  start_time <- Sys.time()
  
  models_tried <- as.data.frame(matrix(
    0, ncol = nrow(covs_tbl), nrow = 1
  ))
  colnames(models_tried) <- paste(covs_tbl$param, covs_tbl$type, sep = "-")
  models_tried$AIC <- NA
  
  accept_det_names <- clean_vec_names(c("duration", "num_observers"), dat_df)
  accept_abd_names <- clean_vec_names(c("elevation", "precip", "tmax"), dat_df)

  if (nCores > 0) {
    cluster.env <- new.env()
    cluster.env$dat_df <- dat_df
    cluster.env$mixture <- mixture
    cluster.env$models_tried <- models_tried
    cluster.env$K <- K
    cluster.env$accept_det_names <- accept_det_names
    cluster.env$accept_abd_names <- accept_abd_names
    cluster.env$verbose <- verbose
    
    in.cl <- makeCluster(nCores)
    capture <- clusterEvalQ(in.cl, source("03_fit_one_ssr_final.R")) 
    clusterExport(in.cl, 
                  varlist = list("dat_df", "models_tried", "mixture", "K",
                                 "accept_det_names", "accept_abd_names", "verbose"), 
                  envir = cluster.env)
  }
  
  start_time <- Sys.time()
  if (!(mixture %in% c("B-P", "B-NB", "BB-P", "BB-NB"))) {
    stop("Unrecognized mixture ", mixture, ".")
  }
  
  # Make the model object
  if (nCores > 0) {
    capture <- clusterEvalQ(in.cl, {
      modellist <- 
        make_nmix_model(dat_df, vars = colnames(models_tried)[-ncol(models_tried)], 
                        mixture, K)
      nmix_model <- modellist$nmix_model
      Cnmix_model <- modellist$Cnmix_model
      Cnmix_model$calculate()
      all_zero_nodes <- c(
        nmix_model$expandNodeNames("abund_coeffs"), 
        nmix_model$expandNodeNames("p_coeffs"),
        if (mixture %in% c("B-NB", "BB-NB")) "ltheta",
        if (mixture %in% c("BB-P", "BB-NB")) "log_s"
      )
      zero_assignments <- nimble::values(Cnmix_model, all_zero_nodes)
    })
  } else {
    modellist <- 
      make_nmix_model(dat_df, vars = colnames(models_tried)[-ncol(models_tried)], 
                      mixture, K)
    nmix_model <- modellist$nmix_model
    Cnmix_model <- modellist$Cnmix_model
    Cnmix_model$calculate()
    all_zero_nodes <- c(
      nmix_model$expandNodeNames("abund_coeffs"), 
      nmix_model$expandNodeNames("p_coeffs"),
      if (mixture %in% c("B-NB", "BB-NB")) "ltheta",
      if (mixture %in% c("BB-P", "BB-NB")) "log_s"
    )
    zero_assignments <- nimble::values(Cnmix_model, all_zero_nodes)
  }
  
  if (verbose) cat("Beginning selection loop.")

  models_tried[1,] <- c(rep(0, nrow(covs_tbl)),
                        NA)
  null_nodes <- row_to_nodes(models_tried[1,1:nrow(covs_tbl)], 
                             mixture = mixture, dat_df = dat_df,
                             accept_det_names, accept_abd_names)
  if (nCores > 0) {
    
    null_fit <- parLapply(in.cl, list(null_nodes), function(x) {
      Cnmix_model_optim(Cnmix_model, nmix_model,
                        x, all_zero_nodes,
                        zero_assignments)
    })[[1]]
  } else {
    null_fit <- Cnmix_model_optim(Cnmix_model, nmix_model,
                                  null_nodes, all_zero_nodes,
                                  zero_assignments)
  }
  
  models_tried$AIC[1] <- AIC_optim(null_fit)
  best_row <- models_tried[1,]
  best_fit <- null_fit
  
  done <- FALSE
  interactions <- FALSE
  ct <- 1
  round <- 0
  
  while (!done) {
    round <- round + 1
    last_best <- best_row
    covs_to_test <- covs_tbl %>% 
      filter(!accepted) %>% 
      filter(grepl(pattern = "\\*", x = param) == interactions)
    if (nrow(covs_to_test) == 0) {
      if (interactions) done <- TRUE
      else interactions <- TRUE
    } else {
      to_test_list <- split(covs_to_test, seq(nrow(covs_to_test)))
      for (i in 1:length(to_test_list)) {
        ct <- ct + 1
        models_tried[ct,] <- c(as.numeric(
          covs_tbl$accepted |
          (covs_tbl$param == to_test_list[[i]]$param & 
           covs_tbl$type  == to_test_list[[i]]$type)), 
          NA)
      }
      
      if (nCores > 0) {
        models_tried$AIC[(ct - length(to_test_list) + 1):ct] <- 
          parLapply(cl = in.cl, X = to_test_list, fun = function(x) {
            if (verbose) cat(paste0("\nTesting round ", round, " add ", x$param))
            this_row <- c(as.numeric(
              covs_tbl$accepted |
                covs_tbl$param == x$param), NA)
            names(this_row) <- 
              c(paste(covs_tbl$param, covs_tbl$type, sep = "-"), "AIC")
            this_nodes <- row_to_nodes(this_row, mixture = mixture,
                                       dat_df = dat_df,
                                       accept_det_names, accept_abd_names)
            
            this_fit <- Cnmix_model_optim(Cnmix_model, nmix_model, this_nodes, 
                                          all_zero_nodes, zero_assignments)
            AIC_optim(this_fit)
          }) %>% unlist()
        
      } else {
        models_tried$AIC[(ct - length(to_test_list) + 1):ct] <- 
          lapply(to_test_list, function(x) {
            if (verbose) cat(paste0("\nTesting round ", round, " add ", x$param, " ", x$type))
            this_row <- c(as.numeric(
              covs_tbl$accepted |
                covs_tbl$param == x$param), NA)
            names(this_row) <- 
              c(paste(covs_tbl$param, covs_tbl$type, sep = "-"), "AIC")
            this_nodes <- row_to_nodes(this_row, mixture = mixture, 
                                       dat_df = dat_df,
                                       accept_det_names, accept_abd_names)
            
            this_fit <- Cnmix_model_optim(Cnmix_model, nmix_model, this_nodes, 
                                          all_zero_nodes, zero_assignments)
            AIC_optim(this_fit)
          }) %>% unlist()
      }
            
      best_row <- models_tried[which.min(models_tried$AIC),]
      
      if (last_best$AIC == best_row$AIC) {
        if (interactions) done <- TRUE
        else interactions <- TRUE
      } else {
        covs_tbl$accepted[as.logical(best_row[1:nrow(covs_tbl)])] <- TRUE
        if (verbose) {
          cat(paste0("\nAccepted: ", 
              paste(covs_tbl$param[covs_tbl$accepted], 
                    collapse = " "),
              ". Best AIC: ", best_row$AIC))
        }
      }
    }
  }
  
  
  
  ### Done with forward AIC
  best_row <- models_tried[which.min(models_tried$AIC), ]
  best_nodes <- row_to_nodes(best_row, mixture, dat_df,
                             accept_det_names, accept_abd_names)
  if (nCores > 0) {
    
    best_fit <- parLapply(in.cl, list(best_nodes), function(x) {
      Cnmix_model_optim(Cnmix_model, nmix_model,
                        x, all_zero_nodes,
                        zero_assignments)
    })[[1]]
  } else {
    best_fit <- Cnmix_model_optim(Cnmix_model, nmix_model,
                                  best_nodes, 
                                  all_zero_nodes, zero_assignments)
  }
  
  terms <- names(best_row)
  
  int_inds <- which(grepl("\\*", terms))
  for (i in int_inds) {
    if (best_row[i] == 1) {
      sp <- unlist(str_split(terms[i], "-"))
      type <- sp[2]
      first_orders <- unlist(str_split(sp[1], "\\*"))
      best_row[grepl(first_orders[1], terms) &
                 grepl(type, terms) & !grepl("\\*", terms)] <- 1
      best_row[grepl(first_orders[2], terms) &
                 grepl(type, terms) & !grepl("\\*", terms)] <- 1
    }
  }
  
  abd_terms <- terms[grepl("abd", terms)]
  det_terms <- terms[grepl("det", terms)]
  

  params_chosen <-
    c("log(p * lambda)", "log(p / lambda)", # Always intercepts
      if (mixture %in% c("B-NB", "BB-NB")) "ltheta",
      if (mixture %in% c("BB-P", "BB-NB")) "log_s",
      accept_det_names,
      accept_abd_names,
      abd_terms[
        which(best_row[colnames(best_row) %in% abd_terms] == 1)
      ],
      det_terms[
        which(best_row[colnames(best_row) %in% det_terms] == 1)
      ])
  
  
  notes <- list()
  tryCatch({
    if (nCores > 0) {
      ses <- parLapply(in.cl, list(best_nodes), x2 = best_fit$par,
                       function(x, x2) {
         Cnmix_model_SEs(Cnmix_model, nmix_model, 
                         x, x2,
                         all_zero_nodes, zero_assignments)
      })[[1]]
    } else {
      ses <- Cnmix_model_SEs(Cnmix_model, nmix_model, 
                             best_nodes, 
                             best_fit$par,
                             all_zero_nodes, zero_assignments)
    }
  }, error = function(err) {} )
  if (!exists("ses")) {
    ses <- rep(NA, length(params_chosen))
  }
  
  if (length(ses) == length(best_fit$par)) {
    coeffs_df <- data.frame(
      param = params_chosen,
      est = c(log(expit(best_fit$par[1])) + best_fit$par[2], 
              log(expit(best_fit$par[1])) - best_fit$par[2],
              best_fit$par[3:length(best_fit$par)]),
      se = ses
    )
  } else if (length(ses) + 1 == length(best_fit$par)) {
    coeffs_df <- data.frame(
      param = params_chosen,
      est = c(log(expit(best_fit$par[1])) + best_fit$par[2], 
              log(expit(best_fit$par[1])) - best_fit$par[2],
              best_fit$par[3:length(best_fit$par)]),
      se = c(ses[1], NA, ses[2:length(ses)])
    )
  }

  end_time <- Sys.time()
  
  best_model_rtn_info <- list(
    species = species,
    subregion = subregion,
    mixture = mixture,
    AIC = AIC_optim(best_fit),
    models_tried = models_tried,
    best_row = best_row,
    fit = best_fit,
    coefficients = coeffs_df,
    K = K,
    time_taken = start_time - end_time,
    notes = notes
  )
  
  if (clearComp && nCores == 0) nimble:::clearCompiled(nmix_model)
  
  best_model_rtn_info
}

# Run forward AIC for all six model flavors on the specified SSR.
fit_one_ssr <- function(subregion, species, checklists, obs, nCores, 
                        overwrite = c("GLMM_Nbin", "GLMM_Pois", 
                                      "Nmix_BP", "Nmix_BNB",
                                      "Nmix_BBP", "Nmix_BBNB"), 
                        verbose = TRUE, glmm_only = FALSE,
                        path_prefix = "output/onemodel_oneyear/") {

  if (!(paste0(subregion, species)) %in% 
      unique(paste0(obs$center, obs$name_clean))) {
    warning("No data for target ssr ", subregion, " ", species, ". Skipping.")
    return(NA)
  }

  dat_df <- make_dat_df(species, subregion)
  
  both_covs <- c("pct_water", "pct_tree", "pct_veg", "pct_ag")
  both_covs <- both_covs[
    unlist(map_lgl(both_covs, function(x) length(unique(unlist(dat_df[,x]))) > 1))]

  det_only_covs <- c("yday", "yday*yday", "tod", "tod*tod", "distance")
  
  nmix_covs_tbl <- data.frame(
    param = c(both_covs, both_covs, det_only_covs),
    accepted = F,
    type = c(rep("abd", length(both_covs)),
             rep("det", length(c(both_covs, det_only_covs)))),
    interaction = grepl("\\*", c(both_covs, both_covs, det_only_covs))
  )
  glmm_covs_tbl <- nmix_covs_tbl %>% 
    select(-type) %>% 
    distinct()
  
  if (!file.exists(paste0(path_prefix, "GLMM_Pois_", 
                          subregion, "_", species, ".RDS")) ||
      "GLMM_Pois" %in% overwrite) {
    if (verbose) cat("\n\n===== Poisson GLMM", species, subregion, "=====")
    glmm_pois_result <- GLMM_fwd_AIC(dat_df = dat_df, 
                                     covs_tbl = glmm_covs_tbl, 
                                     fam = poisson, 
                                     species = species,
                                     subregion = subregion, 
                                     verbose = verbose, 
                                     parallelize = nCores > 1)
    saveRDS(glmm_pois_result, 
            paste0(path_prefix, "GLMM_Pois_", 
                   subregion, "_", species, ".RDS"))
  }
  
  if (!file.exists(paste0(path_prefix, "GLMM_Nbin_", 
                          subregion, "_", species, ".RDS")) ||
      "GLMM_Nbin" %in% overwrite) {
    if (verbose) cat("\n\n===== Negative Binomial GLMM", species, subregion, "=====")
    glmm_nbin_result <- GLMM_fwd_AIC(dat_df = dat_df, 
                                     covs_tbl = glmm_covs_tbl, 
                                     fam = nbinom2, 
                                     species = species,
                                     subregion = subregion, 
                                     verbose = verbose, 
                                     parallelize = nCores > 0)
    saveRDS(glmm_nbin_result, 
            paste0(path_prefix, "GLMM_Nbin_", 
                   subregion, "_", species, ".RDS"))
  }
  
  if (glmm_only) return(NULL)
  
  K <- max(qpois(p = 0.999999, lambda = max(dat_df$total_count) * 50),
           200)
  if (!file.exists(paste0(path_prefix, "Nmix_BP_", 
                          subregion, "_", species, ".RDS")) ||
      "Nmix_BP" %in% overwrite) {
    if (verbose) cat("\n\n===== B-P N-mixture", species, subregion, "=====")
    nmix_bp_result <- Nmix_fwd_AIC(dat_df = dat_df, 
                                   covs_tbl = nmix_covs_tbl, 
                                   mixture = "B-P", 
                                   K = K, 
                                   species = species, 
                                   subregion = subregion, 
                                   verbose = verbose, 
                                   nCores = nCores)
    saveRDS(nmix_bp_result, paste0(path_prefix, "Nmix_BP_", 
                                   subregion, "_", species, ".RDS"))
  }
  
  if (!file.exists(paste0(path_prefix, "Nmix_BNB_", 
                          subregion, "_", species, ".RDS")) ||
      "Nmix_BNB" %in% overwrite) {
    if (verbose) cat("\n\n===== B-NB N-mixture", species, subregion, "=====")
    nmix_bnb_result <- Nmix_fwd_AIC(dat_df = dat_df, 
                                    covs_tbl = nmix_covs_tbl, 
                                    mixture = "B-NB", 
                                    K = K, 
                                    species = species, 
                                    subregion = subregion, 
                                    verbose = verbose, 
                                    nCores = nCores)
    saveRDS(nmix_bnb_result, paste0(path_prefix, "Nmix_BNB_", 
                                    subregion, "_", species, ".RDS"))
  }
  
  if (!file.exists(paste0(path_prefix, "Nmix_BBP_", 
                          subregion, "_", species, ".RDS")) ||
      "Nmix_BBP" %in% overwrite) {
    if (verbose) cat("\n\n===== BB-P N-mixture", species, subregion, "=====")
    nmix_bbp_result <- Nmix_fwd_AIC(dat_df = dat_df, 
                                    covs_tbl = nmix_covs_tbl, 
                                    mixture = "BB-P", 
                                    K = K, 
                                    species = species, 
                                    subregion = subregion, 
                                    verbose = verbose, 
                                    nCores = nCores)
    saveRDS(nmix_bbp_result, paste0(path_prefix, "Nmix_BBP_", 
                                    subregion, "_", species, ".RDS"))
  }
  
  if (!file.exists(paste0(path_prefix, "Nmix_BBNB_", 
                          subregion, "_", species, ".RDS")) ||
      "Nmix_BBNB" %in% overwrite) {
    if (verbose) cat("\n\n===== BB-NB N-mixture", species, subregion, "=====")
    nmix_bbnb_result <- Nmix_fwd_AIC(dat_df = dat_df, 
                                     covs_tbl = nmix_covs_tbl, 
                                     mixture = "BB-NB", 
                                     K = K, 
                                     species = species, 
                                     subregion = subregion, 
                                     verbose = verbose, 
                                     nCores = nCores)
    saveRDS(nmix_bbnb_result, paste0(path_prefix, "Nmix_BBNB_", 
                                     subregion, "_", species, ".RDS"))
  }
  
}


