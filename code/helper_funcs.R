
causal_blb <- function(data, y_formula, prop_formula, y_method, prop_method, r = 100, subsets = 5, b = NULL, m_tune = NULL, prop_tune = NULL,
                       y_args = NULL, prop_args = NULL, cores = 1){
  
  if(data.table::is.data.table(data) == FALSE){
    data <- as.data.table(data)
  }
  
  if(is.null(b)){
    b <- nrow(data)^0.7
  }
  W <- all.vars(prop_formula)[1]
  Y <- all.vars(y_formula)[1]
  
  data$Tr2 <- factor(data[[W]])
  levels(data$Tr2) <- make.names(levels(data$Tr2))
  prop_formula <- update.formula(as.formula(prop_formula), Tr2 ~ .)
  
  
  m <- do.call(caret::train, c(list(form = as.formula(y_formula),
                                    data = data,
                                    method = y_method,
                                    trControl = caret::trainControl(method = 'none'),
                                    tuneGrid = m_tune), y_args))
  
  g <- do.call(caret::train, c(list(form = as.formula(prop_formula),
                                           data = data,
                                           method = prop_method,
                                           trControl = caret::trainControl(method = 'none', classProbs = TRUE),
                                           tuneGrid = prop_tune), prop_args))
  
  # Calculate full sample
  # Full sample
  prop_score <- predict(g, type = 'prob')[, 1]
  pred_data <- data.table::copy(data)
  pred_data[[W]] <- 1
  m1 <- predict(m, pred_data)
  pred_data[[W]] <- 0
  m0 <- predict(m, pred_data)
  pred_data[[W]] <- data[[W]]
  
  phi1_full <- (pred_data[[W]]/prop_score)*(pred_data[[Y]] - m1) + m1
  phi0_full <- (1 - pred_data[[W]])/(1 - prop_score)*(pred_data[[Y]] - m0) + m0
  
  tau_hat_full <- mean(phi1_full) - mean(phi0_full)
  partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = FALSE)
  
  blb_out <- parallel::mclapply(partitions, function(i){
    tmp_dat <- data[i]
    M <- rmultinom(n = r, size = n, prob = rep(1, b))
    
    blb_reps <- sapply(seq_len(r), function(bt){
      m <- do.call(caret::train, c(list(form = as.formula(y_formula),
                                        data = tmp_dat,
                                        method = y_method,
                                        trControl = caret::trainControl(method = 'none'),
                                        tuneGrid = m_tune, weights = M[, bt]), y_args))
      
      g <- do.call(caret::train, c(list(form = as.formula(prop_formula),
                                        data = tmp_dat,
                                        method = prop_method,
                                        trControl = caret::trainControl(method = 'none', classProbs = TRUE),
                                        tuneGrid = prop_tune, weights = M[, bt]), prop_args))
      
      
      # Calculate full sample
      # Full sample
      prop_score <- predict(g, type = 'prob')[, 1]
      pred_data <- data.table::copy(tmp_dat)
      pred_data[[W]] <- 1
      m1 <- predict(m, pred_data)
      pred_data[[W]] <- 0
      m0 <- predict(m, pred_data)
      pred_data[[W]] <- tmp_dat[[W]]

      phi1 <- M[, bt]*((pred_data[[W]]/prop_score)*(pred_data[[Y]] - m1) + m1)
      phi0 <- M[, bt]*((1 - pred_data[[W]])/(1 - prop_score)*(pred_data[[Y]] - m0) + m0)
      sum(phi1)/n - sum(phi0)/n
    })
    
    bias <- mean(blb_reps) - tau_hat_full
    blb_bias_correct <- blb_reps - bias
    perc_ci <- boot:::perc.ci(blb_bias_correct)
    return(data.table::data.table(lower_ci = perc_ci[4],
                      upper_ci = perc_ci[5]))
  }, mc.cores = cores)
  
  
  blb_out <- data.table::rbindlist(blb_out)
  blb_out <- blb_out[, .(lower_ci = mean(lower_ci),
                         upper_ci = mean(upper_ci))]
  blb_out
}

inv_logit <- function(x){
  1/(1 + exp(-x))
}

kangschafer3 <- function(n, te, beta_overlap = 0.5, sigma) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Z1 <- exp(X1/2)
  Z2 <- 0.25*(X1*X2)
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Z1, Z2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

# # Non-independent
# kangschafer3_dep <- function(n, te, beta_overlap = 0.5, sigma, rho = 0.8) {
#   Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
#   X <- mvtnorm::rmvnorm(n, sigma = Sigma)
#   X1 <- X[, 1]
#   X2 <- X[, 2]
#   prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
#   
#   Tr <- rbinom(n, 1, prt)
#   
#   Y0  <- X1 + X2 + rnorm(n, 0, sigma)
#   Y1  <- Y0 + te
#   y   <- Y0*(1-Tr) + Y1*(Tr)
#   out <- cbind(y, Tr, X1, X2, Y1, Y0)
#   out <- data.table::as.data.table(out)
#   return(out)
# }
# 
# # Multiple predictors
# kangschafer3_mult <- function(n, te, beta_overlap = 0.1, sigma, preds = 10) {
#   Sigma <- diag(preds)
#   X <- mvtnorm::rmvnorm(n, sigma = Sigma)
#   X_comb <- rowSums(X)
#   prt <- 1/(1 + exp(-beta_overlap*(X_comb)))
#   
#   Tr <- rbinom(n, 1, prt)
#   
#   Y0  <- X_comb + rnorm(n, 0, sigma)
#   Y1  <- Y0 + te
#   y   <- Y0*(1-Tr) + Y1*(Tr)
#   out <- cbind(y, Tr, X, Y1, Y0)
#   out <- data.table::as.data.table(out)
#   return(out)
# }
# 
# # Heterogeneity
# kangschafer3_het <- function(n, te, beta_overlap = 0.5, sigma) {
#   X1 <- rnorm(n, 0, 1)
#   X2 <- rnorm(n, 0, 1)
#   prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
#   
#   Tr <- rbinom(n, 1, prt)
#   
#   Y0  <- X1 + X2 + rnorm(n, 0, sigma)
#   Y1  <- Y0 + te + X1
#   y   <- Y0*(1-Tr) + Y1*(Tr)
#   out <- cbind(y, Tr, X1, X2, Y0, Y1)
#   out <- data.table::as.data.table(out)
#   return(out)
# }
# 
# # Misspecified logit model (Z1 and Z1 instead of X1 and X2)
# kangschafer3_mis <- function(n, te, beta_overlap = 0.5, sigma) {
#   X1 <- rnorm(n, 0, 1)
#   X2 <- rnorm(n, 0, 1)
#   
#   Z1 <- X1^3/2
#   Z2 <- 0.25*(X1 + X2)^2
#   
#   prt <- 1/(1 + exp(-beta_overlap*(Z1 + Z2)))
#   
#   Tr <- rbinom(n, 1, prt)
#   
#   Y0  <- Z1 + Z2 + rnorm(n, 0, sigma)
#   Y1  <- Y0 + te + Z1
#   y   <- Y0*(1-Tr) + Y1*(Tr)
#   out <- cbind(y, Tr, X1, X2, Y0, Y1)
#   out <- data.table::as.data.table(out)
#   return(out)
# }
# 
# 
# make_crossfit_estim <- function(data, y_method, prop_method, K, y_formula, prop_formula, m_tune, 
#                                 prop_tune, W, prop_args = list(), y_args = list(), crossfit){
#   if(!crossfit){
#     m <- do.call(caret::train, c(list(form = as.formula(y_formula),
#                                       data = data,
#                                       method = y_method,
#                                       trControl = caret::trainControl(method = 'none'),
#                                       tuneGrid = m_tune), y_args))
#     
#     prop_mod <- do.call(caret::train, c(list(form = as.formula(prop_formula),
#                                              data = data,
#                                              method = prop_method,
#                                              trControl = caret::trainControl(method = 'none', classProbs = TRUE),
#                                              tuneGrid = prop_tune), prop_args))
#     
#     pred_dat <- copy(data)
#     pred_dat[[W]] <- 0
#     y0 <- predict(m, pred_dat)
#     pred_dat[[W]] <- 1
#     y1 <- predict(m, pred_dat)
#     
#     data$m0 <- y0
#     data$m1 <- y1
#     ps <- predict(prop_mod, data, 'prob')[, 2]
#     data$ps <- ps
#     return(data)
#   } else{
#     part_idx <- seq_len(nrow(data))
#     folds <- split(part_idx, sample(rep(1:K, length.out = length(part_idx))))
#     cross_dat <- lapply(folds, function(test_idx){
#       train_idx <- setdiff(part_idx, test_idx)
#       train_dat <- data[train_idx]
#       test_dat <- data[-train_idx]
# 
#       m <- do.call(caret::train, c(list(form = as.formula(y_formula),
#                                         data = train_dat,
#                                         method = y_method,
#                                         trControl = caret::trainControl(method = 'none'),
#                                         tuneGrid = m_tune), y_args))
#       
#       prop_mod <- do.call(caret::train, c(list(form = as.formula(prop_formula),
#                                                data = train_dat,
#                                                method = prop_method,
#                                                trControl = caret::trainControl(method = 'none', classProbs = TRUE),
#                                                tuneGrid = prop_tune), prop_args))
#       
#       pred_dat <- copy(test_dat)
#       pred_dat[[W]] <- 0
#       y0 <- predict(m, pred_dat)
#       pred_dat[[W]] <- 1
#       y1 <- predict(m, pred_dat)
#       
#       test_dat$m0 <- y0
#       test_dat$m1 <- y1
#       ps <- predict(prop_mod, test_dat, 'prob')[, 2]
#       test_dat$ps <- ps
#       test_dat
#     })
#     return(rbindlist(cross_dat))
#   }
# 
# 
#   
# 
# }
# 
# process_partition <- function(partition_element, data, part_idx, b, n, y_method, prop_method, K, y_formula, prop_formula,
#                               m_tune, prop_tune, W, y_args, prop_args, boot_stat, r, Y, crossfit){
# 
#   part_dat <- data[partition_element]
#   crossfit_dat <- make_crossfit_estim(data = part_dat, 
#                                       y_method = y_method, 
#                                       prop_method = prop_method, 
#                                       K = K, 
#                                       y_formula = y_formula,
#                                       prop_formula = prop_formula, 
#                                       m_tune = m_tune, 
#                                       prop_tune = prop_tune,
#                                       W = W,
#                                       prop_args = prop_args,
#                                       y_args = y_args,
#                                       crossfit = crossfit)
#   
# 
#   
#   # CHECK OBSERVATIONS WITH EXTREME PS VALUES
#   chk <- data.table::between(crossfit_dat$ps, 0, 1, incbounds = FALSE)
#   assertthat::assert_that(all(chk))
#   # if(any(!chk)){
#   #   crossfit_dat <- crossfit_dat[chk]
#   #   b <- nrow(crossfit_dat)
#   #   part_idx <- seq_len(b)
#   # }
#   theta0 <- boot_stat(data = crossfit_dat, Tr = W, y = Y, return_var = TRUE)
#   boot_reps_n <- blb(data = crossfit_dat, r = r, idx = part_idx, size = n, 
#                    boot_stat = boot_stat, W = W, Y = Y, return_var = TRUE)
#   thetahat <- boot_reps_n$thetahat
#   sdthetahat <- boot_reps_n$sdthetahat
#   stud <- boot_reps_n$stud
#   boot_reps_n <- boot_reps_n$boot_reps
#   assertthat::assert_that(all(!is.na(boot_reps_n)))
#   # PERCENTILE
#   perc_ci <- boot:::perc.ci(t = boot_reps_n)[, 4:5]
#   # BCA
#   a <- jackknife_acceleration(data = crossfit_dat, statistic = boot_stat, Tr = W, y = Y)
#   z <- qnorm(sum(boot_reps_n < thetahat)/r)
#   zalpha <- (1 + c(-0.95, 0.95))/2
#   zalpha <- qnorm(zalpha)
#   adj.alpha <- pnorm(z + (z + zalpha)/(1 - a * (z + zalpha)))
#   bca_ci <- quantile(boot_reps_n, prob = adj.alpha)
#   # BASIC
#   basic_ci <- boot:::basic.ci(t0 = thetahat, t = boot_reps_n)[, 4:5]
#   # NORMAL
#   norm_ci <- boot:::norm.ci(var.t0 = NULL, t0 = thetahat, t = boot_reps_n)[, 2:3]
#   # Studentized
#   quants <- quantile(stud, probs = c(0.025, 0.975))
# 
#   dt <- data.table(estim = thetahat,
#                    estim_sub = theta0[1],
#                    sd_sub = sqrt(theta0[2]),
#                    sdthetahat = sdthetahat,
#                    lower = c(perc_ci[1],
#                              bca_ci[1],
#                              basic_ci[1],
#                              norm_ci[1],
#                              thetahat - quants[2]*sdthetahat,
#                              thetahat - 1.96*sdthetahat),
#                    upper = c(perc_ci[2],
#                              bca_ci[2],
#                              basic_ci[2],
#                              norm_ci[2],
#                              thetahat - quants[1]*sdthetahat,
#                              thetahat + 1.96*sdthetahat),
#                    boot_type = c('percentile', 'bca', 'basic', 'normal', 'studentized', 'asymptotic'))
# 
#   return(dt)
# }


make_partition <- function(n, subsets, b, disjoint = FALSE){
  part_idx <- seq(1, n, by = 1)
  if(disjoint){
    # Generate disjoint sets
    # Permute indices
    partition <- sample(part_idx, n)
    totality <- b*subsets
    stopifnot(totality <= n)
    # Exclude rest of sample
    partition <- partition[1:totality]
    partition <- split(partition, f = rep(1:subsets, each = b))
  } else{
    partition <- replicate(subsets, {
      sample(part_idx, size = b, replace = FALSE)
    }, simplify = FALSE)
  }
  partition
}

# blb <- function(data, r, idx, size, boot_stat, W, Y, return_var = TRUE){
#   boot_reps <- replicate(r, {
#     boot_idx <- sample(idx, size = size, replace = TRUE)
#     boot_dat <- data[boot_idx]
#     boot_stat(data = boot_dat, Tr = W, y = Y, return_var = return_var)
#   })
#   thetahat <- mean(boot_reps[1, ])
#   sdthetahat <- sd(boot_reps[1, ])
#   stud <- (boot_reps[1, ] - thetahat)/sqrt(boot_reps[2, ])
#   boot_reps <- boot_reps[1, ]
#   return(list(boot_reps = boot_reps, stud = stud, thetahat = thetahat, sdthetahat = sdthetahat))
# }
# 
# dbml_blb <- function(data, r, Y, W, K, subsets,
#                      disjoint = FALSE,
#                      boot_stat,
#                      b = NULL, 
#                      prop_method = 'ranger',
#                      prop_formula = NULL, 
#                      prop_grid = NULL,
#                      prop_control = caret::trainControl(),
#                      prop_args = list(),
#                      y_method = 'ranger',
#                      y_formula = NULL,
#                      y_grid = NULL,
#                      y_control = caret::trainControl(),
#                      y_args = list(),
#                      hp_tune = TRUE,
#                      crossfit = TRUE,
#                      aggregate = TRUE){
#   
#   if(data.table::is.data.table(data) == FALSE){
#     data <- data.table::as.data.table(data)
#   }
#   stopifnot(Y %in% names(data))
#   stopifnot(W %in% names(data))
#   assertthat::assert_that(all(data[[W]] %in% 0:1))
#   assertthat::assert_that(length(unique(data[[W]])) > 1)
#   
#   if(is.null(b)){
#     b <- round(n^0.8)
#   }
#   part_idx <- seq_len(b)
#   n <- nrow(data)
# 
#   # Change to factor variable and change names to make compatible with caret
#   data$Tr2 <- factor(data[[W]])
#   levels(data$Tr2) <- make.names(levels(data$Tr2))
#   # Update prop_formula
#   prop_formula <- update.formula(as.formula(prop_formula), Tr2 ~ .)
#   n1 <- sum(data[[W]] == 1)
#   n0 <- n - n1
#   
#   partition <- make_partition(n = n, subsets = subsets, b = b, disjoint = disjoint)
#   if(hp_tune == TRUE){
#     # Hyperparameter tuning on a larger set
#     hyp_p <- make_partition(n = n, subsets = 1, b = b)[[1]]
#     if(b < n){
#       # Use OoB as holdout set
#       prop_control$index <- list(as.integer(hyp_p))
#       y_control$index <- list(as.integer(hyp_p))
#     }
# 
#     # Tune hyperparameters
#     m <- caret::train(form = as.formula(y_formula),
#                       data = data,
#                       trControl = y_control,
#                       method = y_method,
#                       tuneGrid = y_grid,
#                       verbose = FALSE,
#                       ...)
#     
#     prop_mod <- caret::train(form = as.formula(prop_formula),
#                              data = data,
#                              trControl = prop_control,
#                              method = prop_method,
#                              metric = 'logLoss',
#                              tuneGrid = prop_grid,
#                              verbose = FALSE,
#                              ...)
#     
#     m_tune <- m$bestTune
#     prop_tune <- prop_mod$bestTune
#   } else{
#     # assertthat::assert_that(!is.null(prop_grid))
#     # assertthat::assert_that(!is.null(y_grid))
#     # assertthat::assert_that(nrow(y_grid) == 1)
#     # assertthat::assert_that(nrow(prop_grid) == 1)
#     m_tune <- y_grid
#     prop_tune <- prop_grid
#     # Use defaults
#     ## GBM
#     # m_tune <- data.frame(mtry = 2, min.node.size = 10, max.depth = 10)
#     # prop_tune <- data.frame(mtry = 1, min.node.size = 20, max.depth = 10)
#     # m_tune <- data.frame(n.trees = 500, interaction.depth = 1, shrinkage = 0.1, n.minobsinnode = 5)
#     # prop_tune <- data.frame(n.trees = 500, interaction.depth = 1, shrinkage = 0.1, n.minobsinnode = 10)
# 
#   }
#   
#   part_idx <- seq_len(b)
#   
#   blb_estim <- lapply(partition, function(p){
#     process_partition(partition_element = p,
#                       part_idx = part_idx,
#                       b = b,
#                       n = n,
#                       data = data, 
#                       y_method = y_method, 
#                       prop_method = prop_method, 
#                       K = K, 
#                       y_formula = y_formula, 
#                       prop_formula = prop_formula,
#                       m_tune = m_tune, 
#                       prop_tune = prop_tune, 
#                       W = W, 
#                       y_args = y_args, prop_args = prop_args, 
#                       boot_stat = dr, 
#                       r = r, Y = Y,
#                       crossfit = crossfit)
#    
#   })
# 
#   blb_estim <- rbindlist(blb_estim)
#   if(aggregate){
#     blb_estim <- blb_estim[, .(estim = mean(estim),
#                                estim_sub = mean(estim_sub),
#                                sd_sub = mean(sd_sub),
#                                sdthetahat = mean(sdthetahat),
#                                lower = mean(lower),
#                                upper = mean(upper)), by = 'boot_type']
#     
#   }
#   blb_estim
# }
# 
# 
dr <- function(data, Tr, y, return_var = FALSE){
  phi1 <- (data[[Tr]]/data$prop_score)*(data[[y]]-data$m1) + data$m1
  phi0 <- ((1-data[[Tr]])/(1-data$prop_score))*(data[[y]]-data$m0) + data$m0
  if(return_var){
    return(c(mean(phi1) - mean(phi0), var(phi1 - phi0)/length(phi1)))
  } else{
    return(mean(phi1) - mean(phi0))
  }
}
# 
# pc <- function(data, Tr, y){
#   ate <- dr(data = data, Tr = Tr, y = y)
#   return(ate/mean(data$m0))
# }
# 
jackknife_acceleration <- function(data, ...){
  nd <- nrow(data)
  thetas <- sapply(seq_len(nd), function(i){
    dt <- data[-i]
    dr(d = dt, ...)
  })
  mean_theta <- mean(thetas)
  L <- mean_theta - thetas
  sum(L^3)/(6*sum(L^2)^(3/2))
}