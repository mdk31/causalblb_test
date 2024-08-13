library(data.table)
library(pbapply)
library(xgboost)

base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100
K <- 10

base_nm <- 'blb_standard_error_v1_ml'
# Full sample + BLB standard error

image_path <- 'images'
dat_path <- 'data'

# Create the temporary directory
temp_dir <- file.path(base_dir, dat_path, paste0(base_nm, '_tmp'))
if(!file.exists(temp_dir)){
  dir.create(temp_dir, recursive = TRUE)
}

img_tmp_dir <- file.path(base_dir, image_path, paste0(base_nm, '_tmp'))
if(!file.exists(img_tmp_dir)){
  dir.create(img_tmp_dir, recursive = TRUE)
}


hyper_grid <- as.data.table(expand.grid(n = c(50000),
                        gamma = c(0.7),
                        subsets = c(10, 50),
                        prop_form = c('correct', 'wrong'),
                        out_form = c('correct', 'wrong')))

seq_row <- seq_len(nrow(hyper_grid))

# CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'coverage.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- round(n^gamma)
    part_idx <- seq_len(b)
    prop_form <- grid_val$prop_form
    out_form <- grid_val$out_form
    if(prop_form == 'correct'){
      prop_formula <- c('X1', 'X2')
    } else{
      prop_formula <- c('Z1', 'Z2')
    }
    if(out_form == 'correct'){
      out_formula <- c('Tr', 'X1', 'X2')
    } else{
      out_formula <- c('Tr', 'Z1', 'Z2')
    }
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      # Full sample
      m_train <- as.matrix(dat[, out_formula, with = FALSE])
      g_train <- as.matrix(dat[, prop_formula, with = FALSE])
      
      m <- xgboost(data = m_train, label = dat$y, verbose = 0, nrounds = 10,
                   params = list(objective = 'reg:squarederror'))
      g <- xgboost(data = g_train, label = dat$Tr, verbose = 0, nrounds = 10, 
                   params = list(objective = 'binary:logistic'))

      prop_score <- predict(g, newdata = g_train)
      
      newdata <- data.frame(Tr = 1, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
      newdata <- as.matrix(newdata[, out_formula])
      m1 <- predict(m, newdata = newdata)
      
      newdata <- data.frame(Tr = 0, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
      newdata <- as.matrix(newdata[, out_formula])
      m0 <- predict(m, newdata = newdata)
      
      full_dat <- copy(dat)
      full_dat$prop_score <- prop_score
      full_dat$m1 <- m1
      full_dat$m0 <- m0
      
      phi1_full <- (full_dat$Tr/full_dat$prop_score)*(full_dat$y - full_dat$m1) + full_dat$m1
      phi0_full <- (1 - full_dat$Tr)/(1 - full_dat$prop_score)*(full_dat$y - full_dat$m0) + full_dat$m0
      tau_hat_full <- mean(phi1_full) - mean(phi0_full)

      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        m_train <- as.matrix(tmp_dat[, out_formula, with = FALSE])
        g_train <- as.matrix(tmp_dat[, prop_formula, with = FALSE])
        
        M <- rmultinom(n = r, size = n, prob = rep(1, b))

        blb_reps <- sapply(seq_len(r), function(bt){
          m <- xgboost(data = m_train, label = tmp_dat$y, verbose = 0, nrounds = 10,
                       params = list(objective = 'reg:squarederror'), weight = M[, bt])
          g <- xgboost(data = g_train, label = tmp_dat$Tr, verbose = 0, nrounds = 10, 
                       params = list(objective = 'binary:logistic'), weight = M[, bt])
          
          prop_score <- predict(g, newdata = g_train)
          
          newdata <- data.frame(Tr = 1, X1 = tmp_dat$X1, X2 = tmp_dat$X2, Z1 = tmp_dat$Z1, Z2 = tmp_dat$Z2)
          newdata <- as.matrix(newdata[, out_formula])
          m1 <- predict(m, newdata = newdata)
          
          newdata <- data.frame(Tr = 0, X1 = tmp_dat$X1, X2 = tmp_dat$X2, Z1 = tmp_dat$Z1, Z2 = tmp_dat$Z2)
          newdata <- as.matrix(newdata[, out_formula])
          m0 <- predict(m, newdata = newdata)
          
          tmp_dat$prop_score <- prop_score
          tmp_dat$m1 <- m1
          tmp_dat$m0 <- m0
          
          phi1 <- M[, bt]*((tmp_dat$Tr/tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m1) + tmp_dat$m1)
          phi0 <- M[, bt]*((1 - tmp_dat$Tr)/(1 - tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m0) + tmp_dat$m0)
          sum(phi1)/n - sum(phi0)/n
        })

        sd_sub <- sd(blb_reps)
        return(data.table(sd_sub = sd_sub))
      })
      
      
      blb_out <- rbindlist(blb_out)
      blb_out <- blb_out[, .(boot_sd = mean(sd_sub),
                             tau_hat = tau_hat_full)]
      blb_out
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               prop_form = prop_form,
               out_form = out_form)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

cblb[, `:=`(lower_ci = tau_hat - qnorm(0.975)*boot_sd,
            upper_ci = tau_hat + qnorm(0.975)*boot_sd)][
              , .(coverage = mean(lower_ci <= te & upper_ci  >= te)), by = c('n', 'subsets', 'prop_form', 'out_form')]
