library(data.table)
library(pbapply)
library(caret)
base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100

base_nm <- 'penalization'

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

# hyper_grid <- as.data.table(expand.grid(n = c(20000, 5000000),
#                                         gamma = c(0.8),
#                                         subsets = c(25),
#                                         prop_form = c('correct'),
#                                         stringsAsFactors = FALSE))
ns <- c(10000)
hyper_grid <- as.data.table(expand.grid(n = ns,
                        gamma = c(0.7),
                        subsets = c(10, 20),
                        lambda = c(0, 0.001, 0.01),
                        prop_form = c('correct'), stringsAsFactors = FALSE))

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
    if(prop_form == 'correct'){
      form <- Tr ~ X1 + X2
    } else{
      form <- Tr ~ 1
      
    }
    subsets <- grid_val$subsets
    lambda <- grid_val$lambda
    lambda_grid <- expand.grid(alpha = 0, lambda = lambda)
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      dat$Tr2 <- factor(dat$Tr)
      levels(dat$Tr2) <- make.names(levels(dat$Tr2))
      # Update prop_formula
      prop_formula <- update.formula(as.formula(form), Tr2 ~ .)
      # Full sample
      
      # Train the ridge regression model without cross-validation
      m <- train(y ~ Tr + X1 + X2, 
                 data = dat,
                 family = 'gaussian',
                 method = "glmnet",
                 trControl = trainControl(method = "none"), 
                 tuneGrid = lambda_grid)
      
      g <- train(as.formula(prop_formula), 
                 data = dat,
                 family = 'binomial',
                 method = "glmnet",
                 trControl = trainControl(method = "none"), 
                 tuneGrid = lambda_grid)
      
      prop_score <- predict(g, type = 'prob')[, 2]
      m1 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 1))
      m0 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 0))
      
      full_dat <- copy(dat)
      full_dat$prop_score <- prop_score
      full_dat$m1 <- m1
      full_dat$m0 <- m0
      
      phi1_full <- (full_dat$Tr/full_dat$prop_score)*(full_dat$y - full_dat$m1) + full_dat$m1
      phi0_full <- (1 - full_dat$Tr)/(1 - full_dat$prop_score)*(full_dat$y - full_dat$m0) + full_dat$m0
      
      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        m <- train(y ~ Tr + X1 + X2, 
                   data = tmp_dat,
                   family = 'gaussian',
                   method = "glmnet",
                   trControl = trainControl(method = "none"), 
                   tuneGrid = lambda_grid)
        
        g <- train(as.formula(prop_formula), 
                   data = tmp_dat,
                   family = 'binomial',
                   method = "glmnet",
                   trControl = trainControl(method = "none"), 
                   tuneGrid = lambda_grid)
        
        prop_score <- predict(g, type = 'prob')[, 2]
        m1 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 1))
        m0 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 0))
        
        tmp_dat$prop_score <- prop_score
        tmp_dat$m1 <- m1
        tmp_dat$m0 <- m0
        
        phi1_subset <- (tmp_dat$Tr/tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m1) + tmp_dat$m1
        phi0_subset <- (1 - tmp_dat$Tr)/(1 - tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m0) + tmp_dat$m0

        blb_reps <- replicate(r, {
          boot_idx <- sample(part_idx, size = n, replace = TRUE)
          boot_dat <- tmp_dat[boot_idx]
          phi1 <- (boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1
          phi0 <- (1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0
          mean(phi1) - mean(phi0)
        })
        

        data.frame(estim_subset = mean(phi1_subset) - mean(phi0_subset),
                   sd_subset = var(phi1_subset - phi0_subset)/length(phi1_subset),
                   boot_reps = blb_reps,
                   boot_rep_num = seq_len(r))
      })
      
      blb_out <- rbindlist(blb_out)
      blb_out[, `:=`(subset_num = rep(seq_along(partitions), each = r),
                     replication = rp,
                     estim_full = mean(phi1_full) - mean(phi0_full),
                     sd_full = var(phi1_full - phi0_full)/length(phi1_full))]
      blb_out
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               lambda = lambda,
               subsets = subsets,
               prop_form = prop_form)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

tmp <- cblb[, .(lower = boot:::perc.ci(boot_reps)[4], upper = boot:::perc.ci(boot_reps)[5]), 
           by = c('subset_num', 'replication', 'n', 'subsets', 'lambda')][
             ,.(lower = mean(lower), upper = mean(upper)), by = c('replication', 'n', 'subsets', 'lambda')
           ]
tmp[, .(coverage = mean(lower <= te & upper >= te)), by = c('n', 'subsets', 'lambda')]

