library(data.table)
library(pbapply)
library(ranger)
base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100
K <- 10

base_nm <- 'random_forest'

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
# hyper_grid <- data.table(n = ns,
#                         gamma = c(rep(0.7, 1), 1.0),
#                         subsets = c(20, 40, 60, 1),
#                         prop_form = c('correct'))

hyper_grid <- as.data.table(expand.grid(n = ns,
                                        gamma = c(0.7),
                                        subsets = c(10, 20),
                                        prop_form = c('correct'), stringsAsFactors = FALSE))

hyper_grid <- rbindlist(list(hyper_grid, data.table(n = ns, gamma = 1, subsets = 1, prop_form = 'correct')))

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
    cost <- grid_val$cost
    out <- lapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      # Full sample
      m <- ranger(y ~ Tr + X1 + X2, data = dat)
      g <- ranger(as.formula(form), data = dat, probability = TRUE)
      

      prop_score <- predict(g, data = dat)$predictions[, 2]
      m1 <- predict(m, data = data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 1))$predictions
      m0 <- predict(m, data = data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 0))$predictions
      
      full_dat <- copy(dat)
      full_dat$prop_score <- prop_score
      full_dat$m1 <- m1
      full_dat$m0 <- m0
      
      phi1_full <- (full_dat$Tr/full_dat$prop_score)*(full_dat$y - full_dat$m1) + full_dat$m1
      phi0_full <- (1 - full_dat$Tr)/(1 - full_dat$prop_score)*(full_dat$y - full_dat$m0) + full_dat$m0
      
      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        folds <- split(part_idx, sample(rep(1:K, length.out = length(part_idx))))
        cross_dat <- lapply(folds, function(test_idx){
          train_idx <- setdiff(part_idx, test_idx)
          train_dat <- tmp_dat[train_idx]
          test_dat <- tmp_dat[-train_idx]
          
          m <- ranger(y ~ Tr + X1 + X2, data = train_dat)
          g <- ranger(as.formula(form), data = train_dat, probability = TRUE,
                      min.node.size = 25)
          
          prop_score <- predict(g, data = test_dat)$predictions[, 2]
          m1 <- predict(m, data = data.frame(X1 = test_dat$X1, X2 = test_dat$X2, Tr = 1))$predictions
          m0 <- predict(m, data = data.frame(X1 = test_dat$X1, X2 = test_dat$X2, Tr = 0))$predictions
          
          test_dat$prop_score <- prop_score
          test_dat$m1 <- m1
          test_dat$m0 <- m0
          
          return(test_dat)
          
        })
  
        cross_dat <- rbindlist(cross_dat)
        cross_dat <- cross_dat[between(prop_score, 0.001, 0.999)]
        phi1_subset <- (cross_dat$Tr/cross_dat$prop_score)*(cross_dat$y - cross_dat$m1) + cross_dat$m1
        phi0_subset <- (1 - cross_dat$Tr)/(1 - cross_dat$prop_score)*(cross_dat$y - cross_dat$m0) + cross_dat$m0

        blb_reps <- replicate(r, {
          boot_idx <- sample(part_idx, size = n, replace = TRUE)
          boot_dat <- cross_dat[boot_idx]
          phi1 <- (boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1
          phi0 <- (1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0
          mean(phi1) - mean(phi0)
        })
        
        assertthat::assert_that(all(!is.na(blb_reps)) & all(!is.infinite(blb_reps)))
        

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
    })

    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               prop_form = prop_form)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

cblb <- cblb[!is.na(boot_reps) & !is.infinite(boot_reps)]

# Bias
biases <- cblb[subsets == 1, .(bias = mean(boot_reps) - te), by = c('n', 'replication')]
summary(biases$bias)

tmp <- cblb[, .(lower = boot:::perc.ci(boot_reps)[4], upper = boot:::perc.ci(boot_reps)[5]), 
           by = c('subset_num', 'replication', 'n', 'subsets')][
             ,.(lower = mean(lower), upper = mean(upper)), by = c('replication', 'n', 'subsets')
           ]
tmp[, .(coverage = mean(lower <= te & upper >= te)), by = c('n', 'subsets')]

