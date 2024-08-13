library(data.table)
library(pbapply)

base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100

base_nm <- 'weighted_bootstrap_learn_subsetpred'

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
      form <- 'Tr ~ X1 + X2'
    } else{
      form <- 'Tr ~ 1'
    }
    if(out_form == 'correct'){
      out_formula <- 'y ~ X1 + X2'
    } else{
      out_formula <- 'y ~ 1'
    }
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      n1 <- sum(dat$Tr == 1)
      n0 <- sum(dat$Tr == 0)
      # Full sample
      # m <- glm(y ~ Tr + X1 + X2, data = dat, family = 'gaussian')
      # g <- glm(as.formula(form), data = dat, family = 'binomial')
      # 
      # prop_score <- predict(g, type = 'response')
      # m1 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 1))
      # m0 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 0))
      # 
      # full_dat <- copy(dat)
      # full_dat$prop_score <- prop_score
      # full_dat$m1 <- m1
      # full_dat$m0 <- m0
      # 
      # phi1_full <- (full_dat$Tr/full_dat$prop_score)*(full_dat$y - full_dat$m1) + full_dat$m1
      # phi0_full <- (1 - full_dat$Tr)/(1 - full_dat$prop_score)*(full_dat$y - full_dat$m0) + full_dat$m0
      
      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        # m <- glm(y ~ Tr + X1 + X2, data = tmp_dat, family = 'gaussian')
        g <- glm(as.formula(form), data = tmp_dat, family = 'binomial')

        prop_score <- predict(g, type = 'response')
        wts <- (tmp_dat$Tr == 1)*(1/prop_score) + (tmp_dat$Tr == 0)*(1/(1 - prop_score))
        ns <- abs(tapply(wts, tmp_dat$Tr, sum))
        ns <- ns['1']*(tmp_dat$Tr == 1) + (1 - (tmp_dat$Tr == 1))*ns['0']
        wts <- wts/ns
        
        wts <- split(wts, tmp_dat$Tr)
        trt_idx <- which(tmp_dat$Tr == 1)
        crt_idx <- which(tmp_dat$Tr == 0)
        

        blb_reps <- replicate(r, {
          browser()
          trt_boot <- sample(trt_idx, size = n1, replace = TRUE, prob = wts[['1']])
          crt_boot <- sample(crt_idx, size = n0, replace = TRUE, prob = wts[['0']])
          
          trt_boot <- tmp_dat[trt_boot]
          crt_boot <- tmp_dat[crt_boot]
          m1 <- glm(y ~ X1 + X2, data = trt_boot, family = 'gaussian')
          m0 <- glm(y ~ X1 + X2, data = crt_boot, family = 'gaussian')
          
          m1 <- predict(m1, trt_boot)
          m0 <- predict(m0, crt_boot)
          mean(m1) - mean(m0)
        })
        

        data.frame(boot_reps = blb_reps,
                   boot_rep_num = seq_len(r))
      })
      
      blb_out <- rbindlist(blb_out)
      blb_out[, `:=`(subset_num = rep(seq_along(partitions), each = r),
                     replication = rp)]
      blb_out
    }, cl = 1)
    
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

tmp <- cblb[, .(lower = boot:::perc.ci(boot_reps)[4], upper = boot:::perc.ci(boot_reps)[5]), 
           by = c('subset_num', 'replication', 'n', 'subsets')][
             ,.(lower = mean(lower), upper = mean(upper)), by = c('replication', 'n', 'subsets')
           ]
tmp[, .(coverage = mean(lower <= te & upper >= te)), by = c('n', 'subsets')]

