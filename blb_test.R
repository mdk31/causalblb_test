library(data.table)
library(pbapply)

kangschafer3 <- function(n, te, beta_overlap = 0.5, sigma) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

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


source('~/Documents/HW/Research/CI/cblb_dml_dr/paper/code/dbml_funcs.R')

boot_types <- c('percentile', 'bca', 'basic', 'normal', 'studentized', 'asymptotic')
te <- 0.8
sigma <- 1
replications <- 500
r <- 100

setwd('~/Documents/HW/Research/CI/causalblb_test')
base_nm <- 'simple_test'

image_path <- 'images'
dat_path <- 'data'

# Create the temporary directory
temp_dir <- file.path(dat_path, paste0(base_nm, '_tmp'))
if(!file.exists(temp_dir)){
  dir.create(temp_dir, recursive = TRUE)
}

img_tmp_dir <- file.path(image_path, paste0(base_nm, '_tmp'))
if(!file.exists(img_tmp_dir)){
  dir.create(img_tmp_dir, recursive = TRUE)
}

hyper_grid <- as.data.table(expand.grid(n = c(20000),
                                        gamma = c(0.7, 1),
                                        subsets = c(50),
                                        prop_form = c('misspecified'),
                                        stringsAsFactors = FALSE))
hyper_grid[gamma == 1, `:=`(subsets = 1)]
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
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        # m <- glm(y ~ Tr + X1 + X2, data = tmp_dat, family = 'gaussian')
        # g <- glm(as.formula(form), data = tmp_dat, family = 'binomial')
        # 
        # prop_score <- predict(g, type = 'response')
        # m1 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 1))
        # m0 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 0))
        # 
        # tmp_dat$prop_score <- prop_score
        # tmp_dat$m1 <- m1
        # tmp_dat$m0 <- m0

        blb_reps <- replicate(r, {
          boot_idx <- sample(part_idx, size = n, replace = TRUE)
          boot_dat <- tmp_dat[boot_idx]
          m <- glm(y ~ Tr + X1 + X2, data = boot_dat, family = 'gaussian')
          g <- glm(as.formula(form), data = boot_dat, family = 'binomial')
          
          prop_score <- predict(g, type = 'response')
          m1 <- predict(m, data.frame(X1 = boot_dat$X1, X2 = boot_dat$X2, Tr = 1))
          m0 <- predict(m, data.frame(X1 = boot_dat$X1, X2 = boot_dat$X2, Tr = 0))
          
          boot_dat$prop_score <- prop_score
          boot_dat$m1 <- m1
          boot_dat$m0 <- m0
          phi1 <- (boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1
          phi0 <- (1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0
          mean(phi1) - mean(phi0)
        })
        
        ci <- boot:::perc.ci(t = blb_reps)[, 4:5]
        
        data.frame(boot_type = 'percentile',
                   lower = ci[1],
                   upper = ci[2])
      })
      
      blb_out <- rbindlist(blb_out)
      blb_out <- blb_out[, .(lower = mean(lower), upper = mean(upper)), by = 'boot_type']
      blb_out
    }, cl = 4)
    
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


cblb[, .(coverage = mean(lower <= te & upper >= te)), by = c('boot_type', 'n', 'gamma', 'subsets', 'prop_form')]

