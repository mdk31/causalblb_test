library(data.table)
library(pbapply)
library(ggplot2)

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

base_nm <- 'parameter_bias'

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

hyper_grid <- as.data.table(expand.grid(n = c(10000, 20000, 100000),
                                        gamma = c(0.7, 0.75, 0.8)))
seq_row <- seq_len(nrow(hyper_grid))

# CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'simulations.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'simulations.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- round(n^gamma)
    part_idx <- seq_len(b)

    out <- lapply(seq_len(1), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      partitions <- make_partition(n = n, subsets = 500, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- pblapply(partitions, function(i){
        tmp_dat <- dat[i]
        m <- glm(y ~ Tr + X1 + X2, data = tmp_dat, family = 'gaussian')
        data.frame(coef = coef(m)['Tr'])
      }, cl = 4)
      
      blb_out <- rbindlist(blb_out)
      blb_out
    })
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'simulations.rds'))
}

ggplot(cblb, aes(x = coef)) +
  geom_histogram() +
  facet_wrap(n ~ gamma)

