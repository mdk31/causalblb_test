library(data.table)
library(ggplot2)
library(pbapply)

base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 10
r <- 100

base_nm <- 'convergence'
# Try BCa intervals, using subset estimator to estimate acceleration

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


n <- 100000
idx <- seq_len(n)
hyper_grid <- as.data.table(expand.grid(n = c(n),
                        gamma = c(0.5, 0.6, 0.7, 0.8, 0.9),
                        bias_corrected = c(TRUE, FALSE),
                        subsets = c(500)))

seq_row <- seq_len(nrow(hyper_grid))

# Ground truth----

ground_truth <- pblapply(seq_len(replications), function(rp){
  set.seed(rp)
  dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
  
  boot_reps <- replicate(1000, {
    boot_dat <- dat[sample(idx, size = n, replace = TRUE)]
    m <- glm(y ~ Tr + X1 + X2, data = boot_dat, family = 'gaussian')
    g <- glm(Tr ~ X1 + X2, data = boot_dat, family = 'binomial')
    
    prop_score <- predict(g, type = 'response')
    m1 <- predict(m, data.frame(X1 = boot_dat$X1, X2 = boot_dat$X2, Tr = 1))
    m0 <- predict(m, data.frame(X1 = boot_dat$X1, X2 = boot_dat$X2, Tr = 0))

    full_dat <- copy(boot_dat)
    full_dat$prop_score <- prop_score
    full_dat$m1 <- m1
    full_dat$m0 <- m0
    
    phi1_full <- (full_dat$Tr/full_dat$prop_score)*(full_dat$y - full_dat$m1) + full_dat$m1
    phi0_full <- (1 - full_dat$Tr)/(1 - full_dat$prop_score)*(full_dat$y - full_dat$m0) + full_dat$m0
    
    mean(phi1_full) - mean(phi0_full)
  })
  
  perc_ci <- boot:::perc.ci(boot_reps)
  return(data.table(true_lower = perc_ci[4], true_upper = perc_ci[5], rep = rp))
}, cl = 4)

ground_truth <- rbindlist(ground_truth)

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
    subsets <- grid_val$subsets
    bias_correct <- grid_val$bias_corrected
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      
      # Full sample
      m <- glm(y ~ Tr + X1 + X2, data = dat, family = 'gaussian')
      g <- glm(Tr ~ X1 + X2, data = dat, family = 'binomial')

      prop_score <- predict(g, type = 'response')
      m1 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 1))
      m0 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 0))
      #
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
        m <- glm(y ~ Tr + X1 + X2, data = tmp_dat, family = 'gaussian')
        g <- glm(Tr ~ X1 + X2, data = tmp_dat, family = 'binomial')
        
        prop_score <- predict(g, type = 'response')
        m1 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 1))
        m0 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 0))
        #
        tmp_dat$prop_score <- prop_score
        tmp_dat$m1 <- m1
        tmp_dat$m0 <- m0
        
        
        M <- rmultinom(n = r, size = n, prob = rep(1, b))
        blb_reps <- sapply(seq_len(r), function(bt){
          m_boot <- glm(y ~ Tr + X1 + X2, data = tmp_dat, family = 'gaussian', weights = M[, bt])
          g_boot <- glm(Tr ~ X1 + X2, data = tmp_dat, family = 'binomial', weights = M[, bt])
          
          prop_score <- predict(g_boot, newdata = tmp_dat, type = 'response')
          m1 <- predict(m_boot, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 1))
          m0 <- predict(m_boot, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 0))
          
          tmp_dat$prop_score <- prop_score
          tmp_dat$m1 <- m1
          tmp_dat$m0 <- m0

          phi1 <- M[, bt]*((tmp_dat$Tr/tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m1) + tmp_dat$m1)
          phi0 <- M[, bt]*((1 - tmp_dat$Tr)/(1 - tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m0) + tmp_dat$m0)
          sum(phi1)/n - sum(phi0)/n
        })
        
        if(!bias_correct){
          blb_ci <- boot:::perc.ci(blb_reps)
          return(data.table(blb_lower = blb_ci[4],
                            blb_upper = blb_ci[5]))
        } else{
          bias <- mean(blb_reps) - tau_hat_full
          blb_bias_correct <- blb_reps - bias
          blb_ci <- boot:::perc.ci(blb_bias_correct)

          return(data.table(blb_lower = blb_ci[4],
                            blb_upper = blb_ci[5]))
        }
      })
      
      
      blb_out <- rbindlist(blb_out)
      blb_out[, `:=`(rep = rp, subset_num = seq_len(subsets))]
      blb_out
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               bias_correct = bias_correct)]
    out
  })
  cblb <- rbindlist(cblb)
  cblb <- merge(cblb, ground_truth, on = 'rep', all.x = TRUE)
  
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

setorder(cblb, rep, subset_num)
cblb[, `:=`(blb_lower = cumsum(blb_lower)/subset_num, blb_upper = cumsum(blb_upper)/subset_num), 
             by = c('rep', 'n', 'gamma', 'subsets', 'bias_correct')]
cblb[, `:=`(error = abs(true_lower - blb_lower)/2 + abs(true_upper - blb_upper)/2)]
trajectory <- cblb[, .(mean_error = mean(error)), by = c('n', 'gamma', 'subsets', 'bias_correct', 'subset_num')]

trajectory[, `:=`(identifier = paste0('gamma: ', gamma, ' bias_correct: ', bias_correct))]


ggplot(trajectory[gamma != 0.5], aes(x = subset_num, y = mean_error, color = factor(gamma))) + 
  geom_line() +
  facet_grid(~ bias_correct) +
  theme_minimal()

ggplot(cblb[rep==1 & gamma != 0.5], aes(x = subset_num, y = blb_lower)) + 
  geom_line() +
  geom_hline(aes(yintercept= true_lower)) +
  geom_hline(aes(yintercept = true_upper)) +
  facet_wrap(~ gamma + bias_correct)
# tmp <- cblb[, ]
# tmp <- cblb[, .(lower = boot:::perc.ci(boot_reps)[4], upper = boot:::perc.ci(boot_reps)[5]), 
#            by = c('subset_num', 'replication', 'n', 'subsets', 'prop_form', 'out_form')][
#              ,.(lower = mean(lower), upper = mean(upper)), by = c('replication', 'n', 'subsets', 'prop_form', 'out_form')
#            ]
# cblb[, .(coverage = mean(lower_ci <= te & upper_ci  >= te)), by = c('n', 'subsets', 'prop_form', 'out_form')]
