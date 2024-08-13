library(data.table)
library(pbapply)

base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100
K <- 10

base_nm <- 'parametric_weighted_highergamma'

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


# ns <- c(10000)
hyper_grid <- as.data.table(expand.grid(n = c(50000),
                        gamma = c(0.8),
                        subsets = c(30),
                        prop_form = c('correct', 'wrong'),
                        out_form = c('correct', 'wrong')))

# hyper_grid <- hyper_grid[!(subsets == 1 & gamma != 1) & !(gamma == 1 & subsets != 1)]
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
      prop_formula <- 'Tr ~ X1 + X2'
    } else{
      prop_formula <- 'Tr ~ 1'
    }
    if(out_form == 'correct'){
      out_formula <- 'y ~ Tr + X1 + X2'
    } else{
      out_formula <- 'y ~ 1'
    }
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      # Full sample
      m <- glm(y ~ Tr + X1 + X2, data = dat, family = 'gaussian')
      g <- glm(as.formula(prop_formula), data = dat, family = 'binomial')
      
      prop_score <- predict(g, type = 'response')
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
        m <- glm(as.formula(out_formula), data = tmp_dat, family = 'gaussian')
        g <- glm(as.formula(prop_formula), data = tmp_dat, family = 'binomial')

        prop_score <- predict(g, type = 'response')
        m1 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 1))
        m0 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 0))
        
        tmp_dat$prop_score <- prop_score
        tmp_dat$m1 <- m1
        tmp_dat$m0 <- m0
        
        phi1_subset <- (tmp_dat$Tr/tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m1) + tmp_dat$m1
        phi0_subset <- (1 - tmp_dat$Tr)/(1 - tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m0) + tmp_dat$m0
        
        # folds <- split(part_idx, sample(rep(1:K, length.out = length(part_idx))))
        M <- rmultinom(n = r, size = n, prob = rep(1, b))
        
        blb_reps <- sapply(seq_len(r), function(bt){
          m_boot <- glm(as.formula(out_formula), data = tmp_dat, family = 'gaussian', weights = M[, bt])
          g_boot <- glm(as.formula(prop_formula), data = tmp_dat, family = 'binomial', weights = M[, bt])
          
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
        

        # blb_reps <- replicate(r, {
        #   boot_idx <- sample(part_idx, size = n, replace = TRUE)
        #   boot_dat <- tmp_dat[boot_idx]
        # 
        #   boot_dat_m <- aggregate(y ~ Tr + X1 + X2, boot_dat, FUN = function(x) c(count = length(x), mean = mean(x)))
        #   boot_dat_m <- data.frame(X1 = boot_dat_m$X1, X2 = boot_dat_m$X2, Tr = boot_dat_m$Tr, 
        #                            y = boot_dat_m$y[, "mean"], weights = boot_dat_m$y[, "count"])
        #   
        #   boot_dat_g <- aggregate(Tr ~ X1 + X2, boot_dat, FUN = function(x) c(count = length(x), mean = mean(x)))
        #   boot_dat_g <- data.frame(X1 = boot_dat_g$X1, X2 = boot_dat_g$X2,
        #                            Tr = boot_dat_g$Tr[, "mean"], weights = boot_dat_g$Tr[, "count"])
        #   
        #   m_boot <- glm(as.formula(out_formula), data = boot_dat_m, family = 'gaussian', weights = weights)
        #   g_boot <- glm(as.formula(prop_formula), data = boot_dat_g, family = 'binomial', weights = weights)
        #   
        #   prop_score <- predict(g_boot, newdata = boot_dat, type = 'response')
        #   m1 <- predict(m, data.frame(X1 = boot_dat$X1, X2 = boot_dat$X2, Tr = 1))
        #   m0 <- predict(m, data.frame(X1 = boot_dat$X1, X2 = boot_dat$X2, Tr = 0))
        #   
        #   boot_dat$prop_score <- prop_score
        #   boot_dat$m1 <- m1
        #   boot_dat$m0 <- m0
        # 
        #   phi1 <- (boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1
        #   phi0 <- (1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0
        #   mean(phi1) - mean(phi0)
        # })
        

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
    }, cl = 5)
    
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

tmp <- cblb[, .(lower = boot:::perc.ci(boot_reps)[4], upper = boot:::perc.ci(boot_reps)[5]), 
           by = c('subset_num', 'replication', 'n', 'subsets', 'prop_form', 'out_form')][
             ,.(lower = mean(lower), upper = mean(upper)), by = c('replication', 'n', 'subsets', 'prop_form', 'out_form')
           ]
tmp[, .(coverage = mean(lower <= te & upper >= te)), by = c('n', 'subsets', 'prop_form', 'out_form')]

