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

base_nm <- 'blb_standard_error_v4_ml_full'
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
                        prop_form = c('correct', 'wrong'),
                        out_form = c('correct', 'wrong')))

idx <- seq_len(50000)

seq_row <- seq_len(nrow(hyper_grid))

# CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'coverage.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
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
      
      blb_reps <- replicate(r, function(x){
        boot_idx <- sample(idx, size = n, replace = TRUE)
        m_train_boot <- m_train[boot_idx, ]
        g_train_boot <- g_train[boot_idx, ]
        boot_dat <- dat[boot_idx]
        
        m <- xgboost(data = m_train_boot, label = boot_dat$y, verbose = 0, nrounds = 10,
                     params = list(objective = 'reg:squarederror'))
        g <- xgboost(data = g_train_boot, label = boot_dat$Tr, verbose = 0, nrounds = 10, 
                     params = list(objective = 'binary:logistic'))
        
        prop_score <- predict(g, newdata = g_train_boot)
        
        newdata <- data.frame(Tr = 1, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
        newdata <- as.matrix(newdata[, out_formula])
        m1 <- predict(m, newdata = newdata)
        
        newdata <- data.frame(Tr = 0, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = tmp_dat$Z2)
        newdata <- as.matrix(newdata[, out_formula])
        m0 <- predict(m, newdata = newdata)
        
        boot_dat$prop_score <- prop_score
        boot_dat$m1 <- m1
        boot_dat$m0 <- m0
        
        phi1 <- ((boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1)
        phi0 <- ((1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0)
        sum(phi1)/n - sum(phi0)/n
      })

      
      perc_ci <- boot:::perc.ci(blb_reps)
      blb_out <- data.table(lower_ci = perc_ci[4],
                        upper_ci = perc_ci[5])
      blb_out
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               prop_form = prop_form,
               out_form = out_form)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

cblb[, .(coverage = mean(lower_ci <= te & upper_ci  >= te)), by = c('n', 'subsets', 'prop_form', 'out_form')]
