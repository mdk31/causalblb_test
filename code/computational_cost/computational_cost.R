library(data.table)
library(pbapply)
library(ggplot2)

base_dir <- '~/Documents/HW/Research/CI/causalblb_test'

source(file.path(base_dir, 'code/helper_funcs.R'))

te <- 0.8
sigma <- 1
replications <- 25
r <- 100

base_nm <- 'computational_cost'

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

n <- c(20000)
hyper_grid <- as.data.table(expand.grid(n = n,
                          gamma = c(0.7),
                          subsets = c(60, 90, 200),
                          num_cores = c(2, 4, 6),
                          prop_form = c('correct')))

seq_row <- seq_len(nrow(hyper_grid))
set.seed(123)
dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
idx <- seq_len(n)
# Full Bootstrap Timing

full <- lapply(seq_len(replications), function(i){
  times <- system.time({
    m <- glm(y ~ Tr + X1 + X2, data = dat, family = 'gaussian')
    g <- glm(Tr ~ X1 + X2, data = dat, family = 'binomial')
    
    prop_score <- predict(g, type = 'response')
    m1 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 1))
    m0 <- predict(m, data.frame(X1 = dat$X1, X2 = dat$X2, Tr = 0))
    
    full_dat <- copy(dat)
    full_dat$prop_score <- prop_score
    full_dat$m1 <- m1
    full_dat$m0 <- m0
    
    blb_reps <- replicate(r, {
      boot_idx <- sample(idx, size = n, replace = TRUE)
      boot_dat <- full_dat[boot_idx]
      phi1 <- (boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1
      phi0 <- (1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0
      mean(phi1) - mean(phi0)
    })
  })
  data.table(type = 'Full', time = times['elapsed'])
})
full <- rbindlist(full)

# CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'cblb_timing.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'cblb_timing.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    gamma <- grid_val$gamma
    b <- round(n^gamma)
    part_idx <- seq_len(b)
    subsets <- grid_val$subsets
    num_cores <- grid_val$num_cores

    cblb <- lapply(seq_len(replications), function(rp){
      times <- system.time({
        partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)

        blb_out <- parallel::mclapply(partitions, function(i){
          tmp_dat <- dat[i]
          m <- glm(y ~ Tr + X1 + X2, data = tmp_dat, family = 'gaussian')
          g <- glm(Tr ~ X1 + X2, data = tmp_dat, family = 'binomial')

          prop_score <- predict(g, type = 'response')
          m1 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 1))
          m0 <- predict(m, data.frame(X1 = tmp_dat$X1, X2 = tmp_dat$X2, Tr = 0))

          tmp_dat$prop_score <- prop_score
          tmp_dat$m1 <- m1
          tmp_dat$m0 <- m0

          blb_reps <- replicate(r, {
            boot_idx <- sample(part_idx, size = n, replace = TRUE)
            boot_dat <- tmp_dat[boot_idx]
            phi1 <- (boot_dat$Tr/boot_dat$prop_score)*(boot_dat$y - boot_dat$m1) + boot_dat$m1
            phi0 <- (1 - boot_dat$Tr)/(1 - boot_dat$prop_score)*(boot_dat$y - boot_dat$m0) + boot_dat$m0
            mean(phi1) - mean(phi0)
          })
        }, mc.cores = num_cores)

      })
      data.table(type = 'cBLB', time = times['elapsed'])
    })
    cblb <- rbindlist(cblb)
    cblb[, `:=`(subsets = subsets,
                num_cores = num_cores)]
    cblb
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'cblb_timing.rds'))
}

cblb[, `:=`(type = paste0('Subsets ', subsets, ', Cores ', num_cores))]
out <- rbindlist(list(cblb, full), fill = TRUE)

ggplot(out, aes(x = type, y = time)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
