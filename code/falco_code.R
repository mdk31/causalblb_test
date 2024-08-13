base_dir <- '~/Documents/HW/Research/CI/causalblb_test'
source(file.path(base_dir, 'code/helper_funcs.R'))

n <- 10000
te <- 0.8
sigma <- 1

dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)

causal_blb(dat, 
           y_method = 'glm', 
           prop_method = 'glm', 
           y_formula = y ~ Tr + X1 + X2, 
           prop_formula = Tr ~ X1 + X2, cores = 4)


