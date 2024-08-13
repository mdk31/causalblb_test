library(pbapply)
library(data.table)

N <- 50000
r <- 100
gamma <- 0.7
b <- round(N^gamma)
subsets <- 25
replicates <- 1000



out <- pblapply(1:replicates, function(rep){
  x <- rnorm(N)
  
  out <- lapply(1:subsets, function(sub){
    x_sub <- sample(x, size = b, replace = FALSE)
    
    blb <- replicate(r, {
      x_sub_boot <- sample(x_sub, size = N, replace = TRUE)
      mean(x_sub_boot)
    })
    ci <- boot:::perc.ci(blb)
    data.table(lower_ci = ci[4],
               upper_ci = ci[5])
  })
  out <- rbindlist(out)
  out <- out[, .(lower_ci = mean(lower_ci), upper_ci = mean(upper_ci))]
  out[, `:=`(rep = rep)]
}, cl = 4)

out <- rbindlist(out)

out[, .(mean(lower_ci <= 0 & upper_ci >= 0))]
