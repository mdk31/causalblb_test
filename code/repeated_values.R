library(data.table)

# WTS: Does the regression in the subset, using the amount of the subset data, give the same results as the weighted regression. It shoudl because the information is the same but who knows lol

# Simulate data
set.seed(123)
n <- 10000
b <- round(n^0.7)
x <- rnorm(n)
y <- 2 + 3 * x + rnorm(n)

# Create a data frame
data <- data.table(x = x, y = y)
sub_data <- data[sample(1:n, size = b, replace = FALSE)]
rep_data <- sub_data[sample(1:b, size = n, replace = TRUE)]

lm_sub <- lm(y ~ x, data = sub_data)
summary(lm_sub)

# Aggregate data to get unique values and counts
aggregated_data <- aggregate(y ~ x, rep_data, FUN = function(x) c(count = length(x), mean = mean(x)))
aggregated_data <- data.frame(x = aggregated_data$x, y = aggregated_data$y[, "mean"], weights = aggregated_data$y[, "count"])

# Run weighted regression
lm_rep <- lm(y ~ x, data = aggregated_data, weights = weights)
summary(lm_rep)

summary(lm(y ~ x, data = data))
