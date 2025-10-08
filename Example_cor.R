setwd("")
source("Functions_new_ln.txt")
library(MASS) 

set.seed(123)
n = 1500

sigma1 = 1
alpha1 = 1
beta11 = -3
beta12 = 2
beta14 = 1


sigma2 = 1
alpha2 = 1.5
beta21 = 2
beta22 = 2


sigma3 = 1.5
alpha3 = 1
beta31 = -2
beta32 = 3
beta33 = 2


X1 = rnorm(n)
X2 = rnorm(n)
X3 = rnorm(n)
X4 = rnorm(n)

X = cbind(X1, X2, X3, X4)

mu1 = X1*beta11 + X2*beta12 + X4*beta14 + alpha1
mu2 = X1*beta21 + X2*beta22 + alpha2
mu3 = X1*beta31 + X2*beta32 + X3*beta33 + alpha3


# Stack into matrix
mu_mat <- cbind(mu1, mu2, mu3)

# Define correlation and covariance matrix
rho <- 0  # Change as needed
Sigma <- matrix(c(
  sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma3,
  rho*sigma2*sigma1, sigma2^2, rho*sigma2*sigma3,
  rho*sigma3*sigma1, rho*sigma3*sigma2, sigma3^2
), nrow = 3, byrow = TRUE)



# Generate multivariate normal residuals
E <- mvrnorm(n = n, mu = c(0, 0, 0), Sigma = Sigma)

# Generate Y = mu + error
Y_mat <- mu_mat + E

# Extract
Y1 <- Y_mat[, 1]
Y2 <- Y_mat[, 2]
Y3 <- Y_mat[, 3]

Y <- pmin(Y1, Y2, Y3)





T <- exp(Y)


delta = rep(1,n)
## random censorship 
pi <- 0.1 # censorship rate

f <- function(lambda){
  mean(punif(T, min = min(T), max =lambda)-pi)
}

lambda = uniroot(f, interval = c(0.0001, 10^6))$root

set.seed(123)
C <- runif(n, min(T), lambda)
delta = rep(1,n)
for (i in 1:n) {
  if(T[i]>C[i]){
    T[i] = C[i]
    delta[i] = 0
  }
}

sum(delta)/n
sum(delta)
Y=log(T)


label = list(c(1,2,3), c(1,2,3), 3, 1)
### you may run the data_transform(T, delta, X, label) code and then fit this example
## For this example, the fisher information leads to NaN calculation
## here is the result for estimates.

# Print transformed data to check the result
transformed_data = data_transform(Y, delta, X, label)

# with defaul initialization
init_param = param_init(X, label)
init_param


## with initialization
p <- vector("list") 
p$alpha = list(0.1, 1.4, 1.2)
p$beta = list(c(-0.7,0.6,0), c(0.5, 0.5),c(-0.5, 0.7, 0.1))
p$sigma =list(1,1,1.2)
init_param = param_init(X, label, p)

init_param


result <- EM_fit(transformed_data, init_param, lambda1 = 2, lambda2 = 1, maxit = 2000, tolerance = 1e-6)


result$param_new



summary_EM_fit(result,transformed_data)


expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^3)

library(survival)
surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)
c_index



alpha1_hat <- result$param_new[[1]][[1]]
beta_hat <- result$param_new[[1]][[2]] # Only the first 3 correspond to beta11, beta12, beta14
beta_hat
# Compute differences
diff_alpha <- alpha1_hat - alpha1
diff_beta <- beta_hat - c(1,-3, 2, 1)

# Show results
cat("Alpha difference:", diff_alpha, "\n")
cat("Beta differences:", diff_beta, "\n")


sigma2_hat <- result$param_new[[2]][[1]]
beta_hat <- result$param_new[[2]][[2]] # Only the first 3 correspond to beta11, beta12, beta14
beta_hat
# Compute differences
diff_alpha <- sigma2_hat - sigma2
diff_beta <- beta_hat - c(1.5,2, 2)

# Show results
cat("Alpha difference:", diff_alpha, "\n")
cat("Beta differences:", diff_beta, "\n")


sigma3_hat <- result$param_new[[3]][[1]]
beta_hat <- result$param_new[[3]][[2]] # Only the first 3 correspond to beta11, beta12, beta14
beta_hat
# Compute differences
diff_alpha <- sigma3_hat - sigma3
diff_beta <- beta_hat - c(1, -2,3, 2)

# Show results
cat("Alpha difference:", diff_alpha, "\n")
cat("Beta differences:", diff_beta, "\n")
















Ex_sur_time <- function(X, alpha, beta, sigma, method, M = 100, samples) {    ## M is the upper limit for definite integral
  # Nested survival function
  # Nested survival function
  survival_function <- function(t) {
    L <- length(X)
    surv_val <- 1
    for (l in 1:L) {
      mu_l <- as.numeric(alpha[[l]] + X[[l]] %*% beta[[l]])
      z <- (log(t) - mu_l) / sigma[[l]]
      surv_val <- surv_val * (1 - pnorm(z))  # log-normal survival
    }
    return(surv_val)
  }
  # Transformed survival function for importance sampling
  transformed_survival_function <- function(y) {
    t <- exp(y)  # Transform back to the original scale
    L <- length(X)
    sum <- 0
    
    for (l in 1:L) {
      mu_l <- alpha[[l]] + X[[l]] %*% beta[[l]]
      sum <- sum + (t / as.vector(exp(mu_l)))^(1 / sigma[[l]])
    }
    
    S <- exp(-sum)  # Survival function
    return(S * t)  # Include the Jacobian term (e^y)
  }
  
  # Nested Mill's Ratio residual function
  MR_res <- function(M) {
    L <- length(X)
    sum1 <- 0
    sum2 <- 0
    for (l in 1:L) {
      mu_l <- as.numeric(alpha[[l]] + X[[l]] %*% beta[[l]])
      sum1 <- sum1 + (M / exp(mu_l))^(1 / sigma[[l]])
      sum2 <- sum2 + M^((1 / sigma[[l]]) - 1) / (sigma[[l]] * exp(mu_l / sigma[[l]]))
    }
    res <- exp(-sum1) / sum2  # Mill's Ratio residual
    return(res)
  }
  
  # Importance Sampling function
  IS_exp <- function(samples) {
    weights <- survival_function(samples) / dexp(samples, rate = 1)  # Importance weights
    estimate <- mean(weights)  # Monte Carlo estimate of the integral
    return(estimate)
  }
  # Importance Sampling function
  IS_Weibull <- function(samples) {
    weights <- survival_function(samples) / dweibull(samples, shape = 2, scale = 2)  # Importance weights
    estimate <- mean(weights)  # Monte Carlo estimate of the integral
    return(estimate)
  }
  
  
  IS_gumbel <- function(y_samples) {
    # Calculate weights
    weights <- transformed_survival_function(y_samples)/dgumbel(y_samples, loc = 0, scale = 1)
    estimate <- mean(weights)  # Monte Carlo estimate of the integral
    return(estimate)
  }
  
  if (method == "MR") {
    # For Mill's ratio and numerical approximation
    # Integrate the survival function up to M
    definite <- integrate(survival_function, lower = 0, upper = M, subdivisions = 10^8)$value
    
    # Calculate the Mill's ratio residual
    residual <- MR_res(M)
    
    # Expected survival time
    Ex_time <- definite + residual
  } 
  else if (method == "IS_exp") {
    # Use importance sampling for the integral
    Ex_time <- IS_exp(samples)
  } 
  else if (method == "IS_Weibull") {
    # Use importance sampling for the integral
    Ex_time <- IS_Weibull(samples)
  } 
  else if (method == "IS_gumbel") {  # Importance Sampling with Gumbel distribution
    Ex_time <- IS_gumbel(samples)
  } 
  else {
    stop("Invalid method. Use 'MR' for Mills Ratio or 'IS' for Importance Sampling.")
  }
  
  
  
  return(Ex_time)
}




