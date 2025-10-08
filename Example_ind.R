setwd("")
source("Functions_new_ln.txt")

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

Y1 = rep(NA,n)
Y2 = rep(NA,n)
Y3 = rep(NA,n)
Y = rep(NA,n)


for (i in 1:n) {
  Y1[i] = rnorm(1,mean=mu1[i], sd = sigma1)
  Y2[i] = rnorm(1,mean=mu2[i], sd = sigma2)
  Y3[i] = rnorm(1,mean=mu3[i], sd = sigma3)
  Y[i] = min(Y1[i], Y2[i], Y3[i])
}

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


result <- EM_fit(transformed_data, init_param, lambda1 = 0.5, lambda2 = 0.2, maxit = 1000, tolerance = 1e-6)


result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result


expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^4)


surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)

# Print the result
print(c_index$concordance)


expected_time = sur_time_all(result, transformed_data,method = "MR", M = 100, n_samples = 10^4)
expected_time = sur_time_all(result, transformed_data,method = "IS_Weibull", M = 100, n_samples = 10^4)

d <- cbind(T,expected_time,delta)
#write.csv(d,"Exp_epected_EX2_cen_1.csv")

#write.csv(d,"MR_epected_EX2_cen_1.csv")



result$param_new[[1]]




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
diff_alpha <- alpha1_hat - sigma2
diff_beta <- beta_hat - c(1.5,2, 2)

# Show results
cat("Alpha difference:", diff_alpha, "\n")
cat("Beta differences:", diff_beta, "\n")


sigma3_hat <- result$param_new[[3]][[1]]
beta_hat <- result$param_new[[3]][[2]] # Only the first 3 correspond to beta11, beta12, beta14
beta_hat
# Compute differences
diff_alpha <- alpha1_hat - sigma3
diff_beta <- beta_hat - c(1, -2,3, 2)

# Show results
cat("Alpha difference:", diff_alpha, "\n")
cat("Beta differences:", diff_beta, "\n")
