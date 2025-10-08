#setwd("")
source("Functions_new.txt")
source("Example_1_cen.txt")


result <- EM_fit(transformed_data, init_param, lambda1 = 1, lambda2 = 10, maxit = 10^3, tolerance = 1e-6)


result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result


expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^4)
expected_time = sur_time_all(result, transformed_data,method = "MR", M = 100, n_samples = 10^4)
expected_time = sur_time_all(result, transformed_data,method = "IS_Weibull", M = 100, n_samples = 10^4)

d <- cbind(T,expected_time,delta)
#write.csv(d,"MR_epected_EX1_cen.csv")

coordinate_update <- function(j, Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l, lambda1, lambda2) {
  
  # Update sigma_l
  if(j==-1){
    
    initial_sigma <- sigma_l    
    lb <- max(0.001, sigma_l - 0.01)  
    ub <- sigma_l + 0.01             
    
    # Bundle additional parameters into a list
    args_sigma <- list(j = -1, Y = Y, delta = delta, eta_l = eta_l, X_l = X_l,
                       alpha_l = alpha_l, beta_l = beta_l, lambda1 = lambda1, lambda2 = lambda2)
    
    # Use optim with L-BFGS-B method to respect bounds
    result <- optim(
      par = initial_sigma,
      fn = function(sigma_l) {
        with(args_sigma, Q_l(j, Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l, lambda1, lambda2))
      },
      gr = function(sigma_l) {
        with(args_sigma, dQ_l_sigma(Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l))
      },
      method = "L-BFGS-B",
      lower = lb,
      upper = ub,
      control = list(pgtol = 1e-8)
    )
    
    return(result$par)
    
  }
  
  # Update alpha
  else if(j==0){
    initial_alpha <- alpha_l     # Starting value for alpha_l
    
    args_alpha <- list(j = 0, Y = Y, delta = delta, eta_l = eta_l,  X_l = X_l,
                       beta_l = beta_l, sigma_l = sigma_l, lambda1 = lambda1, lambda2 = lambda2)
    
    result_alpha <- optim(
      par = initial_alpha,
      fn = function(alpha_l) {
        with(args_alpha, Q_l(j = 0, Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l, lambda1, lambda2))
      },
      gr = function(alpha_l) {
        with(args_alpha, dQ_l(j = 0, Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l, lambda1, lambda2))
      },
      method = "BFGS",
      control = list(reltol = 1e-8)
    )
    return(result_alpha$par)
    
  }
  
  else{
    
    initial_beta <- beta_l
    
    #dQ_l_0 <- dQ_l(j, Y, delta, eta_l, X_l, alpha_l, construct_beta(0), sigma_l, lambda1 = 0, lambda2 = 0)
    
    # Set beta_lj to 0 if the derivative condition is met
    #if (abs(dQ_l_0) < lambda2) {
    #  return(0)   ## set beta_lj = 0
    #}
    
    # Arguments for optimization
    args_beta <- list(j = 1, Y = Y, delta = delta, eta_l = eta_l, X_l = X_l,
                      alpha_l = alpha_l, sigma_l = sigma_l,
                      lambda1 = lambda1, lambda2 = lambda2)
    
    # Perform optimization using `optim`
    result_beta <- optim(
      par = beta_l,
      fn = function(beta_l) {
        with(args_beta, Q_l(j, Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l, lambda1, lambda2))
      },
      gr = function(beta_l) {
        with(args_beta, dQ_l(j, Y, delta, eta_l, X_l, alpha_l, beta_l, sigma_l, lambda1, lambda2))
      },
      method = "BFGS",
      control = list(reltol = 1e-8)
    )
    
    # Return the optimized beta_lj
    return(result_beta$par)
    
  }
  
}
