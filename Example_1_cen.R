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

