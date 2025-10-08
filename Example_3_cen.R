setwd("")
source("Functions.txt")
source("Example_3_cen1.txt")
source("Example_3_cen3.txt")


#result <- EM_fit(transformed_data, init_param, lambda1 = 0.5, lambda2 = 0.2, maxit = 1000, tolerance = 1e-6)

result <- EM_fit(transformed_data, init_param, lambda1 = 2, lambda2 = 1, maxit = 1000, tolerance = 1e-6)
result <- EM_fit(transformed_data, init_param, lambda1 = 8, lambda2 = 4, maxit = 1000, tolerance = 1e-6)


result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result


expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^4)
expected_time = sur_time_all(result, transformed_data,method = "MR", M = 100, n_samples = 10^4)
expected_time = sur_time_all(result, transformed_data,method = "IS_Weibull", M = 100, n_samples = 10^4)

d <- cbind(T,expected_time,delta)
write.csv(d,"Exp_expected_EX3_cen_1.csv")

#write.csv(d,"MR_epected_EX2_cen_1.csv")















