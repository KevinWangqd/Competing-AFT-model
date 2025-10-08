#setwd("")
source("Functions.txt")
source("Example_1.txt")

result <- EM_fit(transformed_data, init_param, lambda1 = 0.5, lambda2 = 0.2, maxit = 10^3, tolerance = 1e-6)

result$param_new


summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result

##[[1]]$Group 1
#$Alpha         $Beta         $Sigma
# 1.626          1.17         0.917

#[[2]]$Group 2
#$Alpha         $Beta         $Sigma
# 1.219092      2.028185      0.993863

#[[3]]$Group 3

#$Alpha         $Beta         $Sigma
# 2.11          1.041          1.078


expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^3)

expected_time = sur_time_all(result, transformed_data,method = "MR", M = 100, n_samples = 10^4)

#d <- cbind(T,expected_time)

#write.csv(d,"Exp_expected_EX1_new.csv")


time_mat <- data.frame(T1,T2,T3)
time_mat$index <- ifelse(time_mat$T1<time_mat$T2&time_mat$T1<time_mat$T3, 1,
                         ifelse(time_mat$T2<time_mat$T1&time_mat$T2<time_mat$T3,2,3))

eta_m <- data.frame(result$eta_matrix)

eta_m$index <- ifelse(eta_m$X1 > eta_m$X2 & eta_m$X1 > eta_m$X3, 1,
                      ifelse(eta_m$X2 > eta_m$X1 & eta_m$X2 > eta_m$X3,2,3))

write.csv(eta_m,"eta_mat_Ex1.csv")
write.csv(time_mat,"time_Ex1.csv")






