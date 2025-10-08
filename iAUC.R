library(timeROC)
library(survival)

T = MR_expected_EX1_cen$T
delta = MR_expected_EX1_cen$delta
data <- X
#T = MR_expected_EX1$T
#delta = rep(1,1000)
data <- cbind(T,X, delta)
cox_m <- coxph(Surv(T,delta)~X1+X2+X3, data)
summary(cox_m)

w_m <- survreg(Surv(T,delta)~X1+X2+X3, data = data, dist="weibull")
summary(w_m)

concordance(w_m) 


time_points <- quantile(T, probs = seq(0.1, 0.99, by = 0.01))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = T, delta = delta, marker = cox_risk_scores, cause = 1, times = time_points)
#roc_cox <- timeROC(T = T, delta = delta, marker = cox_m$linear.predictors, cause = 1, times = time_points)
mean(roc_cox$AUC)
# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(T))   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = T, delta = delta, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed
mean(roc_weibull$AUC)

roc_comp <- timeROC(T = T, delta = delta, marker = -MR_expected_EX1_cen$expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed

mean(roc_comp$AUC)



jpeg("Ex1.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_cox, title = "Time-dependent ROC (Cox vs Weibull)", time = median(T), col = "blue", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

# 3. Add the ROC curve for the Weibull model at the median survival time
lines(roc_weibull$FP[,45],roc_weibull$TP[,45], col = "red", lwd = 2, lty = 2)

lines(roc_comp$FP[,45],roc_comp$TP[,45], col = "orange", lwd = 2)

# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.82)", "Cox (AUC = 0.78)", "Weibull (AUC = 0.78)"), col = c("orange","blue", "red"), lwd = 2, cex = 1.5)

title("Time-dependent ROCs at median surival time (0.636)")

# Close the jpeg device
dev.off()


roc_cox$AUC
roc_weibull$AUC
roc_comp$AUC
td_AUC <- data.frame(cbind(time_points,roc_comp$AUC,roc_weibull$AUC,roc_cox$AUC))

# Reshape the data to long format
data_long <- td_AUC %>%
  pivot_longer(cols = V2:V4, 
               names_to = "group", 
               values_to = "value")

# Plot the data
ggplot(data_long, aes(x = time_points, y = value, color = group, group = group)) +
  geom_line() +
  theme_minimal() +
  labs(title = "", 
       x = "Time", 
       y = "Value") +
  theme(legend.title = element_blank())


##Ex2

T = MR_epected_EX2_cen_1$T
delta = MR_epected_EX2_cen_1$delta
data <- cbind(T,X_Ex2, delta)
cox_m <- coxph(Surv(T,delta)~X1+X2+X3+X4, data)
summary(cox_m)

median(T)

w_m <- survreg(Surv(T,delta)~X1+X2+X3+X4, data = data, dist="weibull")
summary(w_m)

concordance(w_m) 


time_points <- quantile(T, probs = seq(0.1, 0.99, by = 0.01))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = T, delta = delta, marker = cox_risk_scores, cause = 1, times = time_points)

# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(T))   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = T, delta = delta, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed

roc_comp <- timeROC(T = T, delta = delta, marker = -MR_epected_EX2_cen_1$expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed

mean(roc_cox$AUC)
mean(roc_weibull$AUC)
mean(roc_comp$AUC)

jpeg("Ex2.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_cox, title = "Time-dependent ROC (Cox vs Weibull)", time = median(T), col = "blue", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

# 3. Add the ROC curve for the Weibull model at the median survival time
lines(roc_weibull$FP[,45],roc_weibull$TP[,45], col = "red", lwd = 2, lty = 2)

lines(roc_comp$FP[,45],roc_comp$TP[,45], col = "orange", lwd = 2)

# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.95)", "Cox (AUC = 0.88)", "Weibull (AUC = 0.88)"), col = c("orange","blue", "red"), lwd = 2, cex=1.5)

title("Time-dependent ROCs at median surival time (0.157)")

# Close the jpeg device
dev.off()


## Ex3

T = Exp_expected_EX3_cen_1$T
delta = Exp_expected_EX3_cen_1$delta
data <- cbind(T,X_Ex3, delta)
cox_m <- coxph(Surv(T,delta)~X1+X2+X3+X4+X5+X6, data)
summary(cox_m)

median(T)

w_m <- survreg(Surv(T,delta)~X1+X2+X3+X4+X5+X6, data = data, dist="weibull")
summary(w_m)

concordance(w_m) 


time_points <- quantile(T, probs = seq(0.1, 0.99, by = 0.01))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = T, delta = delta, marker = cox_risk_scores, cause = 1, times = time_points)

# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(T))   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = T, delta = delta, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed

roc_comp <- timeROC(T = T, delta = delta, marker = -Exp_expected_EX3_cen_1$expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed

mean(roc_cox$AUC)
mean(roc_weibull$AUC)
mean(roc_comp$AUC)


jpeg("Ex3.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_cox, title = "Time-dependent ROC (Cox vs Weibull)", time = median(T), col = "blue", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

# 3. Add the ROC curve for the Weibull model at the median survival time
lines(roc_weibull$FP[,45],roc_weibull$TP[,45], col = "red", lwd = 2, lty = 2)

lines(roc_comp$FP[,45],roc_comp$TP[,45], col = "orange", lwd = 2)

# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.92)", "Cox (AUC = 0.83)", "Weibull (AUC = 0.83)"), col = c("orange","blue", "red"), lwd = 5, cex=1.5)

title("Time-dependent ROCs at median surival time (0.098)")

# Close the jpeg device
dev.off()


td_AUC <- data.frame(cbind(time_points,roc_comp$AUC,roc_weibull$AUC,roc_cox$AUC))

# Reshape the data to long format
data_long <- td_AUC %>%
  pivot_longer(cols = V2:V4, 
               names_to = "group", 
               values_to = "value")

# Plot the data
ggplot(data_long, aes(x = time_points, y = value, color = group, group = group)) +
  geom_line() +
  theme_minimal() +
  labs(title = "", 
       x = "Time", 
       y = "Value") +
  theme(legend.title = element_blank())
