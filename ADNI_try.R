library(readxl)
library(survival)
data <- read_excel("final_data3.xlsx")

source("Functions_new.txt")
#source("Functions_new_ln.txt")
data$RAVLT_perc_forgetting <- ifelse(data$RAVLT_perc_forgetting<0,0,data$RAVLT_perc_forgetting)
data <- data[data$survival_time>0,]

summary(data$APOE4)
mean(data$AGE)
range(data$AGE)
summary(as.factor(data$PTEDUCAT))


#增加的变量
data$ADAS13 <- scale(data$ADAS13)
data$ICV <- scale(data$ICV)
data$Hippocampus <- scale(data$Hippocampus)
data$MidTemp <- scale(data$MidTemp)
data$AGE <- scale(data$AGE)

#初次的变量
data$RAVLT_immediate <- scale(data$RAVLT_immediate)
data$RAVLT_learning <- scale(data$RAVLT_learning)
data$FAQ <- scale(data$FAQ)
#data$age <- scale(data$AGE)
#data$PTGENDER <- scale(data$PTGENDER)
data$PTEDUCAT <- scale(data$PTEDUCAT)

#变量和label

T <- data$survival_time
Y <- log(data$survival_time)
delta <- data$censor

#X = cbind(data$RAVLT_immediate,
#          data$PTGENDER, data$RAVLT_perc_forgetting,
#          data$LDELTOTAL, data$TRABSCOR,
#          data$Hippocampus, data$Entorhinal,data$Fusiform,
#          data$ICV)
#label = list(c(1,2), 
#             c(2),c(3),
#             c(1),c(3),
#             c(4),c(1,2,3),c(1,2),
#             c(2,3))



X = cbind(data$AP1, data$AP2, data$ADAS13, data$RAVLT_immediate, 
          data$RAVLT_learning, data$FAQ, data$PTGENDER, data$PTEDUCAT, 
          data$AGE, data$Hippocampus, data$ICV, data$MidTemp )
label = list(c(1), c(1), c(2,3), c(1,3),
             c(1), c(1,2), c(2,3), c(2,3),
             c(1,2), c(1,2), c(2,3), c(1,3))

#init_param$beta <- list(c(0,0,0,0),c(0,0,0),c(0,0,0,0,0,0),c(0),c(0,0),c(0,0))



# Print transformed data to check the result
transformed_data = data_transform(Y, delta, X, label)
#transformed_data
# with defaul initialization
init_param = param_init(X, label)
#init_param
init_param$alpha <- list(8,9,8,8)
init_param$sigma <- list(0.8,0.9,1.1,1)
# Compute the error matrix and eta_matrix (E-step)


result <- EM_fit(transformed_data, init_param, lambda1 = 1, lambda2 = 0.2, maxit = 10^3, tolerance = 1e-6)

result$param_new



summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result
expected_time = sur_time_all(result, transformed_data,method = "IS_Weibull", M = 100, n_samples = 1000)
#expected_time = sur_time_all(result, transformed_data,method = "MR", M = 300, n_samples = 896)
#expected_time


surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)

# Print the result
print(c_index$concordance)





w_m <- survreg(Surv(T, delta)~data$AP1+data$AP2+ data$ADAS13+ data$RAVLT_immediate+ 
                 data$RAVLT_learning+data$FAQ+ data$PTGENDER+ data$PTEDUCAT + data$AGE+
                 data$Hippocampus + data$ICV +data$MidTemp , data, dist="weibull")
summary(w_m)

concordance(w_m) 

g_m <- survreg(Surv(T, delta)~data$AP1+data$AP2+ data$ADAS13+ data$RAVLT_immediate+ 
                 data$RAVLT_learning+data$FAQ+ data$PTGENDER+ data$PTEDUCAT + data$AGE+
                 data$Hippocampus + data$ICV +data$MidTemp , data, dist="lognormal")
summary(g_m)

concordance(g_m) 

cox_m <- coxph(Surv(T, delta)~data$AP1+data$AP2+ data$ADAS13+ data$RAVLT_immediate+ 
                 data$RAVLT_learning+data$FAQ+ data$PTGENDER+ data$PTEDUCAT + data$AGE+
                 data$Hippocampus + data$ICV +data$MidTemp, data)

summary(cox_m)


library(timeROC)
library(survival)
time_points <- quantile(T, probs = seq(0.1, 0.99, by = 0.01))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = T, delta = delta, marker = cox_m$linear.predictors, cause = 1, times = time_points)
mean(roc_cox$AUC)
# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = 365)   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = T, delta = delta, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = 365)  # Adjust time as needed
mean(roc_weibull$AUC)

roc_comp <- timeROC(T = T, delta = delta, marker = -expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = 365)  # Adjust time as needed

mean(roc_comp$AUC)

plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(T))   # Adjust time as needed
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed






jpeg("ADNIoneyear.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_comp, title = "Time-dependent ROC (Cox vs Weibull)", time = 365, col = "orange", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

lines(roc_weibull$FP[,4],roc_weibull$TP[,4], col = "red", lwd = 2, lty = 2)

lines(roc_cox$FP[,4],roc_cox$TP[,4], col = "blue", lwd = 2)


# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.864)", "Cox (AUC = 0.856)", "Weibull (AUC = 0.851)"), col = c("orange","blue", "red"), lwd = 5,cex=1.5)

title("Time-dependent ROCs at 1-year surival time (365)")

# Close the jpeg device
dev.off()






## Eta matrix plot

# Define the range of Y values
Y_values <- log(seq(10, 365*3, length.out = 200))  # 100 values from log(10) to log(730)

# Initialize an empty data frame to store results
results <- data.frame(Y = Y_values, mean_X1 = NA, mean_X2 = NA, mean_X3 = NA)

# Loop over each Y value
for (i in seq_along(Y_values)) {
  # Compute eta_baseline for the current Y value
  eta_baseline <- eta_matrix_fun(Y = rep(Y_values[i], length(transformed_data$Y)), 
                                 delta = transformed_data$delta, 
                                 X = transformed_data$X, 
                                 alpha = lapply(result$param_new, function(x) x[[2]][1]), 
                                 beta = lapply(result$param_new, function(x) x[[2]][-1]), 
                                 sigma = lapply(result$param_new, function(x) x[[1]]))
  
  # Convert to data frame and filter for delta == 1
  output <- data.frame(eta_baseline)
  output$ind <- transformed_data$delta
  output_e <- subset(output, ind == 1)
  
  # Store mean values
  results$mean_X1[i] <- mean(output_e$X1, na.rm = TRUE)
  results$mean_X2[i] <- mean(output_e$X2, na.rm = TRUE)
  results$mean_X3[i] <- mean(output_e$X3, na.rm = TRUE)
}

# Print or visualize results

mean(results$mean_X1)
mean(results$mean_X2)
mean(results$mean_X3)

colnames(results) <- c("Y", 
                       "APOEe4, RAVLT (immediate), RAVLT (learning), \nFAQ, Age, Hippocampus, MidTemp", 
                       "ADAS13, FAQ, GenderM, Education, \nAge, Hippocampus, ICV", 
                       "ADAS13, RAVLT (immediate), GenderM, \nEducation, ICV, MidTemp")



library(ggplot2)
library(reshape2)  # For reshaping data

# Reshape data to long format for ggplot
results_long <- melt(results, id.vars = "Y", variable.name = "Competing_Factors", value.name = "Winning_Probability")





jpeg("WP_time_ADNI_1.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)
ggplot(results_long, aes(x = exp(Y), y = Winning_Probability, color = Competing_Factors, group = Competing_Factors)) +
  geom_line(size = 1) + 
  geom_point(size = 0.5, alpha = 0.7) +  # Optional: add points
  labs(x = "Time", y = "Winning Probability", title = "Winning probability for uncensored sample", color = "Competing Factors") +
  scale_color_manual(values = c("red", "blue", "forestgreen")) +  # Customize colors
  theme_bw() +  ylim(0,0.51)+# Use theme_bw to have more flexibility
  theme(
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.key.height = unit(1, "cm"),  # Increase vertical spacing between legend items
    legend.position = c(0.7, 0.16),  # Adjust legend position closer to the plot
    legend.margin = margin(0, 0, 0, 1),  # Reduce space around the legend
    legend.spacing.x = unit(0.5, "cm"),  # Space between legend items
    plot.title = element_text(hjust = 0.2, size = 14, face = "bold", margin = margin(b = 15)),  # Add space between title and plot
    panel.border = element_rect(color = "black", size = 1),  # Add border around the plot (thicker)
    plot.margin = margin(20, 20, 20, 20),  # Add more margin space to ensure plot fits with border
    axis.title.x = element_text(margin = margin(t = 10)),  # Space between x axis label and plot
    axis.title.y = element_text(margin = margin(r = 10))  # Space between y axis label and plot
  )
dev.off()


## Venn graph

latent_T <- as.data.frame(cbind(data$survival_time/result$eta_matrix,result$eta_matrix, data$survival_time, data$censor))

latent_T <- latent_T[which(latent_T$V8==1),]
colnames(latent_T) <- c("T1","T2","T3","eta1","eta2","eta3","time","status")
cutoff <- median(c(latent_T[,1],latent_T[,2],latent_T[,3]))
latent_T$G1 <- ifelse(latent_T$T1 < cutoff & latent_T$T2 > cutoff & latent_T$T3 > cutoff,1,0)
latent_T$G2 <- ifelse(latent_T$T1 > cutoff & latent_T$T2 < cutoff & latent_T$T3 > cutoff,1,0)
latent_T$G3 <- ifelse(latent_T$T1 > cutoff & latent_T$T2 > cutoff & latent_T$T3 < cutoff,1,0)

latent_T$G12 <- ifelse(latent_T$T1 < cutoff & latent_T$T2 < cutoff & latent_T$T3 > cutoff,1,0)
latent_T$G23 <- ifelse(latent_T$T1 > cutoff & latent_T$T2 < cutoff & latent_T$T3 < cutoff,1,0)
latent_T$G13 <- ifelse(latent_T$T1 < cutoff & latent_T$T2 > cutoff & latent_T$T3 < cutoff,1,0)

latent_T$G123 <- ifelse(latent_T$T1 < cutoff & latent_T$T2 < cutoff & latent_T$T3 < cutoff,1,0)

#write.csv(latent_T,"latent_T_med.csv")

colSums(latent_T[,9:15])






