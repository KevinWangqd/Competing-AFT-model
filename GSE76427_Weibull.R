source("Functions_new.txt")
library(survival)

GSE76427 <- read.csv("surv_GSE76427.csv")

data <- GSE76427[,c(9,10,3,13:20)]
data[,3:10] <- scale(data[,3:10])
data$status <- data$OS_event
colnames(data)

## cox/weibull

cox_m <- coxph(Surv(OS_time, status) ~ Age + GenderM + CCDC107 + CXCL12  + GIGYF1 + GMNN + IFFO1, data = data)
summary(cox_m)
concordance(cox_m) #n= 468 Concordance= 0.621 se= 0.02699

w_m <- survreg(Surv(OS_time, status) ~ Age + GenderM + CCDC107 + CXCL12  + GIGYF1 + GMNN + IFFO1, data = data)
summary(w_m)
concordance(w_m) #Concordance= 0.621 se= 0.02701

#########################
Y <- log(data$OS_time)

delta <- data$status


X = cbind(data$Age, data$GenderM, data$CCDC107,
          data$CXCL12, data$GIGYF1, data$GMNN, data$IFFO1)


label = list(c(1), c(1,2), c(2), 
             c(1), c(2), c(2), c(1))  #0.747 expect

label = list(c(1), c(1,2), c(2), 
             c(1,2), c(2), c(2), c(1))  #748 exp

label = list(c(1), c(1,2), c(2), 
             c(1,2), c(2), c(2), c(1,2))  #755

label = list(c(1), c(1,2), c(2), 
             c(1,2), c(2), c(2), c(1,2)) 

transformed_data = data_transform(Y, delta, X, label)
#transformed_data
# with defaul initialization
init_param = param_init(X, label)
result <- EM_fit(transformed_data, init_param, lambda1 = 1, lambda2 = 0.2, maxit = 10^3, tolerance = 1e-4)
expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^3)
#expected_time = sur_time_all(result, transformed_data,method = "MR", M = 300, n_samples = 1000)
surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)

# Print the result
print(c_index$concordance)

#init_param$alpha <- list(1,2, 1.1)
#init_param$sigma <- list(1.2,1,1.1)
#init_param


#result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result
#expected_time = sur_time_all(result, transformed_data,method = "MR", M = 300, n_samples = 1000)



library(timeROC)
library(survival)
T = data$OS_time
median(T)
time_points <- quantile(T, probs = seq(0.18, 0.9, by = 0.01))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = T, delta = delta, marker = cox_m$linear.predictors, cause = 1, times = time_points)
mean(roc_cox$AUC)
# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = 1.16)   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = T, delta = delta, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = 1.16)  # Adjust time as needed
mean(roc_weibull$AUC)

roc_comp <- timeROC(T = T, delta = delta, marker = -expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = 1.16)  # Adjust time as needed

mean(roc_comp$AUC)

plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(T))   # Adjust time as needed
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed




jpeg("HCC.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_comp, title = "Time-dependent ROC (Cox vs Weibull)", time = 1.16, col = "orange", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

lines(roc_weibull$FP[,33],roc_weibull$TP[,33], col = "red", lwd = 2, lty = 2)

lines(roc_cox$FP[,33],roc_cox$TP[,33], col = "blue", lwd = 2)


# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.785)", "Cox (AUC = 0.739)", "Weibull (AUC = 0.729)"), col = c("orange","blue", "red"), lwd = 5,cex=1.5)

title("Time-dependent ROCs at median surival time (1.16 year)")

# Close the jpeg device
dev.off()






## Eta matrix plot

# Define the range of Y values
Y_values <- log(seq(0.1, 8, length.out = 200))  # 100 values from log(10) to log(730)

# Initialize an empty data frame to store results
results <- data.frame(Y = Y_values, mean_X1 = NA, mean_X2 = NA)

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
  #results$mean_X3[i] <- mean(output_e$X3, na.rm = TRUE)
}

# Print or visualize results

mean(results$mean_X1)
mean(results$mean_X2)
#mean(results$mean_X3)

colnames(results) <- c("Y", 
                       "Age, GenderM, CXCL12, IFFO1", 
                       "GenderM, CCDC107, CXCL12, GIGYF1, GMNN, IFFO1")



library(ggplot2)
library(reshape2)  # For reshaping data

# Reshape data to long format for ggplot
results_long <- melt(results, id.vars = "Y", variable.name = "Competing_Factors", value.name = "Winning_Probability")





jpeg("WP_time_HCC.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)
ggplot(results_long, aes(x = exp(Y), y = Winning_Probability, color = Competing_Factors, group = Competing_Factors)) +
  geom_line(size = 1) + 
  geom_point(size = 0.5, alpha = 0.7) +  # Optional: add points
  labs(x = "Time", y = "Winning Probability", title = "Winning probability for uncensored sample", color = "Competing Factors") +
  scale_color_manual(values = c("red", "blue", "forestgreen")) +  # Customize colors
  theme_bw() +  ylim(0,1)+# Use theme_bw to have more flexibility
  theme(
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.key.height = unit(1, "cm"),  # Increase vertical spacing between legend items
    legend.position = c(0.65, 0.5),  # Adjust legend position closer to the plot
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
latent_T <- as.data.frame(cbind(data$OS_time/result$eta_matrix,result$eta_matrix, data$OS_time, data$OS_event))

latent_T
latent_T <- latent_T[which(latent_T$V6==1),]
latent_T
colnames(latent_T) <- c("T1","T2","eta1","eta2","time","status")
latent_T

cutoff <- median(c(latent_T[,1],latent_T[,2]))
latent_T$G1 <- ifelse(latent_T$T1 < cutoff & latent_T$T2 > cutoff ,1,0)
latent_T$G2 <- ifelse(latent_T$T1 > cutoff & latent_T$T2 < cutoff ,1,0)

latent_T$G12 <- ifelse(latent_T$T1 < cutoff & latent_T$T2 < cutoff,1,0)


write.csv(latent_T,"latent_T_med_HCC.csv")

colSums(latent_T[,9:15])



