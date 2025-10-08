source("Functions_new.txt")
#source("Functions_new_ln.txt")
library(survival)
selected_vars <- c("GenderM", "Age", "status", "time", "ID_REF",
                   "11725985_s_at", "11726583_s_at", "11757883_s_at", '11761102_x_at', 
                   "11746889_a_at", "11722403_a_at", "11743948_at", "11755743_a_at","11729772_at",
                   "11745244_x_at", "11728500_a_at", "11757634_a_at", "11723193_s_at", "11729349_at", "11718130_s_at")

#IGHG1 [1] "11745244_x_at" "11750231_x_at" "11754032_x_at" "11756093_x_at" "11759852_x_at" "11760819_x_at" "11760929_x_at" "11763104_at"   "11763113_at"  
#"11763114_x_at" "11763115_a_at" "11763116_x_at" "11763117_x_at" "11763676_at"  

#IL1R2 "11728500_a_at"

#LCN2 "11757634_a_at"

#LTF "11723193_s_at" "11755996_s_at"

#MMP8 "11729349_at"   "11729350_a_at" "11729351_a_at" "11748210_a_at"

#OLFM4 "11718130_s_at" "11718131_a_at"

GSE65682 <- GSE[, selected_vars]
str(GSE65682)
head(GSE65682)
sum(GSE$status)
## Standardization
data <- GSE65682 # Consider sepsis
data$time <- as.numeric(as.character(data$time))
data$time[ data$time <= 0] <- NA
data <- data[ !is.na(data$time) & !is.na(data$status), ]

# 排除不标准化的列
exclude <- c("status", "time", "ID_REF")
# 数值型变量
data$Age <- scale(data$Age)
data$NONO <- scale(data$'11725985_s_at')
data$CKAP4 <- scale(data$'11726583_s_at')
data$RNF4 <- scale(data$'11757883_s_at')
data$FCAR <- scale(data$'11761102_x_at')
data$PLEKHO1 <- scale(data$"11722403_a_at")
data$BMP6 <- scale(data$"11743948_at")
data$RNASE2 <- scale(data$"11755743_a_at" )

data$IGHG1 <- scale(data$"11745244_x_at" )
data$IL1R2 <- scale(data$"11728500_a_at" )
data$LCN2 <- scale(data$"11757634_a_at" )
data$LTF <- scale(data$"11723193_s_at" )
data$MMP8 <- scale(data$"11729349_at" )
data$OLFM4 <- scale(data$"11718130_s_at" )

data <- read.csv("final_data.csv")
## cox/weibull
cox_m <- coxph(Surv(time, status) ~ Age + GenderM + NONO + CKAP4 + RNF4 + FCAR + PLEKHO1 + BMP6 + RNASE2+
                 IGHG1 + IL1R2 + LCN2 + LTF + MMP8 + OLFM4, data = data)
summary(cox_m)
concordance(cox_m) #n= 468 Concordance= 0.653 se= 0.02699

w_m <- survreg(Surv(time, status) ~ Age + GenderM + NONO + CKAP4 + RNF4 + FCAR+ PLEKHO1 + BMP6 + RNASE2+
                 IGHG1 + IL1R2 + LCN2 + LTF + MMP8 + OLFM4, data = data)
summary(w_m)
concordance(w_m) #Concordance= 0.654 se= 0.02701
exp(0.339)
#########################
Y <- log(data$time)

delta <- data$status



X = cbind(data$Age, data$GenderM, 
          data$NONO, data$CKAP4, data$RNF4, data$FCAR,
          data$PLEKHO1, data$BMP6, data$RNASE2,
          data$IGHG1, data$IL1R2, data$LCN2, 
          data$LTF, data$MMP8, data$OLFM4)



label = list(c(1), c(2),
             c(1), c(3), c(1,3), c(1,3), 
             c(1,2), c(1,2,3), c(1),
             c(2,3), c(3), c(2),
             c(1,2,3), c(1,3), c(1)) #722

transformed_data = data_transform(Y, delta, X, label)

init_param = param_init(X, label)

result <- EM_fit(transformed_data, init_param, lambda1 = 1, lambda2 = 0.2, maxit = 10^3, tolerance = 1e-4)


#expected_time = sur_time_all(result, transformed_data,method = "MR", M = 300, n_samples = 1000)
expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^3)

surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)

# Print the result
print(c_index$concordance)

result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result


library(timeROC)
library(survival)


time_points <- quantile(data[which(data$status==1),]$time, probs = seq(0.1, 0.9, by = 0.1))
time_points[1] <- 1.1
# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = data$time, delta = data$status, marker = cox_risk_scores, cause = 1, times = time_points)
#roc_cox <- timeROC(T = T, delta = delta, marker = cox_m$linear.predictors, cause = 1, times = time_points)
mean(roc_cox$AUC)
# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(time_points))   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = data$time, delta = data$status, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(time_points))  # Adjust time as needed
mean(roc_weibull$AUC)

roc_comp <- timeROC(T = data$time, delta = data$status, marker = -expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(time_points))  # Adjust time as needed

mean(roc_comp$AUC)
roc_weibull$AUC
roc_comp$AUC

T = data$time
median(time_points)


jpeg("sepsis.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_comp, title = "Time-dependent ROC (Cox vs Weibull)", time = 8, col = "orange", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

lines(roc_weibull$FP[,5],roc_weibull$TP[,5], col = "red", lwd = 2, lty = 2)

lines(roc_cox$FP[,5],roc_cox$TP[,5], col = "blue", lwd = 2)


# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.743)", "Cox (AUC = 0.689)", "Weibull (AUC = 0.689)"), 
       col = c("orange","blue", "red"), lwd = 5, cex=1.5)

title("Time-dependent ROCs at median surival time (8 days)")

# Close the jpeg device
dev.off()



## Eta matrix plot

# Define the range of Y values
Y_values <- log(seq(1, 28, length.out = 200))  # 100 values from log(10) to log(730)

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



colnames(results) <- c("Y", "Age, NONO, RNF4, FCAR, PLEKHO1, \nBMP6, RNASE2, LTF, MMP8, OLFM4",
                       "GenderM, PLEKHO1, BMP6, \nIGHG1, LCN2, LTF", 
                       "CKAP4, RNF4, FCAR, BMP6, \nIGHG1, IL1R2, LTF, MMP8"
                       
)




library(ggplot2)
library(reshape2)  # For reshaping data

# Reshape data to long format for ggplot
results_long <- melt(results, id.vars = "Y", variable.name = "Competing_Factors", value.name = "Winning_Probability")



jpeg("WP_time_sepsis.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)
ggplot(results_long, aes(x = exp(Y), y = Winning_Probability, color = Competing_Factors, group = Competing_Factors)) +
  geom_line(size = 1) + 
  geom_point(size = 0.5, alpha = 0.7) +  # Optional: add points
  labs(x = "Time", y = "Winning Probability", title = "Winning probability for uncensored sample", color = "Competing Factors") +
  scale_color_manual(values = c("red", "blue", "forestgreen", "purple")) +  # Customize colors
  theme_bw() + ylim(0,0.65) + # Use theme_bw to have more flexibility
  theme(
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.key.height = unit(0.8, "cm"),  # Increase vertical spacing between legend items
    legend.position = c(0.7, 0.33),  # Adjust legend position closer to the plot
    legend.margin = margin(0, 0, 0, 10),  # Reduce space around the legend
    legend.spacing.x = unit(0.1, "cm"),  # Space between legend items
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 15)),  # Add space between title and plot
    panel.border = element_rect(color = "black", size = 1),  # Add border around the plot (thicker)
    plot.margin = margin(20, 20, 20, 20),  # Add more margin space to ensure plot fits with border
    axis.title.x = element_text(margin = margin(t = 10)),  # Space between x axis label and plot
    axis.title.y = element_text(margin = margin(r = 10))  # Space between y axis label and plot
  )
dev.off()




## Venn graph

latent_T <- as.data.frame(cbind(data$time/result$eta_matrix,result$eta_matrix, data$time, data$status,data$ID_REF))

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

write.csv(latent_T,"latent_T_med_spesis.csv")

colSums(latent_T[,9:15])
