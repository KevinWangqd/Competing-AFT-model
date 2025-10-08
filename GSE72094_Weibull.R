source("Functions_new.txt")
library(survival)


GSE72094 <- read.csv("final_data.csv")
data <- GSE72094[,c(3,4,11,12,13)]
data$CCNA2 <- GSE72094$merck.NM_001237_a_at
data$AURKA <- GSE72094$merck2.BE856617_at
data$AURKB <- GSE72094$merck.NM_004217_at
data$FEN1 <- GSE72094$merck2.XM_937756_a_at
data$CCND3 <- GSE72094$merck2.BQ669293_at
data$NCALD <- GSE72094$merck2.AF251061_at
data$MACF1 <- GSE72094$merck2.AB029290_at
data$LRC4 <- GSE72094$merck2.NM_007360_at
data$NLRC4 <- GSE72094$merck.NM_021209_s_at
data$PLEKHN1 <- GSE72094$merck.NM_032129_at
data$RASIP1 <- GSE72094$merck.NM_017805_at
data$SPP1 <- GSE72094$merck2.CA447290_at
data$GPT2 <- GSE72094$merck2.BX099266_at
data$SGPL1 <- GSE72094$merck.ENST00000299297_at
data$PCOLCE2 <- GSE72094$merck2.AK223633_at

data[,c(2,6:20)] <- scale(data[,c(2,6:20)])

#CCNA2: merck.NM_001237_a_at
#AURKA [1] "merck2.BE856617_at"   "merck.BC050630_at"    "merck.NM_017900_s_at" "merck.NM_198436_s_at" "merck.XM_933637_s_at"
#AURKB "merck.NM_004217_at"
#FEN1 [1] "merck2.XM_937756_a_at" "merck.NM_004111_s_at" 
#CD44 "merck2.BM550721_at"    "merck2.BM792065_at"    "merck2.CR621045_at"    "merck.BC004372_a_at"   "merck.NM_001001389_at"
#CCND3 [1] "merck2.BQ669293_at" "merck.NM_001760_at"
# NCALD [1] "merck2.AF251061_at"    "merck.BC063428_s_at"   "merck.NM_001040630_at"
#MACF1 "merck2.AB029290_at"   "merck2.BP284289_a_at" "merck2.BQ651417_a_at" "merck2.CX870374_a_at" "merck2.NM_033024_at"  "merck.AF141968_a_at"  "merck.AK023406_a_at"  "merck.AK023821_at"   
#LRC4 [1] "merck2.NM_007360_at"   "merck2.NM_013431_at"   "merck2.NM_013431_x_at" "merck.NM_007360_at"    "merck.NM_013431_at"    "merck.NM_021209_s_at"  "merck.NM_176677_at"   
#NLRC4 [1] "merck.NM_021209_s_at"
#PLEKHN1 [1] "merck.NM_032129_at"
#RASIP1 [1] "merck.NM_017805_at"
#SPP1 [1] "merck2.CA447290_at"   "merck2.DQ892544_at"   "merck.BG211014_a_at"  "merck.NM_000582_at"   "merck.NM_024790_s_at"
#GPT2 [1] "merck2.BX099266_at"   "merck.BC062555_s_at"  "merck.NM_001147_at"   "merck.NM_133443_s_at"
#SGPL1 [1] "merck.ENST00000299297_at" "merck.NM_003901_a_at"    
#PCOLCE2 [1] "merck2.AK223633_at" "merck.NM_013363_at"


## cox/weibull
cox_m <- coxph(Surv(time_day, status) ~ age + GenderM + CCNA2 + AURKA + AURKB + FEN1 + CCND3 + NCALD + MACF1 + LRC4 + NLRC4 + PLEKHN1+
                 RASIP1 + SPP1 + GPT2 + SGPL1 + PCOLCE2, data = data)
summary(cox_m)
concordance(cox_m) #n= 468 Concordance= 0.621 se= 0.02699

w_m <- survreg(Surv(time_day, status) ~ age + GenderM + CCNA2 + AURKA + AURKB + FEN1 + CCND3 + NCALD + MACF1 + LRC4 + NLRC4 + PLEKHN1+
                 RASIP1 + SPP1 + GPT2 + SGPL1 + PCOLCE2, data = data)
summary(w_m)
concordance(w_m) #Concordance= 0.621 se= 0.02701

#########################
Y <- log(data$time_day)

delta <- data$status


X = cbind(data$age, data$GenderM, data$CCNA2, data$AURKA, data$AURKB, 
          data$FEN1, data$CCND3, data$NCALD, data$MACF1,  data$LRC4, 
          data$NLRC4, data$PLEKHN1, data$RASIP1, data$SPP1, data$GPT2,
          data$SGPL1, data$PCOLCE2)


label = list(c(1), c(1,2,3), c(1,3), c(1,3), c(2,3),
             c(1), c(2,3), c(2,3), c(2), c(1,2),
             c(1,2), c(2), c(1,2,3), c(3), c(2,3),
             c(3), c(2)) #741


transformed_data = data_transform(Y, delta, X, label)
#transformed_data
# with defaul initialization
init_param = param_init(X, label)
result <- EM_fit(transformed_data, init_param, lambda1 = 1, lambda2 = 0.1, maxit = 10^3, tolerance = 1e-4)
#expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^3)
expected_time = sur_time_all(result, transformed_data,method = "MR", M = 300, n_samples = 1000)
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



library(timeROC)
library(survival)
T <- data$time_day
median(T)
time_points <- quantile(T, probs = seq(0.1, 0.99, by = 0.01))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = T, delta = delta, marker = cox_m$linear.predictors, cause = 1, times = time_points)
mean(roc_cox$AUC)
# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = 824)   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = T, delta = delta, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = 824)  # Adjust time as needed
mean(roc_weibull$AUC)

roc_comp <- timeROC(T = T, delta = delta, marker = -expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = 824)  # Adjust time as needed

mean(roc_comp$AUC)

plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(T))   # Adjust time as needed
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(T))  # Adjust time as needed






jpeg("LUAD.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_comp, title = "Time-dependent ROC (Cox vs Weibull)", time = 824, col = "orange", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

lines(roc_weibull$FP[,41],roc_weibull$TP[,41], col = "red", lwd = 2)

lines(roc_cox$FP[,41],roc_cox$TP[,41], col = "blue", lwd = 2)


# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.771)", "Cox (AUC = 0.700)", "Weibull (AUC = 0.700)"), col = c("orange","blue", "red"), lwd = 5,cex=1.5)

title("Time-dependent ROCs at median surival time (824 days)")

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
                       "Age, GenderM, CCNA2, AURKA \nFEN1, LRC4, NLRC4, RASIP1", 
                       "GenderM, AURKB, CCND3, NCALD, MACF1 \nLRC4, NLRC4, PLEKHN1, RASIP1, GPT2,PCOLCE2", 
                       "GenderM, CCNA2, AURKA, AURKB, CCND3 \nNCALD, RASIP1, SPP1, GPT2, SGPL1")



library(ggplot2)
library(reshape2)  # For reshaping data

# Reshape data to long format for ggplot
results_long <- melt(results, id.vars = "Y", variable.name = "Competing_Factors", value.name = "Winning_Probability")





jpeg("WP_time_LUAD.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)
ggplot(results_long, aes(x = exp(Y), y = Winning_Probability, color = Competing_Factors, group = Competing_Factors)) +
  geom_line(size = 1) + 
  geom_point(size = 0.5, alpha = 0.7) +  # Optional: add points
  labs(x = "Time", y = "Winning Probability", title = "Winning probability for uncensored sample", color = "Competing Factors") +
  scale_color_manual(values = c("red", "blue", "forestgreen")) +  # Customize colors
  theme_bw() +  ylim(0,0.85)+# Use theme_bw to have more flexibility
  theme(
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.key.height = unit(1, "cm"),  # Increase vertical spacing between legend items
    legend.position = c(0.7, 0.8),  # Adjust legend position closer to the plot
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

latent_T <- as.data.frame(cbind(data$time_day/result$eta_matrix,result$eta_matrix, data$time_day, data$status))

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

write.csv(latent_T,"latent_T_med_LUAD.csv")

colSums(latent_T[,9:15])





