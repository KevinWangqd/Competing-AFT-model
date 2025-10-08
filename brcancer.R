library(survival)
library(casebase)
data("brcancer")
attach(brcancer)
str(brcancer)



mean(brcancer$estrec)
min(brcancer$estrec)
max(brcancer$estrec)


sum(brcancer$tgrade3)

brcancer$age <- ifelse(brcancer$age>45,1,0)

brcancer$tgrade2 <- ifelse(brcancer$tgrade=="II",1,0)
brcancer$tgrade3 <- ifelse(brcancer$tgrade=="III",1,0)

cox_m <- coxph(Surv(time,cens)~ hormon + meno + age + tsize+ tgrade2 + tgrade3 + pnodes+ 
                 progrec + estrec, data=brcancer)

summary(cox_m)
concordance(cox_m) 

w_m <- survreg(Surv(time,cens)~ hormon + meno + age + tsize+ tgrade2 + tgrade3 + pnodes+ 
                 progrec + estrec, data=brcancer)
summary(w_m)
concordance(w_m) 


source("Functions.txt")
#source("Functions_new_ln.txt")
data <- brcancer
#data$age <- scale(data$age)
data$tsize <- scale(data$tsize)
data$pnodes <- scale(data$pnodes)
data$progrec <- scale(data$progrec)
data$estrec <- scale(data$estrec)


w_m <- survreg(Surv(time,cens)~ hormon + meno + age + tsize+ tgrade2 + tgrade3 + pnodes+ 
                 progrec + estrec, data=data , dist="weibull")



summary(w_m)
concordance(w_m) 

cox_m <- coxph(Surv(time,cens)~ hormon + meno + age + tsize+ tgrade2 + tgrade3 + pnodes+ 
                 progrec + estrec, data=brcancer)

summary(cox_m)
concordance(cox_m) 


Y <- log(data$time)
delta <- data$cens

X = cbind(data$hormon, data$meno, data$age, 
          data$tsize, data$tgrade2, data$tgrade3, 
          data$pnodes, data$progrec, data$estrec)

#label = list(c(1), c(2), c(2), 
#             c(1,2), c(1), c(1), 
#             c(3), c(3), c(2,3))

label = list(c(1), c(2), c(1,2), 
             c(1,2), c(1), c(1), 
             c(3), c(3), c(2))


transformed_data = data_transform(Y, delta, X, label)
#transformed_data
# with defaul initialization
init_param = param_init(X, label)

#init_param

result <- EM_fit(transformed_data, init_param, lambda1 = 1, lambda2 = 0.2, maxit = 10^3, tolerance = 1e-6)

result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result

expected_time = sur_time_all(result, transformed_data,method = "MR", M = 300, n_samples = 1000)
#expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^3)

surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)

# Print the result
print(c_index$concordance)







library(timeROC)
library(survival)
time_points <- quantile(data$time, probs = seq(0.1, 0.9, by = 0.1))

# 1. For Cox model, get the predicted risk scores (linear predictor)
cox_risk_scores <- predict(cox_m, type = "risk")

# 3. Time-dependent ROC for the Cox model
roc_cox <- timeROC(T = data$time, delta = data$cens, marker = cox_risk_scores, cause = 1, times = time_points)
#roc_cox <- timeROC(T = T, delta = delta, marker = cox_m$linear.predictors, cause = 1, times = time_points)
mean(roc_cox$AUC)
# 5. Plot the ROC curves for both models
plot(roc_cox, main = "Cox Model Time-dependent ROC", time = median(data$time))   # Adjust time as needed

w_risk_scores <- -predict(w_m, type = "lp")

roc_weibull <- timeROC(T = data$time, delta = data$cens, marker = w_risk_scores, cause = 1, times = time_points)
plot(roc_weibull, main = "Weibull Model Time-dependent ROC", time = median(data$time))  # Adjust time as needed
mean(roc_weibull$AUC)

roc_comp <- timeROC(T = data$time, delta = data$cens, marker = -expected_time, cause = 1, times = time_points)
plot(roc_comp, main = "Weibull Model Time-dependent ROC", time = median(data$time))  # Adjust time as needed

mean(roc_comp$AUC)
roc_weibull$AUC
roc_comp$AUC

T = data$time
median(T)

jpeg("brcancer.jpeg", quality = 100, units = "in", width = 8, height = 8, res = 300)
plot(roc_comp, title = "Time-dependent ROC (Cox vs Weibull)", time = 515, col = "orange", lwd = 2, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

lines(roc_weibull$FP[,2],roc_weibull$TP[,2], col = "red", lwd = 2, lty = 2)

lines(roc_cox$FP[,2],roc_cox$TP[,2], col = "blue", lwd = 2)


# 4. Add a legend to differentiate between the models
legend("bottomright", legend = c("Competing Weibull (AUC = 0.78)", "Cox (AUC = 0.75)", "Weibull (AUC = 0.75)"), 
       col = c("orange","blue", "red"), lwd = 5, cex=1.5)

title("Time-dependent ROCs at 20% quantile surival time (515)")

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
  results$mean_X1[i] <- mean(output_e$X3, na.rm = TRUE)
  results$mean_X2[i] <- mean(output_e$X1, na.rm = TRUE)
  results$mean_X3[i] <- mean(output_e$X2, na.rm = TRUE)
}

# Print or visualize results

mean(results$mean_X1)
mean(results$mean_X2)
mean(results$mean_X3)



colnames(results) <- c("Y", "Pnodes, Progrec",
                       "HormonT, Age45, Tsize, Tgrade", 
                       "Monestat, Age45, Tsize, Estrec"
                       
)




library(ggplot2)
library(reshape2)  # For reshaping data

# Reshape data to long format for ggplot
results_long <- melt(results, id.vars = "Y", variable.name = "Competing_Factors", value.name = "Winning_Probability")



jpeg("WP_time_brcancer.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)
ggplot(results_long, aes(x = exp(Y), y = Winning_Probability, color = Competing_Factors, group = Competing_Factors)) +
  geom_line(size = 1) + 
  geom_point(size = 0.5, alpha = 0.7) +  # Optional: add points
  labs(x = "Time", y = "Winning Probability", title = "Winning probability for uncensored sample", color = "Competing Factors") +
  scale_color_manual(values = c("red", "blue", "forestgreen", "purple")) +  # Customize colors
  theme_bw() +  # Use theme_bw to have more flexibility
  theme(
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.key.height = unit(0.8, "cm"),  # Increase vertical spacing between legend items
    legend.position = c(0.8, 0.15),  # Adjust legend position closer to the plot
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

latent_T <- as.data.frame(cbind(data$time/result$eta_matrix,result$eta_matrix, data$time, data$cens))

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

write.csv(latent_T,"latent_T_med_brcancer.csv")

colSums(latent_T[,9:15])

