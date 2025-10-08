load("D:/Lily/Ling/Cox-max/GSE.RData")
source("D:/Lily/Ling/Cox-max/Functions.txt")
library(survival)
selected_vars <- c("GenderM", "Age", "status", "time", "ID_REF",
                   "11725985_s_at", "11726583_s_at", "11757883_s_at", "11761102_x_at")
GSE65682 <- GSE[, selected_vars]
str(GSE65682)
head(GSE65682)

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

## cox/weibull
cox_m <- coxph(Surv(time, status) ~ Age + GenderM + NONO + CKAP4 + RNF4 + FCAR, data = data)
summary(cox_m)
concordance(cox_m) #n= 468 Concordance= 0.6002 se= 0.02699

w_m <- survreg(Surv(time, status) ~ Age + GenderM + NONO + CKAP4 + RNF4 + FCAR, data = data)
summary(w_m)
concordance(w_m) #Concordance= 0.6002 se= 0.02701

#########################
X = cbind(data$Age, data$GenderM, 
          data$NONO, data$CKAP4, 
          data$RNF4, data$FCAR)
label = list(c(1), c(1),
             c(1,2,3), c(2,3), 
             c(2,3), c(1,2,3)) # $cindex [1] 0.6612254
> result$beta
[[1]]
[1]  0.43164937  0.89421338 -0.02579620 -0.04832729

[[2]]
[1]  0.76307341 -0.09637643 -1.28631730  0.13543773

[[3]]
[1] -0.7345900  0.5036846  0.5475775  0.2846583

label = list(c(1), c(1),
             c(2,3), c(3), 
             c(2,3), c(3)) # $cindex [1] 0.6581749

label = list(c(1), c(1),
             c(2,3), c(3), 
             c(2,3), c(1,3)) # $cindex [1] 0.6584458

label = list(c(1), c(1),
             c(1,2,3), c(3), 
             c(2,3), c(1,3)) # $cindex [1] 0.6587049

label = list(c(1), c(1),
             c(1,2,3), c(2,3), 
             c(2,3), c(1,3)) # $cindex [1] 0.6587049

label = list(c(1), c(1),
             c(1,2,3), c(2,3), 
             c(2,3), c(1,2,3)) # $cindex [1] 0.6592938

label = list(c(1), c(1),
             c(1,2,3), c(2,3), 
             c(2,3), c(1,2,3)) # $cindex [1] 0.6612254

transformed_data = data_transform(time=data$time, status=data$status, X, label)
init_param = param_init(X, label)

result <- main_iter(transformed_data$time, transformed_data$status, init_beta = init_param$beta,
                    transformed_data$X, maxit = 1e2, stepsize = 0.8, tol = 5e-2, method = "BFGS", copy = FALSE)
compute_cind_parlik(time=data$time, status=data$status, transformed_data$X, beta = result$beta, group = result$group, component = "both")

result$beta
#-----------------------------------------------
# Bootstrap for Cox-Max 模型系数：点估计 + SE + 5%–95% 置信区间
#-----------------------------------------------

n_obs <- nrow(data)
# —— 1. Bootstrap 参数 —— 
B      <- 2            # 重抽样次数，可视资源调整
set.seed(2025)

# 存储每次 bootstrap 的系数
coef_mat <- matrix(NA, nrow = B, ncol = 17)

# —— 2. 并行/非并行循环都可以 —— 
for (b in seq_len(B)) {
  # (1) 有放回重采样
  idx <- sample.int(n_obs, n_obs, replace = TRUE)
  df_b <- data[idx, ]
  
  # (2) 拟合 Cox PH
  X = cbind(df_b$age, df_b$sex_male, df_b$slos,
            df_b$scoma, df_b$wblc, df_b$pafi, 
            df_b$sps, df_b$aps, df_b$adlsc, df_b$hrt
  )
  
  label = list(c(1,2), c(1,2), c(1,2),
               c(1,3), c(3), c(1), 
               c(1,2,3), c(2), c(3), c(1,2))
  
  transformed_data = data_transform(time=df_b$d.time, status=df_b$death, X, label)
  init_param = param_init(X, label)
  
  fit_b <- main_iter(transformed_data$time, transformed_data$status, init_beta = init_param$beta,
                     transformed_data$X, maxit = 1e2, stepsize = 0.8, tol = 5e-2, method = "BFGS", copy = FALSE)
  # (3) 抽取系数
  coef_mat[b, ] <- unlist(fit_b$beta)
}

# —— 3. 计算 Bootstrap SE & 95% CI —— 
se_boot <- apply(coef_mat, 2, sd, na.rm = TRUE)

ci_mat <- t(apply(coef_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
colnames(ci_mat) <- c("lower95", "upper95")

# —— 4. 汇总结果 —— 
results <- data.frame(
  se_boot   = se_boot,
  lower95   = ci_mat[, "lower95"],
  upper95   = ci_mat[, "upper95"]
)

# （可选）直接计算 HR 及其 CI
results$HR_lower <- exp(results$lower95)
results$HR_upper <- exp(results$upper95)

# 打印
print(round(results, 3))


