# Load the library
library(grid)
library(eulerr)

colnames(latent_T_med)[1] <- "ID"
cutoff <- median(c(latent_T_med$T1,latent_T_med$T2,latent_T_med$T3))

latent_T_med$ID <- as.character(latent_T_med$ID)

A <- unlist(latent_T_med[which(latent_T_med$T1<cutoff),1])
B <- unlist(latent_T_med[which(latent_T_med$T2<cutoff),1])
C <- unlist(latent_T_med[which(latent_T_med$T3<cutoff),1])

all_samples <- unique(latent_T_med$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_ADNI.jpeg", quality = 100, units = "in", width = 10, height = 15, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2.5),
     quantities = list(font = 1, cex = 3),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.18,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()

library(survival)
library(survminer)
latent_T_med$PI_group <- ifelse(latent_T_med$G1==1|latent_T_med$G2==1|latent_T_med$G3==1,"Single CF",
                                       ifelse(latent_T_med$G12==1|latent_T_med$G23==1|latent_T_med$G13==1|latent_T_med$G123==1,"Multiple CFs",NA))


fit <- survfit(Surv(time, status) ~ PI_group, data = latent_T_med)

ggsurvplot(fit,
           data = latent_T_med,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = c("Single CF", "Multiple CFs"),  # explicit labels
           xlab = "Time (day)",
           ylab = "Survival Probability",
           palette = "Set1")


jpeg("KM_ADNI.jpeg", quality = 100, units = "in", width = 14, height = 12, res = 300)
ggsurvplot(
  fit,
  data = latent_T_med,
  risk.table = TRUE,
  pval = TRUE,
  legend.labs = c("Multiple CFs","Single CF"),  # explicit labels
  xlab = "Time (day)",
  ylab = "Survival Probability",
  palette = "Set1",
  #ggtheme = theme_bw(base_size = 16),    # increase base font size for everything
  font.x = 30,                           # x-axis font
  font.y = 30,                           # y-axis font
  font.tickslab = 30,                    # tick labels
  font.legend = 30,
  pval.size = 10,  # Increase p-value font size (default is ~5)
  
  # For risk table
  risk.table.fontsize = 10,  # Increase numbers in risk table
  risk.table.height = 0.25, # Adjust table height if needed
  risk.table.y.text = TRUE, # Keep group labels on y-axis
  risk.table.y.text.col = TRUE, # Color group labels
  
  tables.theme = theme(
    axis.text.y = element_text(size = 30),    # group labels
    axis.text.x = element_text(size = 30),    # time labels
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 30),     # numbers in the table
    plot.title = element_text(size = 30)      # risk table title if any   # numbers in the table
  )
)

dev.off()

#brcancer

colnames(latent_T_med_brcancer)[1] <- "ID"
cutoff <- median(c(latent_T_med_brcancer$T1,latent_T_med_brcancer$T2,latent_T_med_brcancer$T3))

latent_T_med_brcancer$ID <- as.character(latent_T_med_brcancer$ID)

A <- unlist(latent_T_med_brcancer[which(latent_T_med_brcancer$T1<cutoff),1])
B <- unlist(latent_T_med_brcancer[which(latent_T_med_brcancer$T2<cutoff),1])
C <- unlist(latent_T_med_brcancer[which(latent_T_med_brcancer$T3<cutoff),1])

all_samples <- unique(latent_T_med_brcancer$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = C, CF2 = A, CF3 = B))

jpeg("Venn_brcancer.jpeg", quality = 100, units = "in", width = 10, height = 12, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2.5),
     quantities = list(font = 1, cex = 3),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.08,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()



#coma

colnames(latent_T_med_coma)[1] <- "ID"
cutoff <- median(c(latent_T_med_coma$T1,latent_T_med_coma$T2,latent_T_med_coma$T3))

latent_T_med_coma$ID <- as.character(latent_T_med_coma$ID)

A <- unlist(latent_T_med_coma[which(latent_T_med_coma$T1<cutoff),1])
B <- unlist(latent_T_med_coma[which(latent_T_med_coma$T2<cutoff),1])
C <- unlist(latent_T_med_coma[which(latent_T_med_coma$T3<cutoff),1])

all_samples <- unique(latent_T_med_coma$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_coma.jpeg", quality = 100, units = "in", width = 14, height = 18, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2),
     quantities = list(font = 1, cex = 2),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()


#COPD

colnames(latent_T_med_COPD)[1] <- "ID"
cutoff <- median(c(latent_T_med_COPD$T1,latent_T_med_COPD$T2,latent_T_med_COPD$T3))

latent_T_med_COPD$ID <- as.character(latent_T_med_COPD$ID)

A <- unlist(latent_T_med_COPD[which(latent_T_med_COPD$T1<cutoff),1])
B <- unlist(latent_T_med_COPD[which(latent_T_med_COPD$T2<cutoff),1])
C <- unlist(latent_T_med_COPD[which(latent_T_med_COPD$T3<cutoff),1])

all_samples <- unique(latent_T_med_COPD$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_COPD.jpeg", quality = 100, units = "in", width = 14, height = 18, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2),
     quantities = list(font = 1, cex = 2),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()



#coma

colnames(latent_T_med_coma)[1] <- "ID"
cutoff <- median(c(latent_T_med_coma$T1,latent_T_med_coma$T2,latent_T_med_coma$T3))

latent_T_med_coma$ID <- as.character(latent_T_med_coma$ID)

A <- unlist(latent_T_med_coma[which(latent_T_med_coma$T1<cutoff),1])
B <- unlist(latent_T_med_coma[which(latent_T_med_coma$T2<cutoff),1])
C <- unlist(latent_T_med_coma[which(latent_T_med_coma$T3<cutoff),1])

all_samples <- unique(latent_T_med_coma$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_COPD.jpeg", quality = 100, units = "in", width = 14, height = 18, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2),
     quantities = list(font = 1, cex = 2),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()


#cancer

colnames(latent_T_med_cancer)[1] <- "ID"
cutoff <- median(c(latent_T_med_cancer$T1,latent_T_med_cancer$T2,latent_T_med_cancer$T3))

latent_T_med_cancer$ID <- as.character(latent_T_med_cancer$ID)

A <- unlist(latent_T_med_cancer[which(latent_T_med_cancer$T1<cutoff),1])
B <- unlist(latent_T_med_cancer[which(latent_T_med_cancer$T2<cutoff),1])
C <- unlist(latent_T_med_cancer[which(latent_T_med_cancer$T3<cutoff),1])

all_samples <- unique(latent_T_med_cancer$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_cancer.jpeg", quality = 100, units = "in", width = 14, height = 18, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2),
     quantities = list(font = 1, cex = 2),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()



## sepsis
colnames(latent_T_med_sepsis)[1] <- "ID"
cutoff <- median(c(latent_T_med_sepsis$T1,latent_T_med_sepsis$T2,latent_T_med_sepsis$T3))

latent_T_med_sepsis$ID <- as.character(latent_T_med_sepsis$ID)

A <- unlist(latent_T_med_sepsis[which(latent_T_med_sepsis$T1<cutoff),1])
B <- unlist(latent_T_med_sepsis[which(latent_T_med_sepsis$T2<cutoff),1])
C <- unlist(latent_T_med_sepsis[which(latent_T_med_sepsis$T3<cutoff),1])

all_samples <- unique(latent_T_med_sepsis$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_sepsis.jpeg", quality = 100, units = "in", width = 12, height = 15, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 3),
     quantities = list(font = 1, cex = 3.5),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()
library(survival)
library(survminer)
latent_T_med_sepsis$PI_group <- ifelse(latent_T_med_sepsis$G1==1|latent_T_med_sepsis$G2==1|latent_T_med_sepsis$G3==1,"Single CF",
                                       ifelse(latent_T_med_sepsis$G12==1|latent_T_med_sepsis$G23==1|latent_T_med_sepsis$G13==1|latent_T_med_sepsis$G123==1,"Multiple CFs",NA))


fit <- survfit(Surv(time, status) ~ PI_group, data = latent_T_med_sepsis)

ggsurvplot(fit,
           data = latent_T_med_sepsis,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = c("Single CF", "Multiple CFs"),  # explicit labels
           xlab = "Time (day)",
           ylab = "Survival Probability",
           palette = "Set1")


jpeg("KM_sepsis.jpeg", quality = 100, units = "in", width = 14, height = 12, res = 300)
ggsurvplot(
  fit,
  data = latent_T_med_sepsis,
  risk.table = TRUE,
  pval = TRUE,
  legend.labs = c("Multiple CFs","Single CF"),  # explicit labels
  xlab = "Time (day)",
  ylab = "Survival Probability",
  palette = "Set1",
  #ggtheme = theme_bw(base_size = 16),    # increase base font size for everything
  font.x = 30,                           # x-axis font
  font.y = 30,                           # y-axis font
  font.tickslab = 30,                    # tick labels
  font.legend = 30,
  pval.size = 10,  # Increase p-value font size (default is ~5)
  
  # For risk table
  risk.table.fontsize = 10,  # Increase numbers in risk table
  risk.table.height = 0.25, # Adjust table height if needed
  risk.table.y.text = TRUE, # Keep group labels on y-axis
  risk.table.y.text.col = TRUE, # Color group labels
  
  tables.theme = theme(
    axis.text.y = element_text(size = 30),    # group labels
    axis.text.x = element_text(size = 30),    # time labels
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 30),     # numbers in the table
    plot.title = element_text(size = 30)      # risk table title if any   # numbers in the table
  )
)

dev.off()



## sepsis
colnames(latent_T_med_LUAD)[1] <- "ID"
cutoff <- median(c(latent_T_med_LUAD$T1,latent_T_med_LUAD$T2,latent_T_med_LUAD$T3))

latent_T_med_LUAD$ID <- as.character(latent_T_med_LUAD$ID)

A <- unlist(latent_T_med_LUAD[which(latent_T_med_LUAD$T1<cutoff),1])
B <- unlist(latent_T_med_LUAD[which(latent_T_med_LUAD$T2<cutoff),1])
C <- unlist(latent_T_med_LUAD[which(latent_T_med_LUAD$T3<cutoff),1])

all_samples <- unique(latent_T_med_LUAD$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_LUAD.jpeg", quality = 100, units = "in", width = 12, height = 15, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 3),
     quantities = list(font = 1, cex = 3.5),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()

library(survival)
library(survminer)
latent_T_med_LUAD$PI_group <- ifelse(latent_T_med_LUAD$G1==1|latent_T_med_LUAD$G2==1|latent_T_med_LUAD$G3==1,"Single CF",
                                       ifelse(latent_T_med_LUAD$G12==1|latent_T_med_LUAD$G23==1|latent_T_med_LUAD$G13==1|latent_T_med_LUAD$G123==1,"Multiple CFs",NA))

latent_T_med_LUAD$PI_23 <- ifelse(latent_T_med_LUAD$G2==1,"G2",
                                     ifelse(latent_T_med_LUAD$G3==1,"G3",NA))

fit <- survfit(Surv(time, status) ~ PI_group, data = latent_T_med_LUAD)

ggsurvplot(fit,
           data = latent_T_med_LUAD,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = c("Single CF", "Multiple CFs"),  # explicit labels
           xlab = "Time",
           ylab = "Survival Probability",
           palette = "Set1")


jpeg("KM_LUAD.jpeg", quality = 100, units = "in", width = 14, height = 12, res = 300)
ggsurvplot(
  fit,
  data = latent_T_med_LUAD,
  risk.table = TRUE,
  pval = TRUE,
  legend.labs = c("Multiple CFs","Single CF"),  # explicit labels
  xlab = "Time (day)",
  ylab = "Survival Probability",
  palette = "Set1",
  #ggtheme = theme_bw(base_size = 16),    # increase base font size for everything
  font.x = 30,                           # x-axis font
  font.y = 30,                           # y-axis font
  font.tickslab = 30,                    # tick labels
  font.legend = 30,
  pval.size = 10,  # Increase p-value font size (default is ~5)
  
  # For risk table
  risk.table.fontsize = 10,  # Increase numbers in risk table
  risk.table.height = 0.25, # Adjust table height if needed
  risk.table.y.text = TRUE, # Keep group labels on y-axis
  risk.table.y.text.col = TRUE, # Color group labels
  
  tables.theme = theme(
    axis.text.y = element_text(size = 30),    # group labels
    axis.text.x = element_text(size = 30),    # time labels
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 30),     # numbers in the table
    plot.title = element_text(size = 30)      # risk table title if any   # numbers in the table
  )
)

dev.off()










## HCC
colnames(latent_T_med_HCC)[1] <- "ID"
cutoff <- median(c(latent_T_med_HCC$T1,latent_T_med_HCC$T2))

latent_T_med_HCC$ID <- as.character(latent_T_med_HCC$ID)

A <- unlist(latent_T_med_HCC[which(latent_T_med_HCC$T1<cutoff),1])
B <- unlist(latent_T_med_HCC[which(latent_T_med_HCC$T2<cutoff),1])

all_samples <- unique(latent_T_med_HCC$ID)
in_any <- unique(c(A, B))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B))

jpeg("Venn_HCC.jpeg", quality = 100, units = "in", width = 14, height = 18, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange"), alpha = 0.5),
     labels = list(font = 1, cex = 2),
     quantities = list(font = 1, cex = 2),
     main = list(label = c("Categorisation of uncensored samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()



## cox max
colnames(latent_risk_med)[1] <- "ID"
cutoff <- median(c(latent_risk_med$lp1,latent_risk_med$lp2,latent_risk_med$lp3))

latent_risk_med$ID <- as.character(latent_risk_med$ID)

A <- unlist(latent_risk_med[which(latent_risk_med$lp1<cutoff),1])
B <- unlist(latent_risk_med[which(latent_risk_med$lp2<cutoff),1])
C <- unlist(latent_risk_med[which(latent_risk_med$lp3<cutoff),1])

all_samples <- unique(latent_risk_med$ID)
in_any <- unique(c(A, B, C))
complement <- setdiff(all_samples, in_any)
complement_size <- length(complement)

# Step 1: Install and load the eulerr package


fit <- euler(list(CF1 = A, CF2 = B, CF3 = C))

jpeg("Venn_sepsis_coxmax.jpeg", quality = 100, units = "in", width = 14, height = 18, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "orange", "seagreen"), alpha = 0.5),
     labels = list(font = 1, cex = 2),
     quantities = list(font = 1, cex = 2),
     main = list(label = c("Categorisation of samples\n"), cex = 3)
)

grid.text(paste("Complement size:", complement_size),
          x = 0.8, y = 0.06,  # normalized coordinates (0–1)
          gp = gpar(fontsize = 22))

dev.off()

