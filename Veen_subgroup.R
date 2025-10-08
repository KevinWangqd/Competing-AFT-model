data1<- read.csv("latent_T_med_LUAD.csv")
GSE72094 <- read.csv("final_data.csv")
data <- GSE72094[,c(1,3,4,11,12,13)]
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


merged_data <- merge(data1, data, by = "X", all.x = TRUE)
data <- merged_data

colnames(data)
library(dplyr)
library(tidyr)
library(ggplot2)

sp_data <- data %>%
  mutate(
    Group = case_when(
      G1   == 1 ~ "G1",
      G2   == 1 ~ "G2",
      G3   == 1 ~ "G3",
      G12  == 1 ~ "G12",
      G23  == 1 ~ "G23",
      G13  == 1 ~ "G13",
      G123 == 1 ~ "G123",
      TRUE      ~ "Other"
    ),
    Group = factor(Group, levels = c("G1","G2","G3","G12","G23","G13","G123","Other"))
  )

str(sp_data)
colnames(sp_data)
vars <- c("time",  "age","GenderM",  "CCNA2",
          "AURKA", "AURKB",      "FEN1",       "CCND3" ,     "NCALD"  ,    "MACF1"   ,   "LRC4"    ,   "NLRC4"    ,  "PLEKHN1"   , "RASIP1" ,    "SPP1"  ,    
           "GPT2" ,      "SGPL1"   ,   "PCOLCE2")

# mark categorical variables
factorVars <- c("GenderM")

library(tableone)

tab1 <- CreateTableOne(vars = vars, strata = "Group", data = sp_data, factorVars = factorVars)

# Print with all levels shown
print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)


sp_data <- sp_data %>%
  filter(Group != "Other")

# Separate numeric and categorical vars
num_vars <- setdiff(vars, factorVars)
cat_vars <- factorVars

# 1. Numeric: compute group means
df_num <- sp_data %>%
  group_by(Group) %>%
  summarise(across(all_of(num_vars), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-Group, names_to = "Variable", values_to = "Value")

# 2. Categorical: compute percentages
df_cat <- sp_data %>%
  pivot_longer(all_of(cat_vars), names_to = "Variable", values_to = "Level") %>%
  group_by(Group, Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Group, Variable) %>%
  mutate(Value = N / sum(N) * 100) %>%
  ungroup() %>%
  unite("Variable", Variable, Level, sep = "_")  # e.g., GenderM_0, GenderM_1

# 3. Combine summaries
df_combined <- bind_rows(
  df_num %>% select(Group, Variable, Value),
  df_cat %>% select(Group, Variable, Value)
)

# 4. Scale to 0–1 within each variable
df_scaled <- df_combined %>%
  group_by(Variable) %>%
  mutate(ScaledValue = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 5. Complete grid for missing combos
all_groups <- unique(sp_data$Group)
all_vars   <- unique(df_scaled$Variable)
full_grid <- expand.grid(Group = all_groups, Variable = all_vars)

df_heatmap <- full_grid %>%
  left_join(df_scaled %>% select(Group, Variable, ScaledValue), by = c("Group","Variable")) %>%
  mutate(ScaledValue = replace_na(ScaledValue, 0))

# 6. Pretty variable labels
var_labels <- c(
  "time"      = "Survival time",
  "age"       = "Age",
  "GenderM_0" = "GenderF",
  "GenderM_1" = "GenderM",
  "CCNA2"     = "CCNA2",
  "AURKA"      = "AURKA",
  "AURKB"      = "AURKB",
  "FEN1"   = "FEN1",
  "CCND3"      = "CCND3",
  "NCALD"    = "NCALD",
  "MACF1"    = "MACF1",
  "LRC4"    = "LRC4",
  "NLRC4"    = "NLRC4",
  "PLEKHN1"    = "PLEKHN1",
  "RASIP1"    = "RASIP1",
  "SPP1"    = "SPP1",
  "GPT2"    = "GPT2",
  "SGPL1"    = "SGPL1",
  "PCOLCE2"    = "PCOLCE2"
)






df_heatmap <- df_heatmap %>%
  mutate(Variable = recode(Variable, !!!var_labels))

# 7. Define variable order (upside-down: Survival time on top, proteins bottom)
var_order <- rev(c(
  "Survival time",
  "Age",
  "GenderM","GenderF",
  "NONO","CKAP4","RNF4","FCAR","PLEKHO1","BMP6","RNASE2"
))

df_heatmap <- df_heatmap %>%
  mutate(Variable = factor(Variable, levels = var_order))

# 8. Plot heatmap
jpeg("Heatmap_LUAD.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)

ggplot(df_heatmap, aes(x = Group, y = Variable, fill = ScaledValue)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("navy","skyblue","white","pink","red")) +
  labs(x = "Group", y = "Variable", fill = "Relative Level (0–1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
