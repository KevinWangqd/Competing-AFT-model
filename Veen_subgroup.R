data <- read.csv("latent_T_med_sepsis.csv")

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
vars <- c( "time"      , "Age",  "GenderM" ,     "NONO"   ,        "CKAP4"    ,      "RNF4"    ,    
           "FCAR"      ,     "PLEKHO1"      ,  "BMP6"     ,      "RNASE2"   ,      "IGHG1"   ,       "IL1R2"    ,     
           "LCN2"     ,      "LTF"   ,         "MMP8"   ,        "OLFM4"    )

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
  "Age"       = "Age",
  "GenderM_0" = "GenderF",
  "GenderM_1" = "GenderM",
  "NONO"      = "NONO",
  "CKAP4"     = "CKAP4",
  "RNF4"      = "RNF4",
  "FCAR"      = "FCAR",
  "PLEKHO1"   = "PLEKHO1",
  "BMP6"      = "BMP6",
  "RNASE2"    = "RNASE2",
  "IGHG1"    = "IGHG1",
  "IL1R2"    = "IL1R2",
  "LCN2"    = "LCN2",
  "LTF"    = "LTF",
  "MMP8"    = "MMP8",
  "OLFM4"    = "OLFM4"
)

df_heatmap <- df_heatmap %>%
  mutate(Variable = recode(Variable, !!!var_labels))

# 7. Define variable order (upside-down: Survival time on top, proteins bottom)
var_order <- rev(c(
  "Survival time",
  "Age",
  "GenderM","GenderF",
  "NONO","CKAP4","RNF4","FCAR","PLEKHO1","BMP6","RNASE2","IGHG1","IL1R2","LCN2", "LTF","MMP8","OLFM4"
))

df_heatmap <- df_heatmap %>%
  mutate(Variable = factor(Variable, levels = var_order))

# 8. Plot heatmap
jpeg("Heatmap_sep.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)

ggplot(df_heatmap, aes(x = Group, y = Variable, fill = ScaledValue)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("navy","skyblue","white","pink","red")) +
  labs(x = "Group", y = "Variable", fill = "Relative Level (0–1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
