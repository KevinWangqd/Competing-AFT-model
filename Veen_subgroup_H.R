data <- read.csv("latent_T_med.csv")
str(data)

colnames(data)
library(dplyr)

AD_data <- data %>%
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

str(AD_data)
vars <- c("time", "APOE4", "PTGENDER", "PTEDUCAT",
          "AGE", "ADAS13", "RAVLT_immediate", "RAVLT_learning",
          "FAQ", "Hippocampus", "ICV", "MidTemp")

# mark categorical variables
factorVars <- c("APOE4", "PTGENDER")

library(tableone)

tab1 <- CreateTableOne(vars = vars, strata = "Group", data = AD_data, factorVars = factorVars)

# Print with all levels shown
print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
library(dplyr)

library(tidyr)
library(ggplot2)

AD_data <- data %>%
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
    Group = factor(Group, levels = c("G1","G2","G3","G12","G23","G13","G123")) # <- drop "Other" here
  ) %>%
  filter(Group != "Other")   # <- remove "Other" rows completely
# separate numeric and categorical vars
num_vars <- setdiff(vars, factorVars)


cat_vars <- c("APOE4", "PTGENDER")

# 1. Numeric: compute group means
df_num <- AD_data %>%
  group_by(Group) %>%
  summarise(across(all_of(num_vars), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-Group, names_to = "Variable", values_to = "Value")

# 2. Categorical: compute percentages for each level within group
df_cat <- AD_data %>%
  pivot_longer(all_of(cat_vars), names_to = "Variable", values_to = "Level") %>%
  group_by(Group, Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Group, Variable) %>%
  mutate(Value = N / sum(N) * 100) %>%   # percentage
  ungroup() %>%
  unite("Variable", Variable, Level, sep = "_")  # e.g., APOE4_0, APOE4_1, PTGENDER_0, PTGENDER_1

# 3. Combine numeric and categorical summaries
df_combined <- bind_rows(
  df_num %>% select(Group, Variable, Value),
  df_cat %>% select(Group, Variable, Value)
)

# 4. Scale all numeric values to 0–1 range (min-max within variable)
df_scaled <- df_combined %>%
  group_by(Variable) %>%
  mutate(ScaledValue = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 5. Complete grid for missing combinations
all_groups <- unique(AD_data$Group)
all_vars <- unique(df_scaled$Variable)
full_grid <- expand.grid(Group = all_groups, Variable = all_vars)

df_heatmap <- full_grid %>%
  left_join(df_scaled %>% select(Group, Variable, ScaledValue), by = c("Group","Variable")) %>%
  mutate(ScaledValue = replace_na(ScaledValue, 0))

var_labels <- c(
  "APOE4_0"       = "APOEε4 (0)",
  "APOE4_1"       = "APOEε4 (1)",
  "APOE4_2"       = "APOEε4 (2)",
  "ADAS13"        = "ADAS13",
  "RAVLT_immediate" = "RAVLT (immediate)",
  "RAVLT_learning"  = "RAVLT (learning)",
  "FAQ"           = "FAQ",
  "PTGENDER_1"    = "GenderM",
  "PTGENDER_0"    = "GenderF",
  "PTEDUCAT"      = "Education",
  "AGE"           = "Age",
  "Hippocampus"   = "Hippocampus",
  "ICV"           = "ICV",
  "MidTemp"       = "MidTemp",
  "time"          = "Survival time"
)

# Apply renaming
df_heatmap <- df_heatmap %>%
  mutate(Variable = recode(Variable, !!!var_labels))

var_order <- rev(c(  "Survival time",
  "APOEε4 (0)","APOEε4 (1)","APOEε4 (2)",
  "ADAS13",
  "RAVLT (immediate)","RAVLT (learning)",
  "FAQ",
  "GenderM","GenderF",
  "Education",
  "Age",
  "Hippocampus",
  "ICV",
  "MidTemp"
))

df_heatmap <- df_heatmap %>%
  mutate(Variable = factor(Variable, levels = var_order))

# Plot again
ggplot(df_heatmap, aes(x = Group, y = Variable, fill = ScaledValue)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("navy","skyblue","white","pink","red")) +
  labs(x = "Group", y = "Variable", fill = "Scaled 0–1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# 6. Plot heatmap
jpeg("Heatmap_ADNI.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)

ggplot(df_heatmap, aes(x = Group, y = Variable, fill = ScaledValue)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("navy", "skyblue", "white","pink", "red"))+
  labs(title = "",
       x = "Group", y = "Variable", fill = "Relative Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
