library(tidymodels)
library(tidyverse)

Variable_importance <- read_csv2("Variable_Importance.csv")



## Lasso-penalised Logistic Regression

# Binary classification: Rootstock type
Lasso_Rt_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Lasso")) %>%
  filter(Classification_type %in% c("Rootstock_type"))

top_5_Lasso_Rt_vi <- Lasso_Rt_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Binary classification: Homo vs. Absent
Lasso_HA_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Lasso")) %>%
  filter(Classification_type %in% c("Homo_vs_Absent"))

top_5_Lasso_HA_vi <- Lasso_HA_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Binary classification: Homo vs. Galicia
Lasso_HG_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Lasso")) %>%
  filter(Classification_type %in% c("Homo_vs_Galicia"))

top_5_Lasso_HG_vi <- Lasso_HG_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Multiclass classification: Rootstock-Graft
Lasso_RG_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Lasso")) %>%
  filter(Classification_type %in% c("Rootstock_Graft"))

top_5_Lasso_RG_vi <- Lasso_RG_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)



## Random Forest

# Binary classification: Rootstock type
RandomF_Rt_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Random_Forest")) %>%
  filter(Classification_type %in% c("Rootstock_type"))

top_5_RandomF_Rt_vi <- RandomF_Rt_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Binary classification: Homo vs. Absent
RandomF_HA_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Random_Forest")) %>%
  filter(Classification_type %in% c("Homo_vs_Absent"))

top_5_RandomF_HA_vi <- RandomF_HA_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Binary classification: Homo vs. Galicia
RandomF_HG_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Random_Forest")) %>%
  filter(Classification_type %in% c("Homo_vs_Galicia"))

top_5_RandomF_HG_vi <- RandomF_HG_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Multiclass classification: Rootstock-Graft
RandomF_RG_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Random_Forest")) %>%
  filter(Classification_type %in% c("Rootstock_Graft"))

top_5_RandomF_RG_vi <- RandomF_RG_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)



## Boosting

# Binary classification: Rootstock type
Boosting_Rt_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Boosting")) %>%
  filter(Classification_type %in% c("Rootstock_type"))

top_5_Boosting_Rt_vi <- Boosting_Rt_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Binary classification: Homo vs. Absent
Boosting_HA_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Boosting")) %>%
  filter(Classification_type %in% c("Homo_vs_Absent"))

top_5_Boosting_HA_vi <- Boosting_HA_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Binary classification: Homo vs. Galicia
Boosting_HG_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Boosting")) %>%
  filter(Classification_type %in% c("Homo_vs_Galicia"))

top_5_Boosting_HG_vi <- Boosting_HG_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)

# Multiclass classification: Rootstock-Graft
Boosting_RG_vi <- Variable_importance %>%
  filter(Supervised_Learning_model %in% c("Boosting")) %>%
  filter(Classification_type %in% c("Rootstock_Graft"))

top_5_Boosting_RG_vi <- Boosting_RG_vi %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(Importance, na.rm = TRUE)) %>%
  slice_max(mean_importance, n = 5)