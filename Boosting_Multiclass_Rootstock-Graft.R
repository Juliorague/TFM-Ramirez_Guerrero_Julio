library(tidyverse)
library(tidymodels)
library(vip)

# read my data with the read_csv2 function which by default is designed for European data where the 
# delimiter is a ";" and decimals are indicated by ","
metab <- read_csv2("Metabolite_data.csv")

# transform the dataset from wide to long format. All variables related to individuals are grouped in a 
# single column called individual, while abundance values are placed in another column called abundance
metab_1 <- metab %>%
  gather(key = "individual", value = "expression", -c(identifier))

# add the new column sample_type_origin
metab_1 <- metab_1 %>%
  mutate(sample_type_origin = case_when(
    str_detect(individual, "A2") & str_detect(individual, "N") ~ "tolerant_N",
    str_detect(individual, "A2") & str_detect(individual, "R(?!S)") ~ "tolerant_R", # (?!S) ensures that it does not follow a S
    str_detect(individual, "A2") & str_detect(individual, "RS|SS|S") ~ "tolerant_S",
    str_detect(individual, "B16") & str_detect(individual, "N") ~ "sensitive_N",
    str_detect(individual, "B16") & str_detect(individual, "R(?!S)") ~ "sensitive_R", 
    str_detect(individual, "B16") & str_detect(individual, "RS|SS|S") ~ "sensitive_S",
  ))

# add the new column rootstock_graft
metab_1 <- metab_1 %>%
  mutate(rootstock_graft = case_when(
    str_detect(individual, "A2") & str_detect(individual, "Gal") ~ "tolerant_Galicia",
    str_detect(individual, "A2") & str_detect(individual, "homo") ~ "tolerant_homo", 
    str_detect(individual, "A2") & str_detect(individual, "sin") ~ "tolerant_absent",
    str_detect(individual, "B16") & str_detect(individual, "Gal") ~ "sensitive_Galicia",
    str_detect(individual, "B16") & str_detect(individual, "homo") ~ "sensitive_homo", 
    str_detect(individual, "B16") & str_detect(individual, "sin") ~ "sensitive_absent",
  ))

# transform the dataset from long to wide format
metab_1 <- metab_1 %>%
  pivot_wider(names_from = identifier, values_from = expression)

## data splitting
metab_dt <- select(metab_1, -c(individual))
metab_dt$sample_type_origin <- factor(metab_dt$sample_type_origin)
metab_dt$rootstock_graft <- factor(metab_dt$rootstock_graft)

metab_split <- initial_split(metab_dt, strata = sample_type_origin, prop = 0.75)
metab_train <- training(metab_split) %>%
  select(-sample_type_origin)
metab_test <- testing(metab_split) %>%
  select(-sample_type_origin)

## preprocessing
preprocessing_recipe <-
  recipes::recipe(rootstock_graft ~ ., data = metab_train) %>%
  step_corr(all_predictors(), threshold = 0.9) %>% ## remove collinear variables
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_normalize(all_numeric(), -all_outcomes()) %>%
  step_impute_knn(all_numeric(), neighbors = 5) %>%
  prep()

## k-fold cross-validation for tuning
metab_cv <- vfold_cv(metab_train, v=5, repeats = 5, strata = rootstock_graft)

# XGBoost model specification
xgboost_model <-
  boost_tree(
    mode = "classification",
    trees = 100, # B parameter
    min_n = tune(),
    tree_depth = tune(), # d parameter
    learn_rate = tune()
  ) %>%
  set_engine("xgboost", objective = "multi:softprob", lambda=0, alpha=1, verbose=0)

# grid specification
xgboost_params <-
  parameters(
    min_n(),
    tree_depth(),
    learn_rate()
  )

xgboost_grid <-
  grid_max_entropy(
    xgboost_params,
    size = 50
  )

print(xgboost_grid)

## workflow
xgboost_wf <-
  workflows::workflow() %>%
  add_model(xgboost_model) %>%
  add_formula(rootstock_graft ~ .)

# hyperparameter tuning
xgboost_tuned <- tune_grid(
  object = xgboost_wf,
  resamples = metab_cv,
  grid = xgboost_grid,
  # metrics = yardstick::metric_set(rmse, rsq, mae),
  control = control_grid(verbose = FALSE)
)

## explore tuning results
collect_metrics(xgboost_tuned)

options(repr.plot.width=14, repr.plot.height=8)

xgboost_tuned %>%
  collect_metrics() %>%
  filter(.metric == "accuracy") %>%
  dplyr::select(mean, min_n:learn_rate) %>%
  pivot_longer(min_n:learn_rate,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "accuracy")

## select and evaluate the best model
xgboost_tuned %>%
  show_best(metric = "roc_auc")

xgboost_best_params <- xgboost_tuned %>%
  select_best(metric = "roc_auc")         

print(xgboost_best_params)

final_xgb <- finalize_workflow(
  xgboost_wf, # built workflow
  xgboost_best_params # selected best model after finetuning
)

final_res <- last_fit(final_xgb, metab_split)
collect_metrics(final_res)

collect_predictions(final_res) %>%
  metrics(rootstock_graft, .pred_class)

## confusion matrix
cm <- collect_predictions(final_res) %>%
  conf_mat(rootstock_graft, .pred_class)

print(cm)

autoplot(cm, type="heatmap")

## variable importance
variable_importance <- final_xgb %>%
  fit(data = juice(preprocessing_recipe)) %>%
  extract_fit_parsnip() %>%                   
  vip(num_features = 10)

variable_importance

importance_data <- variable_importance$data