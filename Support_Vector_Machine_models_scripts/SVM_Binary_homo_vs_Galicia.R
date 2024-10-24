library(tidyverse)
library(tidymodels)

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

# add the new column graft_type
metab_1 <- metab_1 %>%
  mutate(graft_type = case_when(
    str_detect(individual, "Gal") ~ "Galicia",
    str_detect(individual, "homo") ~ "homo", 
    str_detect(individual, "sin") ~ "absent",
  ))

# filter the data to keep only those individuals with the values I am interested in
homo_vs_Galicia <- metab_1 %>%
  filter(graft_type %in% c("homo", "Galicia"))

# transform the dataset from long to wide format
homo_vs_Galicia <- homo_vs_Galicia %>%
  pivot_wider(names_from = identifier, values_from = expression)

## data splitting
metab_dt <- select(homo_vs_Galicia, -c(individual))
metab_dt$sample_type_origin <- factor(metab_dt$sample_type_origin)
metab_dt$graft_type <- factor(metab_dt$graft_type)
metab_split <- initial_split(metab_dt, strata = sample_type_origin, prop = 0.75)
metab_train <- training(metab_split) %>%
  select(-sample_type_origin)
metab_test <- testing(metab_split) %>%
  select(-sample_type_origin)

## model building
svm_mod <-
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_mode("classification") %>%
  set_engine("kernlab")

## create a workflow that encapsulates a recipe and a model
svm_wflow <- recipe(graft_type ~ ., data = metab_train) %>% 
  # remove any zero variance predictors
  step_zv(all_predictors()) %>%
  # remove any linear combinations
  step_lincomb(all_numeric()) %>%
  step_normalize(all_numeric_predictors()) %>% 
  workflow(svm_mod)
svm_wflow

## k fold cross validation.
metab_rs <- vfold_cv(metab_train, v = 5, repeats = 5)   

## area under the curve (auc)
roc_vals <- metric_set(roc_auc)

ctrl <- control_grid(verbose = FALSE, save_pred = TRUE)

## resampling (fine tuning, we apply these specifications we have just made for the fine tuning workflow 
## resampling dataset using cross validation)
formula_res <- tune_grid(
  object = svm_wflow,
  resamples = metab_rs,
  metrics = roc_vals,
  control = ctrl
)
formula_res

# show best submodel 
formula_res %>% 
  show_best(metric = "roc_auc")

# select best model hyperparameters
best_svm <- formula_res %>% 
  select_best(metric = "roc_auc")

best_svm

## finalize workflow
final_wflow <- svm_wflow %>% 
  finalize_workflow(best_svm)

final_wflow

# the last fit
final_fit <- last_fit(object = final_wflow, split = metab_split)

## collect metrics
final_fit %>% 
  collect_metrics()

## confusion matrix
cm <- final_fit %>% 
  collect_predictions() %>% 
  conf_mat(truth = graft_type, estimate = .pred_class)

print(cm)

autoplot(cm, type = "heatmap")