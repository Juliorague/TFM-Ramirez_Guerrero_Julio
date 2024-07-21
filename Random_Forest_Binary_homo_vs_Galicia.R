library(vip)
library(ggplot2)
library(tidyverse)
library(tidymodels)
library(data.table)
library(randomForest)
library(repr) 
library(dplyr)

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

## preprocessing
metab_recipe <- metab_train %>%
  recipe(graft_type ~ .) %>%
  step_corr(all_predictors(), threshold = 0.9) %>% 
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_normalize(all_numeric(), -all_outcomes()) %>%
  step_impute_knn(all_numeric(), neighbors = 5)

prep_metab <- prep(metab_recipe)
print(prep_metab)

training_set <- juice(prep_metab)
head(training_set)

## model building
tune_spec <- rand_forest(   # when the Random Forest model is specified using the rand_forest() function, 
  mtry = tune(),            # it internally uses bagging to build multiple decision trees
  trees = 100,
  min_n = tune()
) %>%
  set_mode("classification") %>%
  set_engine("randomForest")

tune_wf <- workflow() %>%
  add_formula(graft_type ~ .) %>%
  add_model(tune_spec)

## hyperparameter tuning

# we use k-fold cross-validation to tune the hyperparameters in the training set
trees_folds <- vfold_cv(training_set, v = 5, repeats = 5)
print(trees_folds)

# resampling (cross validation to fine tune the hyperparameters)
doParallel::registerDoParallel()

tune_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = 100 ## n. of tuning combinations 
)

print(tune_res)

options(repr.plot.width=14, repr.plot.height=8)

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  dplyr::select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

m <- round(sqrt(ncol(training_set)-1),0)
print(m)

rf_grid <- grid_regular(
  mtry(range = c(m-10, m+10)),
  min_n(range = c(7, 11)),
  levels = c(6,8) 
)
print(rf_grid)

regular_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = rf_grid
)
print(regular_res)

regular_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, linewidth = 1.5) +  
  geom_point() +
  labs(y = "AUC")

## final model 

# 1. selecting the best model based on AUC
best_auc <- select_best(regular_res, metric = "roc_auc") 

show_best(regular_res, metric = "roc_auc", n=1)

# compare it with what we get with the random grid (leaving random forest select and fine tune the best 
# hyperparameters mtry and min_n)
best_auc <- select_best(tune_res, metric = "roc_auc") 

show_best(tune_res, metric = "roc_auc", n=1)

# 2. finalise the model
final_rf <- finalize_model( # bagging occurs as part of the Random Forest training process
  tune_spec,
  best_auc
)
print(final_rf)

# 3. finalise the workflow and fit it to the initial split (training and test data)
final_wf <- workflow() %>%
  add_recipe(metab_recipe) %>%
  add_model(final_rf)

final_res <- final_wf %>%
  last_fit(metab_split) # bagging occurs as part of the Random Forest training process

# 4. evaluate the fine-tuned RF model in the test set. By evaluating the model we obtain the AUC and the 
# accuracy
print(final_res)

final_res %>%
  collect_metrics()

# 5. get variable importance
variable_importance <- final_res %>%
  pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>%   
  vip(num_features = 10)

variable_importance

importance_data <- variable_importance$data

## predictions
final_res %>%
  collect_predictions()

# confusion matrix
cm <- final_res %>%
  collect_predictions() %>%
  conf_mat(graft_type, .pred_class)

print(cm)

autoplot(cm, type = "heatmap")
