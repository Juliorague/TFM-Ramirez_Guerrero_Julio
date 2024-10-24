library(tidyverse)
library(tidymodels)
library(vip)

### Exploratory analysis ###

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

# transform the dataset from wide to long format to visualise if there are scale differences in the 
# metabolites
mm_metab <- homo_vs_Galicia %>%
  gather(key = "identifier", value = "expression", -c(individual, sample_type_origin, graft_type))

# make a table where it is represented the maximum expression of a metabolite versus the minimum to see 
#the difference that can occur between different metabolites
group_by(mm_metab, identifier) %>% 
  summarise(max_each = max(expression)) %>% 
  summarise(min = min(max_each), max = max(max_each))

# plot the distributions of values for all metabolites
options(repr.plot.width=14, repr.plot.height=8)

mm_metab %>%
  ggplot(aes(identifier, expression, fill = as.factor(identifier))) +
  geom_boxplot(show.legend = FALSE)



### Division of training and test sets ###

# keep only numeric features for convinience
metab_dt <- select(homo_vs_Galicia, -c(individual))
metab_dt$sample_type_origin <- factor(metab_dt$sample_type_origin)
metab_dt$graft_type <- factor(metab_dt$graft_type)
metab_split <- initial_split(metab_dt, strata = sample_type_origin, prop = 0.75)
metab_train <- training(metab_split) %>%
  select(-sample_type_origin)
metab_test <- testing(metab_split) %>%
  select(-sample_type_origin)



### Preprocessing ###

## build a recipe for preprocessing
metab_rec <- recipe(graft_type ~ ., data = metab_train) %>%
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_normalize(all_numeric(), -all_outcomes())

print(metab_rec)

# obtention of preprocessed training set
metab_prep <- metab_rec %>%
  prep(strings_as_factors = FALSE)

print(metab_prep)

metab_train <- juice(metab_prep)

options(repr.plot.width=14, repr.plot.height=8)

mm_train <- metab_train %>%                      
  gather(key = "identifier", value = "expression", -c(graft_type))

group_by(mm_train, identifier) %>% summarise(mean(expression),sd(expression))

mm_train %>%
  ggplot(aes(identifier, expression, fill = as.factor(identifier))) +
  geom_boxplot(show.legend = FALSE)



### Lasso Model ###

# add everything to a workflow object, piecewise, and fit the Lasso model
lasso_spec <- logistic_reg(mode = "classification", penalty = 0.1, mixture = 1) %>%
  set_engine("glmnet")

print(lasso_spec)

# make a workflow with the recipe established in the preprocessing that gives us the 
# normalised data, and the model we just built
wf <- workflow() %>%
  add_recipe(metab_rec) %>%
  add_model(lasso_spec)

print(wf)

# let's fit the workflow to the training data
lasso_fit <- wf %>%
  fit(data = metab_train)

# order the model parameters (coefficients)
lasso_fit %>%
  extract_fit_parsnip() %>%
  tidy()

# filtering of coefficients greater than 0
lasso_fit %>%
  extract_fit_parsnip() %>%
  tidy() %>%
  filter(estimate > 0)



### Tuning the hyperparameters ###

# preparation of the data for cross validation
metab_cv <- vfold_cv(metab_train, v=5, repeats = 10, strata = graft_type)

tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")

lambda_grid <- grid_regular(penalty(), levels = 50, filter = penalty <= .05)

print(lambda_grid)

# add the recipe done before and also the model with the tune specifications (model glmnet with 
# lasso)
wf1 <- workflow() %>%
  add_recipe(metab_rec) %>%
  add_model(tune_spec) # remember: the model equation was specified in the recipe (top of this document)

# fine-tune the model
doParallel::registerDoParallel()

lasso_grid <- tune_grid(
  wf1,
  resamples = metab_cv,
  grid = lambda_grid
)

# results for each value of the penalty parameter that was tried in the fine-tuning process
lasso_grid %>%
  collect_metrics()

# here are the penalties of the model (of each of the lambdas tested). Plotting the results will 
# help us see what happened during fine-tuning of the model, and then select the best value for (the 
# penalty parameter) based on the maximum AUC (binary classification problem)
lasso_grid %>%
  collect_metrics() %>%
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_errorbar(aes(
    ymin = mean - std_err,
    ymax = mean + std_err
  ),
  alpha = 0.5
  ) +
  geom_line(linewidth = 1.5) +                
  facet_wrap(~.metric, scales = "free", nrow = 2) +
  scale_x_log10() +
  theme(legend.position = "none")

# extract the best value of lambda (the one that confers the smallest roc) and specify it according to 
# the value of auc
lowest_roc <- lasso_grid %>%
  select_best(metric = "roc_auc")

print(lowest_roc)

# now let's fit the final model, the final workflow (we start from the workflow we built earlier wf1 
# and add the best value of roc that we have obtained)
final_lasso <- finalize_workflow(
  wf1,
  lowest_roc
)

print(final_lasso)



### Testing the model ###

# we are now ready to test our fine-tuned Lasso model on the test partition
lr_res <- last_fit(
  final_lasso,
  metab_split
)

lr_res %>%
  collect_metrics()

# confusion matrix
cm <- collect_predictions(lr_res) %>%
  conf_mat(graft_type, .pred_class)

print(cm)

autoplot(cm, type="heatmap")

# plot the roc curve in the testing data
lr_auc <- lr_res %>%
  collect_predictions() %>%
  roc_curve(graft_type, .pred_homo) %>%
  mutate(model = "Logistic Regression")

autoplot(lr_auc)



### Variable importance ###

importance_data <- final_lasso %>%
  fit(metab_train) %>%
  extract_fit_parsnip() %>%
  vi(lambda = lowest_roc$penalty)

importance_data %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  filter(Importance > 0) %>%
  top_n(10, Importance) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL)
