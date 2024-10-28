# Supervised Learning Metabolite Analysis of Drought Tolerance and Grafting Effects in *Pinus pinaster* Aiton

### Author: Julio Ramírez Guerrero  

## Overview
This repository contains scripts used for a study focused on understanding the impact of grafting on the metabolomic response of *Pinus pinaster* Aiton with a specific emphasis on drought tolerance. By using supervised machine learning (ML) algorithms, the study identifies key metabolites and molecular markers of grafting, rootstock type, and scion type to infer long-distance metabolic communication within grafted conifers.

## Project Scope
### Objectives
In this study, we aim to explore drought tolerance and long-distance communication between *P. pinaster* organs. To achieve this:
1. **Analyse and integrate the metabolomic profiles from needles, stems, and roots** of both grafted and non-grafted *P. pinaster* plants using supervised learning models.
2. **Assess the influence of grafting, rootstock genotype, scion genotype, and their interactions** on the plant's metabolomic response.

### Approach
Using supervised learning models, this study evaluates the metabolite profiles across different graft configurations (homo-grafted, hetero-grafted, and non-grafted). In the project several classification algorithms have been adapted and tested for the specific models to see wich one is the best for each model:
- Lasso-penalized Logistic Regression
- Random Forest (RF)
- Boosting
- Support Vector Machines (SVM)

## Repository Structure
- **`Boosting_models_scripts`**: Includes 4 R scripts of the models evaluated (Binary_Rootstock_type, Binary_homo_vs_Galicia, Binary_homo_vs_absent and Multiclass_Rootstock-Graft) using Boositng supervised learning method.
- **`Lasso-penalized_Logistic_Regression_models_scripts`**: Includes 4 R scripts of the models evaluated (Binary_Rootstock_type, Binary_homo_vs_Galicia, Binary_homo_vs_absent and Multiclass_Rootstock-Graft) using Lasso-penalized Logistic Regressio supervised learning method.
- **`Random_Forest_models_scripts`**: Includes 4 R scripts of the models evaluated (Binary_Rootstock_type, Binary_homo_vs_Galicia, Binary_homo_vs_absent and Multiclass_Rootstock-Graft) using Random Forest supervised learning method.
- **`Support_Vector_Machine_models_scripts`**: Includes 4 R scripts of the models evaluated (Binary_Rootstock_type, Binary_homo_vs_Galicia, Binary_homo_vs_absent and Multiclass_Rootstock-Graft) using Support Vector Machine supervised learning method.
- **`Variable_Importance_arithmetic_mean`**: Contains a CSV file with the top 5 most important variables, with their importance, in each of the 10 iterations performed in all models (except SVM), and a R script for getting the top 5 variables among the iterations based on the arithmethic mean of the variable importance.

## Requirements
- **R v4.4.0** or higher
- R packages: `tidymodels v1.2.0`, `vip v0.4.1`, `tidyverse v2.0.0`

## Data Availability
The Metabolite_data.csv file required for the supervised learning models scripts is available upon request from the corresponding author, as it is part of an ongoing research project.
