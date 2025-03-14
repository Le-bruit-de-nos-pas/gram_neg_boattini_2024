library(data.table)
library(tidyverse)
library(readxl)
library(purrr)
options(scipen = 999)

db_paulo <- read_excel(path = "data/db_paulo.xlsx")

db_paulo <- db_paulo %>% mutate(group=ifelse(Infections==1,1,0)) %>%
  select(-Infections,  -Carriage_only)

df <- db_paulo 

# Identify categorical (binary) and continuous variables
is_binary <- function(x) length(unique(x)) == 2

binary_vars <- df %>% select(where(is_binary)) %>% colnames()
continuous_vars <- df %>% select(where(~ !is_binary(.))) %>% colnames()

# Define group variable 
group_var <- "group"


# Wilcoxon test for continuous variables
wilcox_results <- map_dfr(continuous_vars, function(var) {
  test <- wilcox.test(df[[var]] ~ df[[group_var]])
  tibble(variable = var, p_value = test$p.value, test = "Wilcoxon")
})



# Fisher's Exact test for categorical variables
fisher_results <- map_dfr(binary_vars, function(var) {
  tbl <- table(df[[var]], df[[group_var]])
  if (all(dim(tbl) == c(2, 2))) { # Only apply if valid 2x2 table
    test <- fisher.test(tbl)
    tibble(variable = var, p_value = test$p.value, test = "Fisher")
  } else {
    tibble(variable = var, p_value = NA, test = "Not applicable")
  }
})


# Combine results
results <- bind_rows(wilcox_results, fisher_results) 

# Print results
print(results)



















library(dplyr)
library(pROC)
library(randomForest)
library(glmnet)
library(caret)

# Ensure 'group' is a factor
df[[group_var]] <- as.factor(df[[group_var]])


# Split into training and testing sets
set.seed(123)  
train_index <- createDataPartition(df[[group_var]], p = 0.7, list = FALSE)
train_data <- df[train_index, ]
test_data <- df[-train_index, ]


# **1. Logistic Regression with LASSO (Feature Selection)**
x_train <- model.matrix(group ~ ., data = train_data)[, -1]  # Remove intercept
y_train <- train_data[[group_var]]

lasso_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)

# Extract important variables
best_lambda <- lasso_model$lambda.min
lasso_coef <- coef(lasso_model, s = best_lambda)

important_vars <- rownames(lasso_coef)[which(lasso_coef != 0)]
print(important_vars)


# **2. Random Forest Feature Importance**
rf_model <- randomForest(group ~ ., data = train_data, importance = TRUE)
varImpPlot(rf_model)


# **3. Model Performance using ROC-AUC**
x_test <- model.matrix(group ~ ., data = test_data)[, -1]
y_test <- as.numeric(test_data[[group_var]])

pred_probs <- predict(lasso_model, newx = x_test, s = best_lambda, type = "response")
roc_curve <- roc(y_test, as.vector(pred_probs))
auc_value <- auc(roc_curve)


# Plot ROC Curve
plot(roc_curve, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))



















# Ensure 'group' is a factor
df[[group_var]] <- as.factor(df[[group_var]])

# Create an empty vector to store predictions
true_labels <- df[[group_var]]
predicted_probs <- rep(NA, nrow(df))

# **Leave-One-Out Cross-Validation (LOOCV)**
for (i in 1:nrow(df)) {
  train_data <- df[-i, ]
  test_data <- df[i, , drop = FALSE]  # Leave-one-out
  
  # Train the random forest model
  rf_model <- randomForest(group ~ ., data = train_data, importance = TRUE, ntree = 500)
  
  # Predict probability for the left-out sample
  predicted_probs[i] <- predict(rf_model, test_data, type = "prob")[, 2]
}

# **Evaluate Model Performance**
roc_curve <- roc(true_labels, predicted_probs)
auc_value <- auc(roc_curve)



# **Feature Importance**
rf_final <- randomForest(group ~ ., data = df, importance = TRUE, ntree = 500)
importance_df <- as.data.frame(importance(rf_final))
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]

# **Print Results**
print(importance_df)
print(paste("AUC-ROC:", round(auc_value, 3)))

# **Plot Feature Importance**
varImpPlot(rf_final)



# Select top 20 important features for better visualization
top_features_gini <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(20)

top_features_accuracy <- importance_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(20)

# **Plot 1: Mean Decrease in Gini**
ggplot(top_features_gini, aes(x = reorder(rownames(top_features_gini), MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#355C7D", alpha = 0.7) +  
  coord_flip() +  # Horizontal bars
  labs(title = "Feature Importance \n [Mean Decrease Gini]", x = "Features \n", y = "\n Importance") +
 theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = -0.5)) +
  theme(plot.title = element_text(hjust = 0.5))  # Centered title


# **Plot 2: Mean Decrease in Accuracy**
ggplot(top_features_accuracy, aes(x = reorder(rownames(top_features_accuracy), MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#6C5B7B", alpha = 0.7) +  
  coord_flip() +
  labs(title = "Feature Importance \n [Mean Decrease Accuracy]", x = "Features \n", y = "\n Importance") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = -0.5)) +
  theme(plot.title = element_text(hjust = 0.5))  # Centered title



# **Plot ROC Curve**
plot(roc_curve, 
     main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"), 
     col = adjustcolor("#C06C84", alpha.f = 0.8), # Firebrick color with 60% opacity
     lwd = 7)  # Line width for better visibility










 # IN HOSPITAL MORTALITY 
db_paulo <- read_excel(path = "data/db_paulo.xlsx")

df <- db_paulo 

# Identify categorical (binary) and continuous variables
is_binary <- function(x) length(unique(x)) == 2

binary_vars <- df %>% select(where(is_binary)) %>% colnames()
continuous_vars <- df %>% select(where(~ !is_binary(.))) %>% colnames()

# Define group variable 
group_var <- "In_hospital_mortality"


# Wilcoxon test for continuous variables
wilcox_results <- map_dfr(continuous_vars, function(var) {
  test <- wilcox.test(df[[var]] ~ df[[group_var]])
  tibble(variable = var, p_value = test$p.value, test = "Wilcoxon")
})



# Fisher's Exact test for categorical variables
fisher_results <- map_dfr(binary_vars, function(var) {
  tbl <- table(df[[var]], df[[group_var]])
  if (all(dim(tbl) == c(2, 2))) { # Only apply if valid 2x2 table
    test <- fisher.test(tbl)
    tibble(variable = var, p_value = test$p.value, test = "Fisher")
  } else {
    tibble(variable = var, p_value = NA, test = "Not applicable")
  }
})


# Combine results
results <- bind_rows(wilcox_results, fisher_results) 

# Print results
print(results)

data.frame(results)















library(dplyr)
library(pROC)
library(randomForest)
library(glmnet)
library(caret)


# Ensure 'group' is a factor
df[[group_var]] <- as.factor(df[[group_var]])

names(df)

df <- df %>% select(-Mortality_14d, -Mortality_30d)

# Split into training and testing sets
set.seed(123)  
train_index <- createDataPartition(df[[group_var]], p = 0.7, list = FALSE)
train_data <- df[train_index, ]
test_data <- df[-train_index, ]


# **1. Logistic Regression with LASSO (Feature Selection)**
x_train <- model.matrix(In_hospital_mortality ~ ., data = train_data)[, -1]  # Remove intercept
y_train <- train_data[[group_var]]

lasso_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)

# Extract important variables
best_lambda <- lasso_model$lambda.min
lasso_coef <- coef(lasso_model, s = best_lambda)

important_vars <- rownames(lasso_coef)[which(lasso_coef != 0)]
print(important_vars)


# **2. Random Forest Feature Importance**
rf_model <- randomForest(In_hospital_mortality ~ ., data = train_data, importance = TRUE)
varImpPlot(rf_model)


# **3. Model Performance using ROC-AUC**
x_test <- model.matrix(In_hospital_mortality ~ ., data = test_data)[, -1]
y_test <- as.numeric(test_data[[group_var]])

pred_probs <- predict(lasso_model, newx = x_test, s = best_lambda, type = "response")
roc_curve <- roc(y_test, as.vector(pred_probs))
auc_value <- auc(roc_curve)


# Plot ROC Curve
plot(roc_curve, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))
















# Ensure 'group' is a factor
df[[group_var]] <- as.factor(df[[group_var]])

# Create an empty vector to store predictions
true_labels <- df[[group_var]]
predicted_probs <- rep(NA, nrow(df))

# **Leave-One-Out Cross-Validation (LOOCV)**
for (i in 1:nrow(df)) {
  train_data <- df[-i, ]
  test_data <- df[i, , drop = FALSE]  # Leave-one-out
  
  # Train the random forest model
  rf_model <- randomForest(In_hospital_mortality ~ ., data = train_data, importance = TRUE, ntree = 500)
  
  # Predict probability for the left-out sample
  predicted_probs[i] <- predict(rf_model, test_data, type = "prob")[, 2]
}

# **Evaluate Model Performance**
roc_curve <- roc(true_labels, predicted_probs)
auc_value <- auc(roc_curve)



# **Feature Importance**
rf_final <- randomForest(In_hospital_mortality ~ ., data = df, importance = TRUE, ntree = 500)
importance_df <- as.data.frame(importance(rf_final))
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]

# **Print Results**
print(importance_df)
print(paste("AUC-ROC:", round(auc_value, 3)))

# **Plot Feature Importance**
varImpPlot(rf_final)



# Select top 20 important features for better visualization
top_features_gini <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(20)

top_features_accuracy <- importance_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(20)

# **Plot 1: Mean Decrease in Gini**
ggplot(top_features_gini, aes(x = reorder(rownames(top_features_gini), MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#355C7D", alpha = 0.7) +  
  coord_flip() +  # Horizontal bars
  labs(title = "Feature Importance \n [Mean Decrease Gini]", x = "Features \n", y = "\n Importance") +
 theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = -0.5)) +
  theme(plot.title = element_text(hjust = 0.5))  # Centered title


# **Plot 2: Mean Decrease in Accuracy**
ggplot(top_features_accuracy, aes(x = reorder(rownames(top_features_accuracy), MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#6C5B7B", alpha = 0.7) +  
  coord_flip() +
  labs(title = "Feature Importance \n [Mean Decrease Accuracy]", x = "Features \n", y = "\n Importance") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = -0.5)) +
  theme(plot.title = element_text(hjust = 0.5))  # Centered title



# **Plot ROC Curve**
plot(roc_curve, 
     main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"), 
     col = adjustcolor("#C06C84", alpha.f = 0.8), # Firebrick color with 60% opacity
     lwd = 7)  # Line width for better visibility




