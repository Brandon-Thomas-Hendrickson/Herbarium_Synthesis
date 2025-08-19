library(car)
library(dplyr)
library(ggplot2)
library(ggforce)
library(randomForest)

set.seed(123) # For reproducibility

n_iter <- 1000
all_summaries <- list()
all_models <- list()
all_errors <- numeric(n_iter)

for (i in 1:n_iter) {
    # --- Sampling ---
    herb_EUR_list$row_id <- seq_len(nrow(herb_EUR_list))
    usa_rows <- herb_EUR_list %>% filter(Country == "USA")
    non_usa_rows <- herb_EUR_list %>% filter(Country != "USA")
    non_usa_sampled <- non_usa_rows %>%
        group_by(Country) %>%
        sample_n(size = min(6, n()), replace = FALSE) %>%
        ungroup()
    final_rows <- bind_rows(usa_rows, non_usa_sampled) %>%
        arrange(row_id)
    
    # --- Subset covariance matrix ---
    cov_matrix_downsampled <- cov_matrix[final_rows$row_id, final_rows$row_id]
    
    # --- PCA ---
    eigen_decomp_downsampled <- eigen(cov_matrix_downsampled)
    pc_axes <- eigen_decomp_downsampled$vectors[, 1:10]
    pca_data <- data.frame(
        PC1 = pc_axes[, 1],
        PC2 = pc_axes[, 2],
        PC3 = pc_axes[, 3],
        PC4 = pc_axes[, 4],
        PC5 = pc_axes[, 5],
        PC6 = pc_axes[, 6],
        PC7 = pc_axes[, 7],
        PC8 = pc_axes[, 8],
        PC9 = pc_axes[, 9],
        PC10 = pc_axes[, 10],
        Country = as.factor(final_rows$Country),
        RANGE = final_rows$RANGE
    )
    
    # --- Split data ---
    input.full.data <- pca_data
    training.set.Brandon <- input.full.data[input.full.data$RANGE == "EUR", ]
    test.set.Brandon <- input.full.data[!(input.full.data$RANGE == "EUR"), ]
    training.set.Brandon$Country <- droplevels(as.factor(training.set.Brandon$Country))
    
    # --- Random forest ---
    Model.10PCs <- randomForest(
        formula = Country ~ PC1 + PC2,
        data = training.set.Brandon,
        ntree = 5000
    )
    
    # --- Prediction and summary ---
    Predictions <- predict(Model.10PCs, test.set.Brandon[, 1:2])
    test.set.Brandon$prediction <- Predictions
    test.set.Brandon$Country <- as.character(test.set.Brandon$Country)
    test.set.Brandon$Count <- 1
    summary_df <- summarySE(test.set.Brandon, measurevar = "Count", groupvars = "prediction")
    summary_df$iteration <- i
    all_summaries[[i]] <- summary_df

    # --- Store model and error rate ---
    all_models[[i]] <- Model.10PCs
    # OOB error rate
    all_errors[i] <- Model.10PCs$err.rate[Model.10PCs$ntree, "OOB"]
}

# Combine all iterations
all_summaries_df <- bind_rows(all_summaries)

# Calculate mean and sd of N (Count) for each prediction group
final_summary <- all_summaries_df %>%
    group_by(prediction) %>%
    summarise(
        mean_N = mean(N, na.rm = TRUE),
        sd_N = sd(N, na.rm = TRUE),
        .groups = "drop"
    )

# ...existing code...

# Find the best model (lowest error rate)
best_iter <- which.min(all_errors)
best_model <- all_models[[best_iter]]
cat("Best model is from iteration:", best_iter, "\n")
print(best_model)

# Get the summarySE result for the best model
best_summarySE <- all_summaries[[best_iter]]
cat("summarySE for best model:\n")
print(best_summarySE)

print(final_summary)


##########################
#All PCs##################
##########################

set.seed(123) # For reproducibility

n_iter <- 1001
all_summaries <- list()
all_models <- list()
all_errors <- numeric(n_iter)

for (i in 1:n_iter) {
    # --- Sampling ---
    herb_EUR_list$row_id <- seq_len(nrow(herb_EUR_list))
    usa_rows <- herb_EUR_list %>% filter(Country == "USA")
    non_usa_rows <- herb_EUR_list %>% filter(Country != "USA")
    non_usa_sampled <- non_usa_rows %>%
        group_by(Country) %>%
        sample_n(size = min(6, n()), replace = FALSE) %>%
        ungroup()
    final_rows <- bind_rows(usa_rows, non_usa_sampled) %>%
        arrange(row_id)
    
    # --- Subset covariance matrix ---
    cov_matrix_downsampled <- cov_matrix[final_rows$row_id, final_rows$row_id]
    
    # --- PCA ---
    eigen_decomp_downsampled <- eigen(cov_matrix_downsampled)
    pc_axes <- eigen_decomp_downsampled$vectors[, 1:10]
    pca_data <- data.frame(
        PC1 = pc_axes[, 1],
        PC2 = pc_axes[, 2],
        PC3 = pc_axes[, 3],
        PC4 = pc_axes[, 4],
        PC5 = pc_axes[, 5],
        PC6 = pc_axes[, 6],
        PC7 = pc_axes[, 7],
        PC8 = pc_axes[, 8],
        PC9 = pc_axes[, 9],
        PC10 = pc_axes[, 10],
        Country = as.factor(final_rows$Country),
        RANGE = final_rows$RANGE
    )
    
    # --- Split data ---
    input.full.data <- pca_data
    training.set.Brandon <- input.full.data[input.full.data$RANGE == "EUR", ]
    test.set.Brandon <- input.full.data[!(input.full.data$RANGE == "EUR"), ]
    training.set.Brandon$Country <- droplevels(as.factor(training.set.Brandon$Country))
    
    # --- Random forest ---
    Model.10PCs <- randomForest(
        formula = Country ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
        data = training.set.Brandon,
        ntree = 5000
    )
    
    # --- Prediction and summary ---
    Predictions <- predict(Model.10PCs, test.set.Brandon[, 1:10])
    test.set.Brandon$prediction <- Predictions
    test.set.Brandon$Country <- as.character(test.set.Brandon$Country)
    test.set.Brandon$Count <- 1
    summary_df <- summarySE(test.set.Brandon, measurevar = "Count", groupvars = "prediction")
    summary_df$iteration <- i
    all_summaries[[i]] <- summary_df

    # --- Store model and error rate ---
    all_models[[i]] <- Model.10PCs
    # OOB error rate
    all_errors[i] <- Model.10PCs$err.rate[Model.10PCs$ntree, "OOB"]
}

# Combine all iterations
all_summaries_df <- bind_rows(all_summaries)

# Calculate mean and sd of N (Count) for each prediction group
final_summary <- all_summaries_df %>%
    group_by(prediction) %>%
    summarise(
        mean_N = mean(N, na.rm = TRUE),
        sd_N = sd(N, na.rm = TRUE),
        .groups = "drop"
    )

# ...existing code...

# Find the best model (lowest error rate)
best_iter <- which.min(all_errors)
best_model <- all_models[[best_iter]]
cat("Best model is from iteration:", best_iter, "\n")
print(best_model)

# Get the summarySE result for the best model
best_summarySE <- all_summaries[[best_iter]]
cat("summarySE for best model:\n")
print(best_summarySE)

print(final_summary)