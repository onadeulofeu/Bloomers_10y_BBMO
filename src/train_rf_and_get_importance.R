# Function to train Random Forest model and calculate variable importance
#' Title
#'
#' @param model_tb_list 
#' @param outcome_colname column name of the response variable
#' @param seed_list list of seeds that we want to evaluate
#'
#' @return
#' @export
#'
#' @examples
train_rf_and_get_importance <- function(model_tb_list, outcome_colname, seed_list) {
  # Initialize an empty list to store results
  importance_list <- list()
  
  # Iterate over each ASV data tibble
  for (i in seq_along(model_tb_list)) {
    asv_data <- model_tb_list[[i]]
    
    # Extract predictor variables
    predictors <- dplyr::select(asv_data, -c(diff_rclr_time))
    
    # Initialize an empty list to store importance for different seeds
    seed_importance_list <- list()
    
    # Iterate over each seed
    for (seed in seed_list) {
      # Train Random Forest model
      trained_model <- train(
        x = predictors,
        y = asv_data[[outcome_colname]],
        method = 'rf',
        importance = TRUE,
        trainControl = trainControl(method = "cv", number = 10, seeds = seed)
      )
      
      # Extract ASV number from the name of the tibble
      asv_num <- names(model_tb_list)[i]
      
      # Create a new column for ASV number in the importance data frame
      importance <- data.frame(asv_num =  asv_num, importance = varImp(trained_model)$importance,
                               perc_IncMSE = trained_model$finalModel$importance[,1],
                               importanceSD = trained_model$finalModel$importanceSD) |>
        rownames_as_column(var = 'env_variable')
      
      # Store importance for the current seed
      seed_importance_list[[seed]] <- importance
    }
    
    # Combine importance for different seeds into a single tibble
    seed_importance_df <- bind_rows(seed_importance_list, .id = "seed")
    
    # Store importance for the current ASV
    importance_list[[i]] <- seed_importance_df
  }
  
  # Combine results for all ASVs into a single tibble
  importance_df <- bind_rows(importance_list)
  
  return(importance_df)
}