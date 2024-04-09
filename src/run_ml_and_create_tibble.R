# Function to run ml and create tibble with the cv information for each seed and the R2 and the R2 of resampling (training)
#' Title 
#'
#' @param model_tb a tibble with the response variable and the explanatory variables in wide format
#' @param asv_num id of the ASV that we want to analyze
#' @param seed seed number
#'
#' @return 
#' @export
#'
#' @examples
#' run_ml_and_create_tibble(model_tb_11, 'asv11', seed_list[i])
#' 
#' 
run_ml_and_create_tibble <- function(model_tb, asv_num, seed) {
  # Run ml
  rf <- run_ml(
    model_tb,
    method = 'rf',
    outcome_colname = 'diff_rclr_time',
    find_feature_importance = FALSE, # it runs faster
    cross_val = control_cv,
    training_frac = 0.8,
    ntree = 1000,
    perf_metric_name = 'RMSE',
    seed = seed
  )
  
  # Extract performance metrics
  rf_performance <- rf$performance %>%
    as_tibble() %>%
    dplyr::mutate(asv_num = asv_num)
  
  # Extract resample results from final model
  resample_results <- rf$trained_model$resample %>%
    as_tibble()
  
  # Compute mean Rsquared from resample results
  mean_Rsquared <- mean(resample_results$Rsquared, na.rm = TRUE)
  
  # Add mean Rsquared to rf_performance
  rf_performance <- rf_performance %>%
    dplyr::mutate(mean_Rsquared_resample = mean_Rsquared)
  
  return(rf_performance)
}