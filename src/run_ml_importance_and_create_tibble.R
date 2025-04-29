# 
#' Title Function to run ml and create tibble with the feature importance
#'
#' @param model_tb tibble with the response variable and the explanatory ones in wide format
#' @param asv_num asv name 
#' @param seed seed number
#'
#' @return
#' @export
#'
#' @examples
run_ml_importance_and_create_tibble <- function(model_tb, asv_num, seed) {
  # Run ml
  rf <- mikropml::run_ml(
    model_tb,
    method = 'rf',
    outcome_colname = 'diff_rclr_time',
    find_feature_importance = TRUE, 
    cross_val = control_cv,
    training_frac = 0.8,
    ntree = 100,
    perf_metric_name = 'RMSE',
    seed = seed
  )
  
  # Extract performance metrics
  rf_importance <- rf$feature_importance %>%
    as_tibble() %>%
    dplyr::mutate(asv_num = asv_num)
  
  return(rf_importance)
}
