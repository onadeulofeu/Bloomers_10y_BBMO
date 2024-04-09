### 
#' Title Function to train model and compute the partial dependence for each ASVs at different seeds. 
#'
#' @param data tibble with the response variables and explanatory variables in wide format
#' @param outcome_colname name of the response column
#' @param predictors names of the explanatory variables
#' @param seed_list number of seeds that we want to explore
#' @param grid_resolution number of data points for the partial analysis form the pdp package: Integer giving the number of equally spaced points to use for the continuous variables listed in pred.var when pred.grid is not supplied. If left NULL, it will default to the minimum between 51 and the number of unique data points for each of the continuous independent variables listed in pred.var.
#' @param asv_num the asv name that we are analyzing
#'
#' @return
#' @export
#'
#' @examples
train_rf_and_save_pdp <- function(data, outcome_colname, predictors, seed_list = 1:2,  grid_resolution = 100, asv_num) {
  # Initialize a list to store partial dependence data
  pdp_data_list <- list()
  
  # Loop through each seed
  for (seed in seed_list) {
    # Set seed for reproducibility
    set.seed(seed)
    
    # Train Random Forest model
    trained_model <- train(as.formula(paste(outcome_colname, "~ .")), data = data, method = "rf")
    
    # Initialize a list to store partial dependence data for this seed
    seed_pdp_data_list <- list()
    
    # Loop through each predictor variable
    for (var in predictors) {
      # Compute partial dependence for the current variable
      partial_dep <- pdp::partial(trained_model, pred.var = var, train = data, grid.resolution = grid_resolution)
      
      # Extract the predictor variable values
      pdp_data <- partial_dep |>
        as_tibble()
      
      # Store partial dependence data for this variable in the list
      seed_pdp_data_list[[var]] <- pdp_data
    }
    
    # Store partial dependence data for this seed in the main list
    pdp_data_list[[paste("seed_", seed, sep = "")]] <- seed_pdp_data_list
  }
  
  # Save partial dependence data
  saveRDS(pdp_data_list, file = paste0('results/figures/random_forest/pdp/pdp_', asv_num,"_data.RDS"))
}
