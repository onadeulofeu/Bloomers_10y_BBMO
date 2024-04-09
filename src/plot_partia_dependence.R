# Function to plot partial dependence for each predictor variable
#' Title
#'
#' @param pdp_data_list 
#' @param importance_df 
#' @param asv_num_id 
#' @param num_plots 
#'
#' @return
#' @export
#'
#' @examples
#' 
plot_partial_dependence <- function(pdp_data_list, importance_df, asv_num_id, num_plots = 6) {
  # Initialize a list to store all partial plots
  all_plots <- list()
  
  # Sort predictor names based on their importance
  sorted_predictor_names <- importance_df |>
    dplyr::filter(asv_num == asv_num_id) |>
    slice_tail(n = num_plots) |>
    arrange(desc(abs(mean_diff))) |>
    dplyr::select(feat) |>
    as_vector()
  
  for (predictor_name in sorted_predictor_names) {
    # Combine data from all seeds for the current predictor
    plot_data <- data.frame()
    for (i in seq_along(pdp_data_list)) {
      seed_name <- paste("Seed_", i, sep = "")
      seed_data <- pdp_data_list[[i]][[predictor_name]]
      seed_data <- cbind(seed_data, Seed = seed_name)
      plot_data <- rbind(plot_data, seed_data)
    }
    
    ## extract axis x title
    axis_x_title <- paste0("'",predictor_name, "'", collapse = "")
    
    labs_env_models_no_nas_ed <- labs_env_models_no_nas |>
      as.data.frame() |>
      rownames_to_column() |>
      as_tibble() |>
      dplyr::filter(rowname == predictor_name)
    
    # Plot partial dependence for the current predictor
    partial_plot <- ggplot(plot_data, aes_string(x = predictor_name, y = "yhat")) +
      geom_line(aes(group = Seed)) +
      labs(x = labs_env_models_no_nas_ed$labs_env_models_no_nas, y = expression("Predicted rCLR change (t"^"-1"~")"),
           title = asv_num_id) +
      geom_hline(yintercept = 0, linetype = 'dashed') +
      theme_bw() +
      theme(panel.grid = element_blank(),
            text = element_text(size = 5))
    
    # Add the partial plot to the list
    all_plots[[predictor_name]] <- partial_plot
  }
  
  # Arrange all the partial plots in a single grid
  all_plots_grid <- do.call(gridExtra::grid.arrange, c(all_plots, ncol = 3))
  
  # Save the combined plot as a PDF
  ggsave(all_plots_grid, filename = paste0('pdp_',asv_num_id, '_', num_plots,'.pdf'),
         path = 'Results/Figures/random_forest/',
         width = 188, height = (num_plots/3*66), units = 'mm')
  
  # Save the all_plots_grid object in the R environment
  assign(paste0("pdp_all_plots_grid_", asv_num_id), all_plots_grid, envir = .GlobalEnv)
}