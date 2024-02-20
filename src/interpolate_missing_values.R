#' Interpolate missing values
#'
#' @param data tibble in a long format with a column containing the variables name, values and the date in decimal format
#' @param variable_to_inter name of the variable to interpolate the missing values
#' @param span_value value used in the loess function which controls the degree of smoothing.
#'
#' @return
#' @export
#'
#' @examples
#' interpolated_data <- interpolate_missing(bbmo_env_long, variable_to_inter = "day_length", span_value = 0.1)
#' 
#' 
interpolate_missing <- function(data, variable_to_inter, span_value) {
  
  env_to_fit <- data |>
    dplyr::filter(variable == {{variable_to_inter}}) |>
    dplyr::select(value, decimal_date) |>
    dplyr::filter(!is.na(value))
  
  env_fitted <- data %>%
    dplyr::filter(variable == {{variable_to_inter}}) |>
    dplyr::select(value, decimal_date)
  
  fit <- loess(value ~ decimal_date, data = env_to_fit, span = span_value)
  
  env_fitted[[paste0(variable_to_inter, "_interpolated")]] <- predict(fit, newdata = env_fitted)
  
  env_fitted <- env_fitted |>
    rename(!!variable_to_inter := value)
  
  return(env_fitted)
}