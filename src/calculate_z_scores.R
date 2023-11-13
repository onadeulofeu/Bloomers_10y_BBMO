#' Calculate z scores
#' 
#' This function is used to calculate z scores of variables, with the option to work in a long format tibble and be able to normalize different groups
#' of variables.
#'
#' @param data a tibble with a column containing the values we want to normalize, without NAs
#' @param col the column of the tibble containing the values 
#' @param name name of the output column 
#' @param group in case we have the dataset in a long format the environmental groups we want to calculate the mean & sd for normalizing by z scores.
#'
#' @return
#' @export
#'
#' @examples
calculate_z_score <- function(data, col, name = NULL, group = NULL) {
  stopifnot(is.numeric(data[[col]]))
  
  # Check for NAs in the specified column
  if (anyNA(data[[col]])) {
    warning("The specified column contains NA values.")
  }
  
  col_name <- ifelse(!is.null(name), paste0("z_score_", name), 'z_score')
  
  if (!is.null(group)) {
    data <- data |>
      dplyr::group_by(!!sym(group)) |>
      dplyr::mutate(!!col_name := (!!sym(col) - base::mean(!!sym(col), na.rm = TRUE)) / stats::sd(!!sym(col), na.rm = TRUE))
  } else {
    data <- data |>
      dplyr::mutate(!!col_name := (!!sym(col) - base::mean(!!sym(col), na.rm = TRUE)) / stats::sd(!!sym(col), na.rm = TRUE))
  }
  
  return(data)
}