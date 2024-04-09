#' Create tibble to perform the randomforest
#'
#' @param data_diff a tibble with the differnces between one point and the next one in the timeseries
#' @param data_previous_ab a tibble with the abundances, to add the previous timepoint as a explantory variable
#' @param asv_num the asv number that we are analyzing
#' @param fraction the fraction at which the asv belong to
#'
#' @return this function will return a tibble ready to perform the random forest with the caret package
#' @export
#'
#' @examples
create_model_tb <- function(data_diff, data_previous_ab, asv_num, fraction = c('0.2', '3')) {
  # Filter the differences between the timepoint and the previous one
  # Check if fraction is valid
  if (!(fraction %in% c('0.2', '3'))) {
    stop("Invalid fraction value. Fraction must be either '0.2' or '3'.")
  }
  
  model_tb <- data_diff |>
    dplyr::filter(!is.na(diff_time)) |>
    dplyr::select(-sample_id_num) |>
    dplyr::filter(asv_num == {{asv_num}} & 
                    fraction == {{fraction}}) |>
    ungroup() |>
    dplyr::select(-asv_num, -diff_rclr, -diff_time, -date, -fraction, -fitted_rclr) |>
    bind_cols(env_data_interpolated_values_all_red) |>
    dplyr::select(-decimal_date)
  
  # # Filter abundances to add the previous one as a explanatory variable
  asvs_sm_asv <- data_previous_ab |>
    dplyr::filter(asv_num == {{asv_num}} &
                    fraction == {{fraction}}) |>
    arrange(decimal_date) |>
    dplyr::select(-asv_num, -fraction, -decimal_date)
  # 
  # Remove the last row from asvs_sm_asv
  asvs_sm_asv <- asvs_sm_asv[-nrow(asvs_sm_asv),] # the last rCLR is not needed since it won't have the next value to explain the change in abundance.
  # 
  # Combine data to create model_tb
  model_tb <- bind_cols(model_tb, asvs_sm_asv)  #|>
    #dplyr::select(-sample_id_num, -BP_FC1.55_no_nas)
  
  # Check if all columns are numeric
  if (!all(sapply(model_tb, is.numeric))) {
    stop("Not all columns in the model_tb are numeric.")
  }
  
  return(model_tb)
}
