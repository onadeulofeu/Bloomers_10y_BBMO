#' Count number of potential bloomers at different relative abunance () threshold
#'
#' @param threshold relative abundance at which you would like to set the threshold to consider a blooming event
#' @param fraction size fraction at which you are working on
#' @param asv_tab_pseudo table with relative abunances
#' @param asv_tab_zclr table with z scores for all the taxa at each point of the timeseries
#'
#' @return 
#' @export
#'
#' @examples
#' n_0.2_0 <- conunt_num_bloomers(threshold = 0, fraction = '0.2', asv_tab_pseudo =  asv_tab_10y_02_pseudo,
#' asv_tab_zclr = asv_tab_10y_02_zclr)
#' 

conunt_num_bloomers <- function(threshold, fraction, asv_tab_pseudo, asv_tab_zclr) {
  result <- asv_tab_pseudo |>
    inner_join(asv_tab_zclr, by = c('sample_id', 'asv_num')) |>
    group_by(asv_num) |>
    dplyr::filter(any(relative_abundance >= threshold)) |>
    dplyr::filter(zclr >= 1.96) |>
    dplyr::distinct(asv_num) |>
    reframe( n = n()) |>
    ungroup() |>
    reframe(num = sum(n)) |>
    dplyr::mutate(threshold = threshold, fraction = fraction)
  
  return(result)
}
