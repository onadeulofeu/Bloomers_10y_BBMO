#' Count number of potential bloomers at different relative abunance () threshold
#'
#' @param threshold relative abundance at which you would like to set the threshold to consider a blooming event
#' @param fraction size fraction at which you are working on
#' @param asv_tab_pseudo table with relative abundances
#' @param z_score_tb table with the z_scores for each asv and sample (columns needed asv_num and sample_id)
#'
#' @return 
#' @export
#'
#' @examples
#' n_0.2_0 <- conunt_num_bloomers(threshold = 0, fraction = '0.2', asv_tab_pseudo =  asv_tab_10y_02_pseudo,
#' asv_tab_zclr = asv_tab_10y_02_zclr)
#' 

count_num_bloomers <- function(threshold, fraction, asv_tab_pseudo, z_scores_tb) {
 
  result <- asv_tab_pseudo |>
    dplyr::left_join(z_scores_tb, by= c('sample_id', 'asv_num')) |>
    dplyr::filter(relative_abundance >= threshold &
                    z_score_ra >= 1.96) |>
    dplyr::distinct(asv_num) |>
    reframe( n = n()) |>
    ungroup() |>
    reframe(num = sum(n)) |>
    dplyr::mutate(threshold = threshold, fraction = fraction)
  
  return(result)
}
