#' Function to calculate residuals from mean abundance and plot them
#' 
#' What we intent to do here is to observe if some blooming events present a relative abundance lower than the 
#' mean abundance of that specific ASV in the whole dataset and additionally observe if the relative abundance of
#' these events is over the variance of their relative abundance or not.
#'
#' @param data_anom_abund tibble with the information of the relative abundance during a potential blooming event
#' @param data_mean_abund_asv tibble with the mean abundance of each ASV for the whole dataset
#' @param asv_num ASV number from the ASV we want to plot and observe its residuals.
#' @param community_fraction if we want to analyze the free living or particle attached community
#'
#' @return
#' @export
#'
#' @examples
#' Example of ASV17 individual plot with residuals:
# plot_residuals(data_anom_abund = anom_rel_abund,
#                data_mean_abund_asv = mean_abund_f_asv,
#                asv_num = 'asv17',
#                community_fraction = '0.2')

plot_residuals <- function(data_anom_abund, data_mean_abund_asv, asv_num, community_fraction){
  
  community_fraction <- match.arg(community_fraction, c('0.2', '3')) # in case the fraction don't match this values an error will appear
  
  data_mean_abund_f <- data_mean_abund_asv |>
    dplyr::filter(asv_num == {{asv_num}} &
                    fraction == {{community_fraction}})
  
  title_plot <- data_anom_abund |>
    dplyr::filter(asv_num == {{asv_num}} &
                    fraction == {{community_fraction}}) |>
    dplyr::mutate(asv_num_fa = paste0(asv_num, ' ', family)) |>
    distinct(asv_num_fa)
    
  plot_residuals <- data_anom_abund |>
    left_join(data_mean_abund_asv, by = c('asv_num', 'fraction')) |>
    dplyr::mutate(residual = abundance_value - mean_abund) |>
    dplyr::filter(asv_num == {{asv_num}} & 
                    fraction == {{community_fraction}}) |>
    ggplot(aes(date, residual))+
    geom_rect(data = timeseries_limits, mapping = aes(xmin = date_min, xmax = date_max,  x=NULL, y=NULL,
                                                      ymin = (data_mean_abund_f$mean_abund - data_mean_abund_f$sd_abund), 
                                                      ymax = (data_mean_abund_f$mean_abund + data_mean_abund_f$sd_abund)),
              fill = '#94969E',
              alpha = 0.4)+
    geom_hline(yintercept = data_mean_abund_f$mean_abund)+
    geom_point(aes(color = bloom_event))+
    scale_color_manual(values = c('bloom' = 'darkred', 'no-bloom' = 'white'))+
    geom_segment(aes(x = date, y = residual, xend = date, yend = data_mean_abund_f$mean_abund), linetype = "dashed")+ 
    labs(title = title_plot)+
    scale_x_datetime(date_labels = "%Y", expand = c(0,0), limits = c(min(timeseries_limits$date_min), max(timeseries_limits$date_max)))+
    labs(x = 'Time', y = 'Relative abundance (%)')+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          text = element_text(size = 7),
          axis.title.y = element_text(size = 5),
          legend.key.size = unit(4, 'mm'))
  
  return(plot_residuals)
  
}
