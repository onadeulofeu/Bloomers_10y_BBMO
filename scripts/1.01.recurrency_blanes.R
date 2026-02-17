# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                     data analysis pipeline                  ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                    BBMO timeseries 10-Y data                ++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++                         metabarcoding                       ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++             Code developed by Ona Deulofeu-Capo 2024        ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Context for the BBMO time series 
### Edit the function so that we extract all comparisons between samples -----
dissimilarity_matrix_lag <- function(data, sample_id_col, values_cols_prefix) {
  
  # Extract rownames to mantain them at the output table
  sample_id_unique <- data %>%
    group_by({{sample_id_col}}) %>% ##sample id that identifies uniquely each sample
    #dplyr::arrange(.by_group = TRUE)  %>% ## reorder so that it is ordered equally
    dplyr::distinct({{sample_id_col}})
  
  # Index samples to filter for only consecutive comparisons
  # check name of the columns
  #if ({{sample_id_col}} %in% colnames(data)==FALSE) {stop("There is no sample_id_col column in your data tibble")}
  
  samples_index <- data %>%
    dplyr::group_by({{sample_id_col}}) %>%
    #dplyr::arrange(.by_group = TRUE) %>%
    dplyr::select({{sample_id_col}}) %>%
    dplyr::distinct({{sample_id_col}}) %>%
    as_tibble() %>%
    dplyr::mutate('row_index_2' := dplyr::row_number()) %>%
    dplyr::select({{sample_id_col}}, row_index_2)
  
  # Compute pairwise Bray-Curtis distances between all rows of data
  
  dissim_mat <- data %>%
    pivot_wider(id_cols = {{sample_id_col}},  names_from = asv_num, values_from = relative_abundance) %>%
    group_by({{sample_id_col}}) %>%
    #dplyr::arrange(.by_group = TRUE) %>%
    tibble::column_to_rownames('sample_id') %>%
    vegan::vegdist(method = 'bray', upper = T) %>%
    as.matrix() %>%
    as_data_frame() %>%
    cbind(sample_id_unique) %>%
    dplyr::mutate(row_index := row_number()) %>%
    rowid_to_column() %>%
    as_tibble() %>%
    pivot_longer(cols = starts_with({{values_cols_prefix}}), values_to = 'bray_curtis_result', names_to = 'samples') %>%
    left_join(samples_index, by = c('samples' = 'sample_id'))
  
  # # Chech that diagonal elements to zero (i.e., each sample is identical to itself)
  # diag(dissim_mat) <- 0
  
  # Return the dissimilarity matrix
  return(dissim_mat)
}

# Time lag recurrency Bray Curtis ----
## ASV tab rarefied (with EcoUtils package)
asv_10y_l_rel_rar <- asv_tab_bbmo_10y_w_rar |>
  pivot_longer(cols = starts_with('asv'), values_to = 'reads', names_to = 'asv_num') |>
  dplyr::mutate(reads = as.numeric(reads)) |>
  calculate_rel_abund(group_cols = sample_id) 

asv_tab_10y_3_rel_rar <- asv_tab_bbmo_10y_l_rel_rar |>
  dplyr::filter(sample_id %in% m_3$sample_id)

asv_tab_10y_filt_02_rel_rar <- asv_tab_bbmo_10y_l_rel_rar %>%
  dplyr::filter(sample_id %in% m_02$sample_id)

bray_curtis_02_rar <- dissimilarity_matrix_lag(data = asv_tab_10y_filt_02_rel_rar, 
                                           sample_id_col = sample_id,
                                           values_cols_prefix = 'BL')

bray_curtis_3_rar <- dissimilarity_matrix_lag(data = asv_tab_10y_3_rel_rar,
                                          sample_id_col = sample_id,
                                          values_cols_prefix = 'BL')

bray_curtis_02_lag <- bray_curtis_02_rar |>
  dplyr::mutate(lag = row_index_2 - row_index) |>
  dplyr::mutate(fraction = '0.2') 

m_02_ed <- m_02 |>
  dplyr::select(sample_id,  date, sample_id_num) |>
  separate(sample_id, sep = '_', into = c('sample_code', 'fraction', 'seqid'), remove = F) |>
  dplyr::select(-fraction, -seqid, -date, -sample_id)

bray_curtis_3_lag <- bray_curtis_3_rar |> ## distances won't be correctly defined since we have 3 missing samples in PA. I need to change number of samples correctly
  separate(sample_id, sep = '_', into = c('sample_code', 'fraction', 'seqid'), remove = F) |>
  left_join(m_02_ed, by = c('sample_code')) |>
  separate(samples, sep = '_', into = c('samplesid', 'fraction', 'seqid'), remove = F) |>
  left_join(m_02_ed, by = c('samplesid' = 'sample_code')) |>
  dplyr::mutate(lag = as.numeric(sample_id_num.x) - as.numeric(sample_id_num.y)) |>
  dplyr::mutate(fraction = '3') |>
  dplyr::select(-row_index, -row_index_2) |>
  dplyr::select(-sample_code, -samplesid, -seqid, row_index = sample_id_num.x, row_index_2 = sample_id_num.y)

## 
bray_curtis_lag_all <- bray_curtis_02_lag |>
  dplyr::mutate(row_index = as.character(row_index),
                row_index_2 = as.character(row_index_2)) |>
  bind_rows(bray_curtis_3_lag) |>
  dplyr::filter(lag > 0)

bray_curtis_lag_plot <- bray_curtis_lag_all |>
  ggplot(aes(lag, bray_curtis_result))+
  scale_x_continuous(breaks = c(12, 24, 36, 48, 12*5, 12*6, 12*7, 12*8, 12*9, 120), labels = c(1, 2, 3,4, 5, 6, 7, 8,  9, 10))+
  geom_point(aes(shape = fraction), alpha = 0.02)+
  scale_shape_discrete(labels = labs_fraction_env)+
  labs(y = 'Bray Curtis Dissimilarity', x = 'Number of years between samples', shape = '', color = '', fill = '')+
  geom_line(data = bray_curtis_lag_all |>
              group_by(fraction, lag) |>
              dplyr::reframe(mean = mean(bray_curtis_result)), aes(lag, mean, group = fraction, color = fraction), linewidth = 1)+
  scale_color_manual(values = palette_fraction_env, labels = labs_fraction_env)+
  scale_fill_manual(values = palette_fraction_env, labels = labs_fraction_env)+
  guides(shape = guide_legend(ncol = 2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom')

bray_curtis_lag_plot 

# ggsave(filename = 'bray_curtis_lag_plot.pdf',
#        plot = bray_curtis_lag_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 120, units = 'mm')

## UNIFRAC recurrence -----
m_02_ed <- m_02 |>
  dplyr::select(sample_id,  date, sample_id_num) |>
  separate(sample_id, sep = '_', into = c('sample_code', 'fraction', 'seqid'), remove = F) |>
  dplyr::select(-fraction, -seqid, -date, -sample_id)

unifrac_tibble_02 <- unifrac_tibble |>
  dplyr::filter(str_detect(sample_id_1, '_0.2_') &
                  str_detect(sample_id_2, '_0.2_')) |>
  separate(sample_id_1, sep = '_', into = c('sample_code_1', 'fraction', 'seqid'), remove = F) |>
  separate(sample_id_2, sep = '_', into = c('sample_code_2', 'fraction', 'seqid'), remove = F) |>
  left_join(m_02_ed, by = c('sample_code_1' = 'sample_code')) |>
  dplyr::select(sample_id_num_1 = sample_id_num, sample_code_1, sample_code_2, fraction, wunifrac_distance ) |>
  left_join(m_02_ed, by = c('sample_code_2' = 'sample_code')) |>
  dplyr::select(sample_id_num_1, sample_code_1,sample_code_2, fraction, wunifrac_distance, sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(lag = as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))

unifrac_tibble_3 <- unifrac_tibble |>
  dplyr::filter(str_detect(sample_id_1, '_3_') &
                  str_detect(sample_id_2, '_3_')) |>
  separate(sample_id_1, sep = '_', into = c('sample_code_1', 'fraction', 'seqid'), remove = F) |>
  separate(sample_id_2, sep = '_', into = c('sample_code_2', 'fraction', 'seqid'), remove = F) |>
  left_join(m_02_ed, by = c('sample_code_1' = 'sample_code')) |>
  dplyr::select(sample_id_num_1 = sample_id_num, sample_code_1, sample_code_2, fraction, wunifrac_distance ) |>
  left_join(m_02_ed, by = c('sample_code_2' = 'sample_code')) |>
  dplyr::select(sample_id_num_1, sample_code_1,sample_code_2, fraction, wunifrac_distance, sample_id_num_2 = sample_id_num) |>
  dplyr::mutate(lag = as.numeric(sample_id_num_1) - as.numeric(sample_id_num_2))

unifrac_tibble_02 %$%
  range(lag)

unifrac_tibble_3 %$%
  range(lag)

wunifrac_lag_all <- unifrac_tibble_02 |>
  bind_rows(unifrac_tibble_3) |>
  dplyr::filter(lag > 0)

wunifrac_lag_plot <- wunifrac_lag_all |>
  ggplot(aes(lag, wunifrac_distance))+
  scale_x_continuous(breaks = c(12, 24, 36, 48, 12*5, 12*6, 12*7, 12*8, 12*9, 120), labels = c(1, 2, 3,4, 5, 6, 7, 8,  9, 10))+
  geom_point(aes(shape = fraction), alpha = 0.08)+
  scale_shape_discrete(labels = labs_fraction)+
  labs(y = 'wUNIFRAC', x = 'Number of years between samples', shape = '', color = '', fill = '')+
  geom_line(data = wunifrac_lag_all |>
              group_by(fraction, lag) |>
              dplyr::reframe(mean = mean(wunifrac_distance)), aes(lag, mean, group = fraction, color = fraction), linewidth = 1)+
  scale_color_manual(values = palette_fraction, labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction, labels = labs_fraction)+
  guides(shape = guide_legend(ncol = 2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom')

wunifrac_lag_plot

# Save the composition as a PDF
# ggsave(filename = 'wUNIFRAC_recurrence_bbmo_plot.pdf',
#        plot = wunifrac_lag_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 120, units = 'mm')

## environmental recurrency based on euclidean distance ----
euclidean_distance

m_02_ed <- m_02 |>
  dplyr::select(sample_id,  decimal_date, sample_id_num)

euclidean_distance_lag <- euclidean_distance |>
  dplyr::mutate(lag = as.numeric(sample_id_num) - as.numeric(sample_id_num_2)) |>
  dplyr::filter(as.numeric(lag) > 0)

euclidean_env_distance_plot <- euclidean_distance_lag |>
  dplyr::filter(as.numeric(lag) > 0) |>
  ggplot(aes(lag, euclidean_distance))+
  scale_x_continuous(breaks = c(12, 24, 36, 48, 12*5, 12*6, 12*7, 12*8, 12*9, 120), labels = c(1, 2, 3,4, 5, 6, 7, 8,  9, 10))+
  geom_point(aes(), alpha = 0.08)+
  labs(y = 'Euclidean distance', x = 'Number of years between samples', shape = '', color = '', fill = '')+
  geom_line(data = euclidean_distance_lag |>
              group_by( lag) |>
              dplyr::reframe(mean = mean(euclidean_distance)), aes(lag, mean), linewidth = 1)+
  scale_color_manual(values = palette_fraction_env, labels = labs_fraction)+
  scale_fill_manual(values = palette_fraction_env, labels = labs_fraction)+
  guides(shape = guide_legend(ncol = 2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        legend.position = 'bottom')

euclidean_env_distance_plot

# Save the composition as a PDF
# ggsave(filename = 'euclidean_env_distance_plot.pdf',
#        plot = euclidean_env_distance_plot,
#        path = 'Results/Figures/',
#        width = 180, height = 120, units = 'mm')