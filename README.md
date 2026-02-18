# Bloomers_10y_BBMO

Studying blooming events on a 10years timeseries [Blanes Bay Microbial Observatory](http://bbmo.icm.csic.es/) which corresponds to the following manuscript:

The natural history of bacterial bloomers in a decade-long time series

## Authors:

Ona Deulofeu-Capo, Carmen García-Comas, Xavier Rey-Velasco, Adrià Auladell, Ramiro Logares, Esther Garcés, Isabel Ferrera, Olga Sánchez, Josep M Gasol and Marta Sebastián

# Abstract
Prokaryotic bloomers—populations that experience rapid and significant increases in abundance in response to environmental triggers—briefly dominate marine microbial communities, potentially impacting the ecosystem by channeling large amounts of nutrients and affecting carbon fluxes. Due to their ephemeral nature, prokaryotic bloomers are challenging to capture, and it remains unknown whether they are restricted to specific taxonomic groups or whether they exhibit recurrent patterns. Here, we analyzed a decade-long time series from the Blanes Bay Microbial Observatory (BBMO, NW Mediterranean Sea) to investigate bacterial bloomers in two size fractions. We identified 57 Amplicon Sequence Variants (ASVs), less than 1% of the total bacterial richness, exhibiting either recurrent or chaotic blooming-like behavior. The ability to bloom was observed across diverse phyla, though some taxonomic coherence appeared within families containing multiple blooming taxa. In our monthly sampling, bloom events were detected on average 4.6±1.9 times per year (both size fractions considered), and weakly related to biological and physicochemical variables once seasonality was accounted for, likely because of the low sampling resolution. Yet, a notable shift in the blooming community in the 3-20 μm fraction reflected ecosystem disturbances associated to the nearby harbor restoration, suggesting bloomers may serve as sentinels for ecosystem disturbances. Combined analyses with metagenomic data showed that blooms led to marked shifts in the community functional potential. Overall, our findings underscore the importance of investigating bloom dynamics to understand microbial contributions to biogeochemical cycles and stress the need for higher-frequency sampling to accurately capture these transient but ecologically relevant events.

### Additional Notes
We updated the taxonomy to the 138.1 SILVA database, the new correct file for the taxonomy
is located in 03_tax_assignation/ folder.

## Reproducibility and Data Availability
Code for the MDR-S map https://github.com/biozoo/MDR_S-map
The raw metabarcoding data (PRJEB38773) and raw metagenomic data (PRJEB48035) are available in the European Nucleotide Archive. 
Genes tables are available in: https://zenodo.org/records/17183573 
Find bloomer's funciton: https://github.com/EcologyR/Bloomers

### -----------  Packages version information -------
# sessionInfo()
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS 15.3.2
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] grateful_0.2.11      zCompositions_1.4.1  truncnorm_1.0-9      NADA_1.6-1.1         survival_3.5-7       MASS_7.3-60          EcolUtils_0.1        vegan_2.6-4          lattice_0.22-5      
# [10] permute_0.9-7        speedyseq_0.5.3.9018 scales_1.3.0         magrittr_2.0.3       Bloomers_0.0.0.9000  janitor_2.2.0        lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1       
# [19] dplyr_1.1.4          purrr_1.0.2          readr_2.1.4          tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0      readxl_1.4.3         phyloseq_1.42.0   
