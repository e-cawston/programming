#!/usr/bin/env Rscript
# Install missing packages

if (!require('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

BiocManager::install(c('ggplot2', 'sleuth', 'UpSetR', 'RankProd', 'openxlsx'),
                     update = FALSE, ask = FALSE)
