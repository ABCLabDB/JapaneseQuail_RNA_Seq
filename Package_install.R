setRepositories(ind=1:7)
packages <- c("ggplot2", "ggsci", "ggpubr", "data.table", "edgeR", "DESeq2", 
              "ggrepel", "dplyr", "RColorBrewer", "GenomicRanges", "Gviz", 
              "tidyverse", "viridis", "randomcoloR", "ggalt", "stringr", "tidyr", 
              "ggradar")


installed_packages <- rownames(installed.packages())

for (package in packages) {
  if (!(package %in% installed_packages)) {
    install.packages(package, dependencies = TRUE)
  }
}
