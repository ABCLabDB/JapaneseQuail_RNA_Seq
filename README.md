## Repository description

This repository includes source codes for the study **"Exploring sexual dimorphism of Japanese quail based on RNA sequencing analysis"** submitted to the journal.

## Source code

- **Package_install.R**: Script to install all the packages needed for RNA-Seq analysis. Please run this first.
- **Figure2.R**: R code for identifying the characterization of brain tissue in Japanese quail.
- **Figure3.R**: R code for identifying the characterization of gonadal tissues.
- **Figure4.R**: R code for investigating tissue-specific expression of sex-biased genes in Japanese quail.
- **Figure5.R**: R code for investigating sex determination during embryonic development of Japanese quail.

## Data

- **mart_export.txt (tab-separated format)**: Detailed information on the corresponding gene annotation obtained from Ensembl Biomart of Coturnix_japonica_2.0.
- **metadata.txt (tab-separated format)**: Detailed information of RNA-Seq samples of Japanese quail.

## Folder explanation

- **Figure2_data**: 
  - **David_result_brain.txt**: The result of functional analysis of sex-biased genes in brain tissue.
  - **metadata_FB_MB.txt**: Detailed information of RNA-Seq samples of brain tissue in Japanese quail.
  - **quail_annotation.txt**: Gene annotation file of Japanese quail.
  - **Galgal6_annotation.txt**: Gene annotation file of chicken.
  - **Result_FB_MB.txt**: Information on sex-biased genes in brain tissue of Japanese quail.
  - **WD.txt**: Gene list related to WDR function.

- **Figure3_data**: 
  - **Female_up.txt**: Clustering of functional analysis of female-biased genes in the reproductive organ of Japanese quail.
  - **Female_up_bar_plot_data.txt**: Data for plotting the bar plot of clustering of functional analysis with female-biased genes.
  - **Male_up.txt**: Clustering of functional analysis of male-biased genes in the reproductive organ of Japanese quail.
  - **Male_up_bar_plot_data.txt**: Data for plotting the bar plot of clustering of functional analysis with male-biased genes.
  - **metadata_O_T.txt**: Detailed information of RNA-Seq samples of gonadal tissues in Japanese quail.
  - **Result_O_T.txt**: Information on sex-biased genes in gonadal tissues of Japanese quail.
  - **O_T.RData**: R data for saving the result of sex-biased genes analysis in gonadal tissues of Japanese quail.

- **Figure4_data**: 
  - **Figure4_data.txt**: Data for finding tissue-specific genes.

- **Figure5_data**: 
  - **Sex_biased_genes**: Folder containing the results of statistical tests for finding sexual dimorphism in each embryonic stage of chicken embryo.
  - **Figure5A.RData**: R data for plotting Figure 5A, showing overall gene expression pattern in the embryonic stage of chicken embryo.
  - **Figure5B.RData**: R data for plotting Figure 5B, showing overall gene expression pattern in the embryonic stage of quail embryo.
  - **Figure5C_chicken.RData**: Showing DMRT1 gene expression pattern in the embryonic stage of chicken embryo.
  - **Figure5C_quail.RData**: Showing DMRT1 gene expression pattern in the embryonic stage of quail embryo.
  - **Figure5D_data.txt**: Heatmap data for drawing the sex-determined genes expression pattern in chicken embryo.
  - **Figure5E_data.txt**: Heatmap data for drawing the sex-determined genes expression pattern in quail embryo.
  - **Figure5F_data.txt**: Data for comparison of gene expression in the Z chromosome between female and male chicken embryos.
  - **Figure5G_data.txt**: Data for comparison of gene expression in the Z chromosome between female and male quail embryos.
  - **Figure5I.RData**: Data for investigating the gene expression of the Z chromosome in quail embryo.
  - **Figure5K_L_data.csv**: Data for evaluation of gene annotation of quail and chicken embryo.

- **Quantification_quail**: Folder containing the quantification result of Japanese quail RNA-Seq data.

