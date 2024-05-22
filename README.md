## Repository Description

This repository includes source codes for the study **"Exploring sexual dimorphism of Japanese quail based on RNA sequencing analysis"** submitted to the journal.

## Source Code

- **Package_install.R**: Script to install all the packages needed for RNA-Seq analysis. Please run this first.
- **Characterization_quail_brain.R**: Chapter 1: identifying the characterization of brain tissue in Japanese quail.
- **Characterization_quail_gonadal_tissues.R**: Chapter 2: identifying the characterization of gonadal tissues.
- **Investigation_tissue_specificity_genes.R**: Chapter 3: for investigating tissue-specific expression of sex-biased genes in Japanese quail.
- **Sex_determined_genes_during_embryonic_stages.R**: Chapter 4: investigating sex-determined genes during embryonic development of Japanese quail.

## Data

- **mart_export.txt (tab-separated format)**: Detailed information on the corresponding gene annotation obtained from Ensembl Biomart of Coturnix_japonica_2.0.
- **metadata.txt (tab-separated format)**: Detailed information of RNA-Seq samples of Japanese quail.

## Folder Explanation

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
  - **Figure5A.RData**: R data for plotting Figure 5A, showing the overall gene expression pattern in the embryonic stage of chicken embryo.
  - **Figure5B.RData**: R data for plotting Figure 5B, showing the overall gene expression pattern in the embryonic stage of quail embryo.
  - **Figure5C_chicken.RData**: Showing DMRT1 gene expression pattern in the embryonic stage of chicken embryo.
  - **Figure5C_quail.RData**: Showing DMRT1 gene expression pattern in the embryonic stage of quail embryo.
  - **Figure5D_data.txt**: Heatmap data for drawing the sex-determined genes expression pattern in chicken embryo.
  - **Figure5E_data.txt**: Heatmap data for drawing the sex-determined genes expression pattern in quail embryo.
  - **Figure5F_data.txt**: Data for comparison of gene expression in the Z chromosome between female and male chicken embryos.
  - **Figure5G_data.txt**: Data for comparison of gene expression in the Z chromosome between female and male quail embryos.
  - **Figure5I.RData**: Data for investigating the gene expression of the Z chromosome in quail embryo.
  - **Figure5K_L_data.csv**: Data for evaluation of gene annotation of quail and chicken embryos.
  - **mart_export_galgal6.txt**: Data for annotation of the chicken GRCg6a genome.
  - **Metadata.txt**: Information of chicken and quail embryo RNA-Seq data.
  - **Orthologue_filtered_galgal6.txt**: Orthologue gene list between Coturnix_japonica_2.0 and GRCg6a.

- **Quantification_quail**: Folder containing the quantification result of Japanese quail RNA-Seq data.
- **Quantification_chicken_embryo**: Folder containing the quantification result of chicken embryo RNA-Seq data.

