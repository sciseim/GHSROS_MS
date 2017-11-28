# GHSROS_MS
Supplemental information, code, and data for the MS: 
Patrick B. Thomas, Penny L. Jeffery, Manuel D. Gahete, Eliza Whiteside, Carina Walpole, Michelle Maugham, Lidija Jovanovic, Jennifer H. Gunter, Elizabeth D. Williams, Colleen C. Nelson, Adrian C. Herington, Raúl M. Luque, Rakesh Veedu, Lisa K. Chopin, Inge Seim. *The long non-coding RNA ¬¬¬¬¬¬GHSROS mediates expression of genes associated with aggressive prostate cancer and promotes tumor progression*.

*Summary*:
See individual scripts for additional information.

Scripts 1A-C is the code for the interrogation of your gene of interest using Affymetrix exon array data. Raw data (CEL) are input into the script `1A-SCAN.UPC.R`. [SCAN.UPC](http://www.pnas.org/content/110/44/17778.long) outputs standardised expression values (UPC value), ranging from 0 to 1, which indicate whether a gene is actively transcribed in a sample of interest: higher values indicate that a gene is ‘active’. UPC scores are platform-independent and allow cross-experimental and cross-platform integration. `1B-UPC_parse.R` parses the output from SCAN.UPC, and `1C-UPC_scatter_plot.R` draws a figure showing the expression of your transcript of interest in the data sets examined.



 
 
 
 

 
Scripts 2-4 will interrogate [TCGA](https://cancergenome.nih.gov) data sets:
- `2A-H, and 2x`. Reveals whether your genes of interest (e.g. a gene signature) are associated with overall survival (OS) and/or disease-free survival (DFS). Also assesses the performance of a gene signature compared to random gene signatures of the same size (`2G-gene signature VS random signatures.R`) and known gene signatures (`2H-gene signature VS known signatures.R`).
- `3-survival_analysis_TCGA-datasets.R`. Performs survival analysis of a gene signature in 32 TCGA cancer data sets. Clinical information from [cBioPortal](http://www.cbioportal.org/data_sets.jsp) are stored in ./data/cBioPortal (last accessed 08.05.16).
- `4-TCGA-PRAD heatmap.R.R`. Draws a heatmap of your gene signature in a TCGA data set. Requires [heatmap.3.R](https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R).







[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.163883.svg)](https://doi.org/10.5281/zenodo.163883)
