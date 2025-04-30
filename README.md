# Breast-cancer-research
Identification of Hub Genes, Pathways, and Potential Drug Repurposing for Breast Cancer Through WGCNA and Integrative Bioinformatics Analysis

Dependencies:                                                 
R/4.4.2                                              
RStudio/2024.12.0+467

File Descriptions

Files in /datasets:

GSE42568.csv, GSE61304.csv, and GSE29431.csv: Processed breast cancer gene expression
datasets in CSV format for differential gene expression analysis.

missing_free_phenoData_GSE42568.csv: Phenotypic dataset after excluding samples with
missing information of GSE42568.

missing_free_data_GSE42568.csv: Prepared gene expression dataset according to
phenotypic dataset of GSE42568.

pred_accuracy.csv: Filtered and formatted dataset to evaluate the predictive accuracy of ten
HubGs from these three datasets.

Files in /R_Code:

All source code for principal component analysis (PCA), cluster dendrogram (CD),
differentially expressed genes (DEGs) identification, weighted gene co-expression network
analysis (WGCNA), and predictive accuracy analysis of hub genes (HubGs) through SVM, RF,
and XGBoost models.
