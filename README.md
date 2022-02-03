GWAS pipeline is divided into three components.

1. Initial Quality control steps - initial_qc.R
2. Ancestry matching and Population outliers using IBD and PCA - ibd_pca_after_initial_qc.R
3. Phasing and Imputation - imputation.R

All these three steps are supported using custom functions from gwasfunctions.r and common_snps_1kg_othercohort.pl scripts.
