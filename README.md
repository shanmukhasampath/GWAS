GWAS pipeline is divided into three components.

1. Initial Quality control steps - Initial_Quality_Control.r
2. Ancestry matching and Population outliers using IBD and PCA - AncestryMatching_OutlierDectection.r
3. Phasing and Imputation - Imputation.r

All these three steps are supported using custom functions from gwasfunctions.r and common_snps_1kg_othercohort.pl scripts.
