# Code for differential gene correlation analysis of GTEx Human thyroid RNA-Seq data
# Author: Ameya Kulkarni
# Date: October 31st 2018

# Load libraries
library(DGCA, quietly = T)
library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(grex, quietly = T)
library(repmis, quietly = T)

thyroidurl <- "https://github.com/ameya225/Aging_DGCA/blob/master/Data/thyroid.RData?raw=true"

# Load expression data, sample data and age-stratified data
source_data(thyroidurl)

# Stratify expression data by age groups
thyroid_young_edat <- cbind(thyroid_edat[,which(colnames(thyroid_edat)%in%rownames(thyroid_splitbyage$`20-29`))], thyroid_edat[,which(colnames(thyroid_edat)%in%rownames(thyroid_splitbyage$`30-39`))])
thyroid_mid_edat <- cbind(thyroid_edat[,which(colnames(thyroid_edat)%in%rownames(thyroid_splitbyage$`40-49`))], thyroid_edat[,which(colnames(thyroid_edat)%in%rownames(thyroid_splitbyage$`50-59`))])
thyroid_old_edat <- cbind(thyroid_edat[,which(colnames(thyroid_edat)%in%rownames(thyroid_splitbyage$`60-69`))], thyroid_edat[,which(colnames(thyroid_edat)%in%rownames(thyroid_splitbyage$`70-79`))])
thyroid <- cbind(thyroid_young_edat, thyroid_mid_edat, thyroid_old_edat)

# Filter expression data for low expression and low dispersion (3030x350)
thyroid_filt = DGCA::filterGenes(thyroid, filterTypes = "central", "dispersion", sequential = T, filterCentralType = "median", filterCentralPercentile = 0.90, filterDispersionType = "cv", filterDispersionPercentile = 0.90)
rownames(thyroid_filt) <- unlist(lapply(strsplit(rownames(thyroid_filt), split = "[.]"), function(x){x[1]}))
thyroid_filt <- cbind(rownames(thyroid_filt), thyroid_filt)
thyroid_genes <- grex(cleanid(rownames(thyroid_filt)))
thyroid_filt <- as.data.frame(thyroid_filt)
colnames(thyroid_filt)[1] <- "ensembl_id"
thyroid_genesym <- dplyr::left_join(thyroid_filt, thyroid_genes, by="ensembl_id")
thyroid_genesym <- thyroid_genesym[complete.cases(thyroid_genesym),]
rownames(thyroid_genesym) <- thyroid_genesym$hgnc_symbol
thyroid_filt_mat <- thyroid_genesym[,-c(1,(ncol(thyroid_genesym)-5):ncol(thyroid_genesym))]
thyroid_filt_mat <- data.matrix(thyroid_filt_mat)

# Characterize subjects into age groups
colnames(thyroid_filt_mat)[1:42] <- paste(colnames(thyroid_filt_mat)[1:42], "1", sep = "_")
colnames(thyroid_filt_mat)[43:221] <- paste(colnames(thyroid_filt_mat)[43:221], "2", sep = "_")
colnames(thyroid_filt_mat)[222:323] <- paste(colnames(thyroid_filt_mat)[222:323], "3", sep = "_")
thyroid_groups <- factor(unlist(lapply(strsplit(colnames(thyroid_filt_mat), split = "_"), function(x){x[2]})))

thyroid_groups <- gsub(pattern = "1", replacement = "Y", x = thyroid_groups)
thyroid_groups <- gsub(pattern = "2", replacement = "M", x = thyroid_groups)
thyroid_groups <- gsub(pattern = "3", replacement = "O", x = thyroid_groups)

# Make a design matrix
thyroid_des <- matrix(nrow = ncol(thyroid_filt_mat), ncol = 3)
thyroid_des <- as.data.frame(thyroid_des)
thyroid_des[grep("_1", colnames(thyroid_filt_mat)),1] <- rep(1,length(grep("_1", colnames(thyroid_filt_mat))))
thyroid_des[grep("_2", colnames(thyroid_filt_mat)),2] <- rep(1,length(grep("_2", colnames(thyroid_filt_mat))))
thyroid_des[grep("_3", colnames(thyroid_filt_mat)),3] <- rep(1,length(grep("_3", colnames(thyroid_filt_mat))))
thyroid_des[is.na(thyroid_des)] <- 0
thyroid_des <- as.matrix(thyroid_des)
colnames(thyroid_des) <- c("Y", "M", "O")

# Correlations and significance in each condition 
thyroid_cor_res = getCors(inputMat = thyroid_filt_mat, design = thyroid_des, corrType = "spearman")

# Differential correlations between conditions
thyroid_dcpairs_YO = pairwiseDCor(thyroid_cor_res, compare = c("Y", "O"))
thyroid_dcpairs_YM = pairwiseDCor(thyroid_cor_res, compare = c("Y", "M"))
thyroid_dcpairs_MO = pairwiseDCor(thyroid_cor_res, compare = c("M", "O"))

# Extracting top differentially correlated pairs between conditions
thyroid_dd_pairs_YO = na.omit(dcTopPairs(thyroid_dcpairs_YO, nPairs = (((nrow(thyroid_filt_mat)*(nrow(thyroid_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
thyroid_dd_pairs_YO <- thyroid_dd_pairs_YO[thyroid_dd_pairs_YO$Classes!="NonSig",] 
thyroid_dd_pairs_YO_rel <- thyroid_dd_pairs_YO[thyroid_dd_pairs_YO$Classes!="+/+",] 
thyroid_dd_pairs_YO_rel <- thyroid_dd_pairs_YO_rel[thyroid_dd_pairs_YO_rel$Classes!="-/-",] 
thyroid_dd_pairs_YO_rel <- thyroid_dd_pairs_YO_rel[thyroid_dd_pairs_YO_rel$Classes!="0/0",]

thyroid_dd_pairs_MO = na.omit(dcTopPairs(thyroid_dcpairs_MO, nPairs = (((nrow(thyroid_filt_mat)*(nrow(thyroid_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
thyroid_dd_pairs_MO <- thyroid_dd_pairs_MO[thyroid_dd_pairs_MO$Classes!="NonSig",] 
thyroid_dd_pairs_MO_rel <- thyroid_dd_pairs_MO[thyroid_dd_pairs_MO$Classes!="+/+",] 
thyroid_dd_pairs_MO_rel <- thyroid_dd_pairs_MO_rel[thyroid_dd_pairs_MO_rel$Classes!="-/-",] 
thyroid_dd_pairs_MO_rel <- thyroid_dd_pairs_MO_rel[thyroid_dd_pairs_MO_rel$Classes!="0/0",]

thyroid_dd_pairs_YM = na.omit(dcTopPairs(thyroid_dcpairs_YM, nPairs = (((nrow(thyroid_filt_mat)*(nrow(thyroid_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
thyroid_dd_pairs_YM <- thyroid_dd_pairs_YM[thyroid_dd_pairs_YM$Classes!="NonSig",] 
thyroid_dd_pairs_YM_rel <- thyroid_dd_pairs_YM[thyroid_dd_pairs_YM$Classes!="+/+",] 
thyroid_dd_pairs_YM_rel <- thyroid_dd_pairs_YM_rel[thyroid_dd_pairs_YM_rel$Classes!="-/-",] 
thyroid_dd_pairs_YM_rel <- thyroid_dd_pairs_YM_rel[thyroid_dd_pairs_YM_rel$Classes!="0/0",]

# sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     
# 
# other attached packages:
#   [1] repmis_0.5    grex_1.8      dplyr_0.7.6   ggplot2_3.0.0
# [5] DGCA_1.0.1   
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.3.1            Biobase_2.40.0       
# [3] dynamicTreeCut_1.63-1 bit64_0.9-7          
# [5] splines_3.5.1         foreach_1.4.4        
# [7] R.utils_2.7.0         Formula_1.2-3        
# [9] assertthat_0.2.0      stats4_3.5.1         
# [11] latticeExtra_0.6-28   blob_1.1.1           
# [13] fit.models_0.5-14     yaml_2.2.0           
# [15] robustbase_0.93-2     impute_1.54.0        
# [17] pillar_1.3.0          RSQLite_2.1.1        
# [19] backports_1.1.2       lattice_0.20-35      
# [21] glue_1.3.0            digest_0.6.16        
# [23] RColorBrewer_1.1-2    checkmate_1.8.5      
# [25] colorspace_1.3-2      R.oo_1.22.0          
# [27] htmltools_0.3.6       preprocessCore_1.42.0
# [29] Matrix_1.2-14         plyr_1.8.4           
# [31] pcaPP_1.9-73          pkgconfig_2.0.2      
# [33] purrr_0.2.5           GO.db_3.6.0          
# [35] mvtnorm_1.0-8         scales_1.0.0         
# [37] htmlTable_1.12        tibble_1.4.2         
# [39] IRanges_2.14.10       withr_2.1.2          
# [41] fastcluster_1.1.25    nnet_7.3-12          
# [43] BiocGenerics_0.26.0   lazyeval_0.2.1       
# [45] survival_2.42-6       magrittr_1.5         
# [47] crayon_1.3.4          memoise_1.1.0        
# [49] R.methodsS3_1.7.1     R.cache_0.13.0       
# [51] doParallel_1.0.11     MASS_7.3-50          
# [53] foreign_0.8-71        tools_3.5.1          
# [55] data.table_1.11.4     matrixStats_0.54.0   
# [57] stringr_1.3.1         S4Vectors_0.18.3     
# [59] munsell_0.5.0         cluster_2.0.7-1      
# [61] AnnotationDbi_1.42.1  bindrcpp_0.2.2       
# [63] compiler_3.5.1        rlang_0.2.2          
# [65] grid_3.5.1            iterators_1.0.10     
# [67] rstudioapi_0.7        htmlwidgets_1.2      
# [69] WGCNA_1.63            robust_0.4-18        
# [71] base64enc_0.1-3       gtable_0.2.0         
# [73] codetools_0.2-15      curl_3.2             
# [75] DBI_1.0.0             rrcov_1.4-4          
# [77] R6_2.2.2              gridExtra_2.3        
# [79] knitr_1.20            bit_1.1-14           
# [81] bindr_0.1.1           Hmisc_4.1-1          
# [83] stringi_1.2.4         parallel_3.5.1       
# [85] Rcpp_0.12.18          rpart_4.1-13         
# [87] acepack_1.4.1         DEoptimR_1.0-8       
# [89] tidyselect_0.2.4    


