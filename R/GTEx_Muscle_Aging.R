# Code for differential gene correlation analysis of GTEx Human muscle RNA-Seq data
# Author: Ameya Kulkarni
# Date: October 31st 2018

# Load libraries
library(DGCA, quietly = T)
library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(grex, quietly = T)
library(repmis, quietly = T)

# Load expression data, sample data and age-stratified data
muscleurl <- "https://github.com/ameya225/Aging_DGCA/blob/master/Data/muscle.RData?raw=true"
source_data(muscleurl)

# Stratify expression data by age groups
# Young = 20-39 years old
# Middle = 40-59 years old
# Old = 60-79 years old
muscle_young_edat <- cbind(muscle_edat[,which(colnames(muscle_edat)%in%rownames(muscle_splitbyage$`20-29`))], muscle_edat[,which(colnames(muscle_edat)%in%rownames(muscle_splitbyage$`30-39`))])
muscle_mid_edat <- cbind(muscle_edat[,which(colnames(muscle_edat)%in%rownames(muscle_splitbyage$`40-49`))], muscle_edat[,which(colnames(muscle_edat)%in%rownames(muscle_splitbyage$`50-59`))])
muscle_old_edat <- cbind(muscle_edat[,which(colnames(muscle_edat)%in%rownames(muscle_splitbyage$`60-69`))], muscle_edat[,which(colnames(muscle_edat)%in%rownames(muscle_splitbyage$`70-79`))])
muscle <- cbind(muscle_young_edat, muscle_mid_edat, muscle_old_edat)

# Filter expression data for low expression and low dispersion (3030x350)
muscle_filt = DGCA::filterGenes(muscle, filterTypes = "central", "dispersion", sequential = T, filterCentralType = "median", filterCentralPercentile = 0.90, filterDispersionType = "cv", filterDispersionPercentile = 0.90)
rownames(muscle_filt) <- unlist(lapply(strsplit(rownames(muscle_filt), split = "[.]"), function(x){x[1]}))
muscle_filt <- cbind(rownames(muscle_filt), muscle_filt)
muscle_genes <- grex(cleanid(rownames(muscle_filt)))
muscle_filt <- as.data.frame(muscle_filt)
colnames(muscle_filt)[1] <- "ensembl_id"
muscle_genesym <- dplyr::left_join(muscle_filt, muscle_genes, by="ensembl_id")
muscle_genesym <- muscle_genesym[complete.cases(muscle_genesym),]
rownames(muscle_genesym) <- muscle_genesym$hgnc_symbol
muscle_filt_mat <- muscle_genesym[,-c(1,(ncol(muscle_genesym)-5):ncol(muscle_genesym))]
muscle_filt_mat <- data.matrix(muscle_filt_mat)

# Characterize subjects into age groups
colnames(muscle_filt_mat)[1:70] <- paste(colnames(muscle_filt_mat)[1:70], "1", sep = "_")
colnames(muscle_filt_mat)[71:288] <- paste(colnames(muscle_filt_mat)[71:288], "2", sep = "_")
colnames(muscle_filt_mat)[288:430] <- paste(colnames(muscle_filt_mat)[288:430], "3", sep = "_")
muscle_groups <- factor(unlist(lapply(strsplit(colnames(muscle_filt_mat), split = "_"), function(x){x[2]})))

muscle_groups <- gsub(pattern = "1", replacement = "Y", x = muscle_groups)
muscle_groups <- gsub(pattern = "2", replacement = "M", x = muscle_groups)
muscle_groups <- gsub(pattern = "3", replacement = "O", x = muscle_groups)

# Make a design matrix
muscle_des <- matrix(nrow = ncol(muscle_filt_mat), ncol = 3)
muscle_des <- as.data.frame(muscle_des)
muscle_des[grep("_1", colnames(muscle_filt_mat)),1] <- rep(1,length(grep("_1", colnames(muscle_filt_mat))))
muscle_des[grep("_2", colnames(muscle_filt_mat)),2] <- rep(1,length(grep("_2", colnames(muscle_filt_mat))))
muscle_des[grep("_3", colnames(muscle_filt_mat)),3] <- rep(1,length(grep("_3", colnames(muscle_filt_mat))))
muscle_des[is.na(muscle_des)] <- 0
muscle_des <- as.matrix(muscle_des)
colnames(muscle_des) <- c("Y", "M", "O")


# Sanity check
# boxplot(muscle_filt_mat[,1:100])

# Correlations and significance in each condition 
muscle_cor_res = getCors(inputMat = muscle_filt_mat, design = muscle_des, corrType = "spearman")

# Differential correlations between conditions
muscle_dcpairs_YO = pairwiseDCor(muscle_cor_res, compare = c("Y", "O"))
muscle_dcpairs_YM = pairwiseDCor(muscle_cor_res, compare = c("Y", "M"))
muscle_dcpairs_MO = pairwiseDCor(muscle_cor_res, compare = c("M", "O"))

# Extracting top differentially correlated pairs between conditions
muscle_dd_pairs_YO = na.omit(dcTopPairs(muscle_dcpairs_YO, nPairs = (((nrow(muscle_filt_mat)*(nrow(muscle_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
muscle_dd_pairs_YO <- muscle_dd_pairs_YO[muscle_dd_pairs_YO$Classes!="NonSig",] 
muscle_dd_pairs_YO_rel <- muscle_dd_pairs_YO[muscle_dd_pairs_YO$Classes!="+/+",] 
muscle_dd_pairs_YO_rel <- muscle_dd_pairs_YO_rel[muscle_dd_pairs_YO_rel$Classes!="-/-",] 
muscle_dd_pairs_YO_rel <- muscle_dd_pairs_YO_rel[muscle_dd_pairs_YO_rel$Classes!="0/0",]

muscle_dd_pairs_MO = na.omit(dcTopPairs(muscle_dcpairs_MO, nPairs = (((nrow(muscle_filt_mat)*(nrow(muscle_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
muscle_dd_pairs_MO <- muscle_dd_pairs_MO[muscle_dd_pairs_MO$Classes!="NonSig",] 
muscle_dd_pairs_MO_rel <- muscle_dd_pairs_MO[muscle_dd_pairs_MO$Classes!="+/+",] 
muscle_dd_pairs_MO_rel <- muscle_dd_pairs_MO_rel[muscle_dd_pairs_MO_rel$Classes!="-/-",] 
muscle_dd_pairs_MO_rel <- muscle_dd_pairs_MO_rel[muscle_dd_pairs_MO_rel$Classes!="0/0",]

muscle_dd_pairs_YM = na.omit(dcTopPairs(muscle_dcpairs_YM, nPairs = (((nrow(muscle_filt_mat)*(nrow(muscle_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
muscle_dd_pairs_YM <- muscle_dd_pairs_YM[muscle_dd_pairs_YM$Classes!="NonSig",] 
muscle_dd_pairs_YM_rel <- muscle_dd_pairs_YM[muscle_dd_pairs_YM$Classes!="+/+",] 
muscle_dd_pairs_YM_rel <- muscle_dd_pairs_YM_rel[muscle_dd_pairs_YM_rel$Classes!="-/-",] 
muscle_dd_pairs_YM_rel <- muscle_dd_pairs_YM_rel[muscle_dd_pairs_YM_rel$Classes!="0/0",]

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

