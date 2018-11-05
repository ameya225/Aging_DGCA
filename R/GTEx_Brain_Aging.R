# Code for differential gene correlation analysis of GTEx Human brain RNA-Seq data
# Author: Ameya Kulkarni
# Date: October 31st 2018

# Load libraries
library(DGCA, quietly = T)
library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(grex, quietly = T)
library(repmis, quietly = T)

brainurl <- "https://github.com/ameya225/Aging_DGCA/blob/master/Data/brain.RData?raw=true"

# Load expression data, sample data and age-stratified data
source_data(brainurl)

# brain0 = Cortex + others
# brain1 = cerebellum
# brain2 = basal ganglia

##########  BRAIN0  ##############

brain0_young_edat <- cbind(brain0_edat[,which(colnames(brain0_edat)%in%rownames(brain0_splitbyage$`20-29`))], brain0_edat[,which(colnames(brain0_edat)%in%rownames(brain0_splitbyage$`30-39`))])
brain0_mid_edat <- cbind(brain0_edat[,which(colnames(brain0_edat)%in%rownames(brain0_splitbyage$`40-49`))], brain0_edat[,which(colnames(brain0_edat)%in%rownames(brain0_splitbyage$`50-59`))])
brain0_old_edat <- cbind(brain0_edat[,which(colnames(brain0_edat)%in%rownames(brain0_splitbyage$`60-69`))], brain0_edat[,which(colnames(brain0_edat)%in%rownames(brain0_splitbyage$`70-79`))])
brain0 <- cbind(brain0_young_edat, brain0_mid_edat, brain0_old_edat)

# Filter expression data for low expression and low dispersion (3030x350)
brain0_filt = DGCA::filterGenes(brain0, filterTypes = "central", "dispersion", sequential = T, filterCentralType = "median", filterCentralPercentile = 0.90, filterDispersionType = "cv", filterDispersionPercentile = 0.90)
rownames(brain0_filt) <- unlist(lapply(strsplit(rownames(brain0_filt), split = "[.]"), function(x){x[1]}))
brain0_filt <- cbind(rownames(brain0_filt), brain0_filt)
brain0_genes <- grex(cleanid(rownames(brain0_filt)))
brain0_filt <- as.data.frame(brain0_filt)
colnames(brain0_filt)[1] <- "ensembl_id"
brain0_genesym <- dplyr::left_join(brain0_filt, brain0_genes, by="ensembl_id")
brain0_genesym <- brain0_genesym[complete.cases(brain0_genesym),]
rownames(brain0_genesym) <- brain0_genesym$hgnc_symbol
brain0_filt_mat <- brain0_genesym[,-c(1,(ncol(brain0_genesym)-5):ncol(brain0_genesym))]
brain0_filt_mat <- data.matrix(brain0_filt_mat)

# Characterize subjects into age groups
colnames(brain0_filt_mat)[1:36] <- paste(colnames(brain0_filt_mat)[1:36], "1", sep = "_")
colnames(brain0_filt_mat)[37:346] <- paste(colnames(brain0_filt_mat)[37:346], "2", sep = "_")
colnames(brain0_filt_mat)[347:702] <- paste(colnames(brain0_filt_mat)[347:702], "3", sep = "_")
brain0_groups <- factor(unlist(lapply(strsplit(colnames(brain0_filt_mat), split = "_"), function(x){x[2]})))

brain0_groups <- gsub(pattern = "1", replacement = "Y", x = brain0_groups)
brain0_groups <- gsub(pattern = "2", replacement = "M", x = brain0_groups)
brain0_groups <- gsub(pattern = "3", replacement = "O", x = brain0_groups)

# Make a design matrix
brain0_des <- matrix(nrow = ncol(brain0_filt_mat), ncol = 3)
brain0_des <- as.data.frame(brain0_des)
brain0_des[grep("_1", colnames(brain0_filt_mat)),1] <- rep(1,length(grep("_1", colnames(brain0_filt_mat))))
brain0_des[grep("_2", colnames(brain0_filt_mat)),2] <- rep(1,length(grep("_2", colnames(brain0_filt_mat))))
brain0_des[grep("_3", colnames(brain0_filt_mat)),3] <- rep(1,length(grep("_3", colnames(brain0_filt_mat))))
brain0_des[is.na(brain0_des)] <- 0
brain0_des <- as.matrix(brain0_des)
colnames(brain0_des) <- c("Y", "M", "O")

# Correlations and significance in each condition 
brain0_cor_res = getCors(inputMat = brain0_filt_mat, design = brain0_des, corrType = "spearman")
  
# Differential correlations between conditions
brain0_dcpairs_YO = pairwiseDCor(brain0_cor_res, compare = c("Y", "O"))
brain0_dcpairs_YM = pairwiseDCor(brain0_cor_res, compare = c("Y", "M"))
brain0_dcpairs_MO = pairwiseDCor(brain0_cor_res, compare = c("M", "O"))

# Extracting top differentially correlated pairs between conditions
brain0_dd_pairs_YO = na.omit(dcTopPairs(brain0_dcpairs_YO, nPairs = (((nrow(brain0_filt_mat)*(nrow(brain0_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain0_dd_pairs_YO <- brain0_dd_pairs_YO[brain0_dd_pairs_YO$Classes!="NonSig",] 
brain0_dd_pairs_YO_rel <- brain0_dd_pairs_YO[brain0_dd_pairs_YO$Classes!="+/+",] 
brain0_dd_pairs_YO_rel <- brain0_dd_pairs_YO_rel[brain0_dd_pairs_YO_rel$Classes!="-/-",] 
brain0_dd_pairs_YO_rel <- brain0_dd_pairs_YO_rel[brain0_dd_pairs_YO_rel$Classes!="0/0",]

brain0_dd_pairs_MO = na.omit(dcTopPairs(brain0_dcpairs_MO, nPairs = (((nrow(brain0_filt_mat)*(nrow(brain0_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain0_dd_pairs_MO <- brain0_dd_pairs_MO[brain0_dd_pairs_MO$Classes!="NonSig",] 
brain0_dd_pairs_MO_rel <- brain0_dd_pairs_MO[brain0_dd_pairs_MO$Classes!="+/+",] 
brain0_dd_pairs_MO_rel <- brain0_dd_pairs_MO_rel[brain0_dd_pairs_MO_rel$Classes!="-/-",] 
brain0_dd_pairs_MO_rel <- brain0_dd_pairs_MO_rel[brain0_dd_pairs_MO_rel$Classes!="0/0",]

brain0_dd_pairs_YM = na.omit(dcTopPairs(brain0_dcpairs_YM, nPairs = (((nrow(brain0_filt_mat)*(nrow(brain0_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain0_dd_pairs_YM <- brain0_dd_pairs_YM[brain0_dd_pairs_YM$Classes!="NonSig",] 
brain0_dd_pairs_YM_rel <- brain0_dd_pairs_YM[brain0_dd_pairs_YM$Classes!="+/+",] 
brain0_dd_pairs_YM_rel <- brain0_dd_pairs_YM_rel[brain0_dd_pairs_YM_rel$Classes!="-/-",] 
brain0_dd_pairs_YM_rel <- brain0_dd_pairs_YM_rel[brain0_dd_pairs_YM_rel$Classes!="0/0",]


##########  BRAIN1  ##############

brain1_young_edat <- cbind(brain1_edat[,which(colnames(brain1_edat)%in%rownames(brain1_splitbyage$`20-29`))], brain1_edat[,which(colnames(brain1_edat)%in%rownames(brain1_splitbyage$`30-39`))])
brain1_mid_edat <- cbind(brain1_edat[,which(colnames(brain1_edat)%in%rownames(brain1_splitbyage$`40-49`))], brain1_edat[,which(colnames(brain1_edat)%in%rownames(brain1_splitbyage$`50-59`))])
brain1_old_edat <- cbind(brain1_edat[,which(colnames(brain1_edat)%in%rownames(brain1_splitbyage$`60-69`))], brain1_edat[,which(colnames(brain1_edat)%in%rownames(brain1_splitbyage$`70-79`))])
brain1 <- cbind(brain1_young_edat, brain1_mid_edat, brain1_old_edat)

# Filter expression data for low expression and low dispersion (3030x350)
brain1_filt = DGCA::filterGenes(brain1, filterTypes = "central", "dispersion", sequential = T, filterCentralType = "median", filterCentralPercentile = 0.90, filterDispersionType = "cv", filterDispersionPercentile = 0.90)
rownames(brain1_filt) <- unlist(lapply(strsplit(rownames(brain1_filt), split = "[.]"), function(x){x[1]}))
brain1_filt <- cbind(rownames(brain1_filt), brain1_filt)
brain1_genes <- grex(cleanid(rownames(brain1_filt)))
brain1_filt <- as.data.frame(brain1_filt)
colnames(brain1_filt)[1] <- "ensembl_id"
brain1_genesym <- dplyr::left_join(brain1_filt, brain1_genes, by="ensembl_id")
brain1_genesym <- brain1_genesym[complete.cases(brain1_genesym),]
brain1_genesym <- brain1_genesym[-c(2045,2046),] # Remove duplicated rows with NPIPA5
rownames(brain1_genesym) <- brain1_genesym$hgnc_symbol
brain1_filt_mat <- brain1_genesym[,-c(1,(ncol(brain1_genesym)-5):ncol(brain1_genesym))]
brain1_filt_mat <- data.matrix(brain1_filt_mat)

# Characterize subjects into age groups
colnames(brain1_filt_mat)[1:16] <- paste(colnames(brain1_filt_mat)[1:16], "1", sep = "_")
colnames(brain1_filt_mat)[17:120] <- paste(colnames(brain1_filt_mat)[17:120], "2", sep = "_")
colnames(brain1_filt_mat)[120:230] <- paste(colnames(brain1_filt_mat)[120:230], "3", sep = "_")
brain1_groups <- factor(unlist(lapply(strsplit(colnames(brain1_filt_mat), split = "_"), function(x){x[2]})))

brain1_groups <- gsub(pattern = "1", replacement = "Y", x = brain1_groups)
brain1_groups <- gsub(pattern = "2", replacement = "M", x = brain1_groups)
brain1_groups <- gsub(pattern = "3", replacement = "O", x = brain1_groups)

# Make a design matrix
brain1_des <- matrix(nrow = ncol(brain1_filt_mat), ncol = 3)
brain1_des <- as.data.frame(brain1_des)
brain1_des[grep("_1", colnames(brain1_filt_mat)),1] <- rep(1,length(grep("_1", colnames(brain1_filt_mat))))
brain1_des[grep("_2", colnames(brain1_filt_mat)),2] <- rep(1,length(grep("_2", colnames(brain1_filt_mat))))
brain1_des[grep("_3", colnames(brain1_filt_mat)),3] <- rep(1,length(grep("_3", colnames(brain1_filt_mat))))
brain1_des[is.na(brain1_des)] <- 0
brain1_des <- as.matrix(brain1_des)
colnames(brain1_des) <- c("Y", "M", "O")

# Correlations and significance in each condition 
brain1_cor_res = getCors(inputMat = brain1_filt_mat, design = brain1_des, corrType = "spearman")

# Differential correlations between conditions
brain1_dcpairs_YO = pairwiseDCor(brain1_cor_res, compare = c("Y", "O"))
brain1_dcpairs_YM = pairwiseDCor(brain1_cor_res, compare = c("Y", "M"))
brain1_dcpairs_MO = pairwiseDCor(brain1_cor_res, compare = c("M", "O"))

# Extracting top differentially correlated pairs between conditions
brain1_dd_pairs_YO = na.omit(dcTopPairs(brain1_dcpairs_YO, nPairs = (((nrow(brain1_filt_mat)*(nrow(brain1_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain1_dd_pairs_YO <- brain1_dd_pairs_YO[brain1_dd_pairs_YO$Classes!="NonSig",] 
brain1_dd_pairs_YO_rel <- brain1_dd_pairs_YO[brain1_dd_pairs_YO$Classes!="+/+",] 
brain1_dd_pairs_YO_rel <- brain1_dd_pairs_YO_rel[brain1_dd_pairs_YO_rel$Classes!="-/-",] 
brain1_dd_pairs_YO_rel <- brain1_dd_pairs_YO_rel[brain1_dd_pairs_YO_rel$Classes!="0/0",]

brain1_dd_pairs_MO = na.omit(dcTopPairs(brain1_dcpairs_MO, nPairs = (((nrow(brain1_filt_mat)*(nrow(brain1_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain1_dd_pairs_MO <- brain1_dd_pairs_MO[brain1_dd_pairs_MO$Classes!="NonSig",] 
brain1_dd_pairs_MO_rel <- brain1_dd_pairs_MO[brain1_dd_pairs_MO$Classes!="+/+",] 
brain1_dd_pairs_MO_rel <- brain1_dd_pairs_MO_rel[brain1_dd_pairs_MO_rel$Classes!="-/-",] 
brain1_dd_pairs_MO_rel <- brain1_dd_pairs_MO_rel[brain1_dd_pairs_MO_rel$Classes!="0/0",]

brain1_dd_pairs_YM = na.omit(dcTopPairs(brain1_dcpairs_YM, nPairs = (((nrow(brain1_filt_mat)*(nrow(brain1_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain1_dd_pairs_YM <- brain1_dd_pairs_YM[brain1_dd_pairs_YM$Classes!="NonSig",] 
brain1_dd_pairs_YM_rel <- brain1_dd_pairs_YM[brain1_dd_pairs_YM$Classes!="+/+",] 
brain1_dd_pairs_YM_rel <- brain1_dd_pairs_YM_rel[brain1_dd_pairs_YM_rel$Classes!="-/-",] 
brain1_dd_pairs_YM_rel <- brain1_dd_pairs_YM_rel[brain1_dd_pairs_YM_rel$Classes!="0/0",]


##########  BRAIN2  ##############

brain2_young_edat <- cbind(brain2_edat[,which(colnames(brain2_edat)%in%rownames(brain2_splitbyage$`20-29`))], brain2_edat[,which(colnames(brain2_edat)%in%rownames(brain2_splitbyage$`30-39`))])
brain2_mid_edat <- cbind(brain2_edat[,which(colnames(brain2_edat)%in%rownames(brain2_splitbyage$`40-49`))], brain2_edat[,which(colnames(brain2_edat)%in%rownames(brain2_splitbyage$`50-59`))])
brain2_old_edat <- cbind(brain2_edat[,which(colnames(brain2_edat)%in%rownames(brain2_splitbyage$`60-69`))], brain2_edat[,which(colnames(brain2_edat)%in%rownames(brain2_splitbyage$`70-79`))])
brain2 <- cbind(brain2_young_edat, brain2_mid_edat, brain2_old_edat)

# Filter expression data for low expression and low dispersion (3030x350)
brain2_filt = DGCA::filterGenes(brain2, filterTypes = "central", "dispersion", sequential = T, filterCentralType = "median", filterCentralPercentile = 0.90, filterDispersionType = "cv", filterDispersionPercentile = 0.90)
rownames(brain2_filt) <- unlist(lapply(strsplit(rownames(brain2_filt), split = "[.]"), function(x){x[1]}))
brain2_filt <- cbind(rownames(brain2_filt), brain2_filt)
brain2_genes <- grex(cleanid(rownames(brain2_filt)))
brain2_filt <- as.data.frame(brain2_filt)
colnames(brain2_filt)[1] <- "ensembl_id"
brain2_genesym <- dplyr::left_join(brain2_filt, brain2_genes, by="ensembl_id")
brain2_genesym <- brain2_genesym[complete.cases(brain2_genesym),]
rownames(brain2_genesym) <- brain2_genesym$hgnc_symbol
brain2_filt_mat <- brain2_genesym[,-c(1,(ncol(brain2_genesym)-5):ncol(brain2_genesym))]
brain2_filt_mat <- data.matrix(brain2_filt_mat)

# Characterize subjects into age groups
colnames(brain2_filt_mat)[1:15] <- paste(colnames(brain2_filt_mat)[1:15], "1", sep = "_")
colnames(brain2_filt_mat)[16:162] <- paste(colnames(brain2_filt_mat)[16:162], "2", sep = "_")
colnames(brain2_filt_mat)[163:327] <- paste(colnames(brain2_filt_mat)[163:327], "3", sep = "_")
brain2_groups <- factor(unlist(lapply(strsplit(colnames(brain2_filt_mat), split = "_"), function(x){x[2]})))

brain2_groups <- gsub(pattern = "1", replacement = "Y", x = brain2_groups)
brain2_groups <- gsub(pattern = "2", replacement = "M", x = brain2_groups)
brain2_groups <- gsub(pattern = "3", replacement = "O", x = brain2_groups)

# Make a design matrix
brain2_des <- matrix(nrow = ncol(brain2_filt_mat), ncol = 3)
brain2_des <- as.data.frame(brain2_des)
brain2_des[grep("_1", colnames(brain2_filt_mat)),1] <- rep(1,length(grep("_1", colnames(brain2_filt_mat))))
brain2_des[grep("_2", colnames(brain2_filt_mat)),2] <- rep(1,length(grep("_2", colnames(brain2_filt_mat))))
brain2_des[grep("_3", colnames(brain2_filt_mat)),3] <- rep(1,length(grep("_3", colnames(brain2_filt_mat))))
brain2_des[is.na(brain2_des)] <- 0
brain2_des <- as.matrix(brain2_des)
colnames(brain2_des) <- c("Y", "M", "O")

# Correlations and significance in each condition 
brain2_cor_res = getCors(inputMat = brain2_filt_mat, design = brain2_des, corrType = "spearman")

# Differential correlations between conditions
brain2_dcpairs_YO = pairwiseDCor(brain2_cor_res, compare = c("Y", "O"))
brain2_dcpairs_YM = pairwiseDCor(brain2_cor_res, compare = c("Y", "M"))
brain2_dcpairs_MO = pairwiseDCor(brain2_cor_res, compare = c("M", "O"))

# Extracting top differentially correlated pairs between conditions
brain2_dd_pairs_YO = na.omit(dcTopPairs(brain2_dcpairs_YO, nPairs = (((nrow(brain2_filt_mat)*(nrow(brain2_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain2_dd_pairs_YO <- brain2_dd_pairs_YO[brain2_dd_pairs_YO$Classes!="NonSig",] 
brain2_dd_pairs_YO_rel <- brain2_dd_pairs_YO[brain2_dd_pairs_YO$Classes!="+/+",] 
brain2_dd_pairs_YO_rel <- brain2_dd_pairs_YO_rel[brain2_dd_pairs_YO_rel$Classes!="-/-",] 
brain2_dd_pairs_YO_rel <- brain2_dd_pairs_YO_rel[brain2_dd_pairs_YO_rel$Classes!="0/0",]

brain2_dd_pairs_MO = na.omit(dcTopPairs(brain2_dcpairs_MO, nPairs = (((nrow(brain2_filt_mat)*(nrow(brain2_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain2_dd_pairs_MO <- brain2_dd_pairs_MO[brain2_dd_pairs_MO$Classes!="NonSig",] 
brain2_dd_pairs_MO_rel <- brain2_dd_pairs_MO[brain2_dd_pairs_MO$Classes!="+/+",] 
brain2_dd_pairs_MO_rel <- brain2_dd_pairs_MO_rel[brain2_dd_pairs_MO_rel$Classes!="-/-",] 
brain2_dd_pairs_MO_rel <- brain2_dd_pairs_MO_rel[brain2_dd_pairs_MO_rel$Classes!="0/0",]

brain2_dd_pairs_YM = na.omit(dcTopPairs(brain2_dcpairs_YM, nPairs = (((nrow(brain2_filt_mat)*(nrow(brain2_filt_mat)-1)))/2), classify = TRUE, adjust = "BH", sigThresh = 0.05, corSigThresh = 0.05))
brain2_dd_pairs_YM <- brain2_dd_pairs_YM[brain2_dd_pairs_YM$Classes!="NonSig",] 
brain2_dd_pairs_YM_rel <- brain2_dd_pairs_YM[brain2_dd_pairs_YM$Classes!="+/+",] 
brain2_dd_pairs_YM_rel <- brain2_dd_pairs_YM_rel[brain2_dd_pairs_YM_rel$Classes!="-/-",] 
brain2_dd_pairs_YM_rel <- brain2_dd_pairs_YM_rel[brain2_dd_pairs_YM_rel$Classes!="0/0",]


brain0_dd_pairs_YO_pairs <- as.character(unlist(paste(brain0_dd_pairs_YO$Gene1, brain0_dd_pairs_YO$Gene2, sep = "_")))
brain0_dd_pairs_YO_pathway <- calculatePathwayEnrichment(pathways, intersect(brain0_dd_pairs_YO_pairs, unlist(pathways)), universeNum = 77754, sigCutoff = 0.1)
brain0_YO_pathway <- brain0_dd_pairs_YO_pathway$pathway

brain0_dd_pairs_YM_pairs <- as.character(unlist(paste(brain0_dd_pairs_YM$Gene1, brain0_dd_pairs_YM$Gene2, sep = "_")))
brain0_dd_pairs_YM_pathway <- calculatePathwayEnrichment(pathways, intersect(brain0_dd_pairs_YM_pairs, unlist(pathways)), universeNum = 77754, sigCutoff = 0.1)
brain0_YM_pathway <- brain0_dd_pairs_YM_pathway$pathway

brain0_dd_pairs_MO_pairs <- as.character(unlist(paste(brain0_dd_pairs_MO$Gene1, brain0_dd_pairs_MO$Gene2, sep = "_")))
brain0_dd_pairs_MO_pathway <- calculatePathwayEnrichment(pathways, intersect(brain0_dd_pairs_MO_pairs, unlist(pathways)), universeNum = 77754, sigCutoff = 0.1)
brain0_MO_pathway <- brain0_dd_pairs_MO_pathway$pathway

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
