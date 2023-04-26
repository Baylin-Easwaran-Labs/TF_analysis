## Download the expression and methylation data using getTCGA (ELMER)
## and Prepare the MAE data object
# "ELMER v.2: an R/Bioconductor package to reconstruct gene regulatory networks from 
# DNA methylation and transcriptome profiles"
# Bioinformatics 2019 Jun 1;35(11):1974-1977. doi: 10.1093/bioinformatics/bty902

## install and load ELMER v2
library(stringr)
library(TCGAbiolinks)
library(dplyr)
library(ELMER)
library(MultiAssayExperiment)

## specify your working directory
wdir = "relative_path_to_code"

####  Download the colon data from GDC data commons  ####
#########################################################

# specify the path where you want to download the data at the "basedir" variable)
getTCGA(disease = "COAD",
        basedir = wdir,
        Meth = TRUE, RNA = TRUE,
        genome = "hg19")

#### create the ELMER mae object with all probes  #####
#######################################################

# exp and met ate the expression and methylation data downloaded by the getTCGA subroutine
mae <- createMAE(exp = paste0(wdir,"/COAD/COAD_RNA_hg19.rda"), 
                 met = paste0(wdir,"/COAD/COAD_meth_hg19.rda"), 
                 met.platform = "450K",
                 genome = "hg19",
                 linearize.exp = TRUE,
                 filter.probes = NULL,
                 met.na.cut = 0.2,
                 save = FALSE,
                 TCGA = TRUE) 
# Remove FFPE samples from the analysis
mae <- mae[,!mae$is_ffpe]
save(mae, file = paste0(wdir,"/mae_allprobes.Rda"))


