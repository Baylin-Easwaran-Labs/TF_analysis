## This script present a workflow for identifying the TFs using the ELMER method
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
## specify which cis regulatory regions the analysis corresponds to: options("PLS", "pELS", "dELS")
cre.type = "dELS"  


## Load ELMER MAE object
load(paste0(wdir,"/mae_allprobes.Rda"))

#### load the differentially methylated cCRE probes: options: PLS, pELS, dELS #####
##########################################################

## load hypermethylated probes
load(paste0(wdir, "/diff_meth/unmethNormal_dELS_Probes.Rda"))
hyper.probes = as.character(na.omit(diffMethProbesAllComparisons$Tumor.vs.Normal$HyperMethProbes))
## load hypomethylated probes
load(paste0(wdir, "/diff_meth/methNormal_dELS_Probes.Rda"))
hypo.probes = as.character(na.omit(diffMethProbesAllComparisons$Tumor.vs.Normal$HypoMethProbes))
length(hyper.probes); length(hypo.probes)

meth.tcga = assays(mae)[[1]]
dim(meth.tcga)
common.hyper = unique(intersect(hyper.probes, rownames(meth.tcga)))  
common.hypo = unique(intersect(hypo.probes, rownames(meth.tcga)))  
length(common.hyper); length(common.hypo)  

pheno.elmer = colData(mae)
#pheno.elmer$TN

group.col <- "TN"
group1 <-  "Tumor"
group2 <- "Normal"
probes.meth = list("hyper" = common.hyper, "hypo" = common.hypo)

## change to working directory if you have wandered
dir.create(paste(wdir, "overall_ana", cre.type, sep = "/"))
setwd(paste(wdir, "overall_ana", cre.type, sep = "/"))


for (meth.status in list("hyper", "hypo")) {
# print(probes.meth[[meth.status]])
print(paste0("Finding enriched motifs for ", meth.status))
#### Identify enriched motifs ####
enriched.motif <- get.enriched.motif(data = mae,
                                     min.motif.quality = "DS",
                                     probes = probes.meth[[meth.status]], 
                                     dir.out = meth.status, 
                                     label = meth.status,
                                     min.incidence = 10,
                                     lower.OR = 1.1)
#-------------------------------------------------------------
# Identifying regulatory TFs                        |
#-------------------------------------------------------------
print(paste0("Finding TFs for ", meth.status))
TF <- get.TFs(data           = mae,
              enriched.motif = enriched.motif,
              dir.out        = meth.status,
              mode           = "unsupervised",
              group.col      = group.col,
              group1         = group1,
              group2         = group2,
              diff.dir       = meth.status,
              cores          = 2,
              label          = meth.status)

}








