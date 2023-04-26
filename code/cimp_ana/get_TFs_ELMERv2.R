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

#### load the differentially methylated cCRE probes  #####
##########################################################

load(paste0(wdir, "/cimp_ana/diff_meth/unmethNormal_dELS_Probes.Rda"))
hyper.probes = as.character(na.omit(diffMethProbesAllComparisons$cimpH.vs.Normal$HyperMethProbes))
load(paste0(wdir, "/cimp_ana/diff_meth/methNormal_dELS_Probes.Rda"))
hypo.probes = as.character(na.omit(diffMethProbesAllComparisons$cimpH.vs.Normal$HypoMethProbes))
length(hyper.probes); length(hypo.probes)

meth.tcga = assays(mae)[[1]]
dim(meth.tcga)
common.hyper = unique(intersect(hyper.probes, rownames(meth.tcga)))  
common.hypo = unique(intersect(hypo.probes, rownames(meth.tcga)))  
length(common.hyper); length(common.hypo)  

pheno.elmer = colData(mae)
#pheno.elmer$TN

## Load CIMP classes  ##
load(paste0(wdir, "/data/cimpClassSamples_tcga_Elmer.Rda"))
pheno.elmer$cimp = ifelse(pheno.elmer$sample %in% cimpH.samples & pheno.elmer$TN == "Tumor" , "cimpH",
                          ifelse(pheno.elmer$sample %in% cimpInt.samples & pheno.elmer$TN == "Tumor", "cimpInt",
                          ifelse(pheno.elmer$sample %in% cimpL1.samples & pheno.elmer$TN == "Tumor", "cimpL1",
                          ifelse(pheno.elmer$sample %in% cimpL2.samples & pheno.elmer$TN == "Tumor", "cimpL2", "Normal"
                          ))))       
table(pheno.elmer$cimp)
## reassign the pheno.elmer with cimp class into the elmer MAE object
colData(mae) = pheno.elmer

group.col <- "cimp"
#### Specify which CIMP group you want to compare against the normal ####
group1 <-  "cimpH"
group2 <- "Normal"
probes.meth = list("hyper" = common.hyper, "hypo" = common.hypo)

## change to working directory if you have wandered
dir.create(paste(wdir, "cimp_ana", cre.type, sep = "/"))
setwd(paste(wdir, "cimp_ana", cre.type, sep = "/"))


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


####################### Comments on STAD, ESAC, GBM_LGG  #########################################

# Gastric cancer (STAD: stomach adenocarcinoma), esophagus cancer (ESAC: esophagus adenocarcinoma) 
# and glioma (GBM: glioblastoma multiforme and LGG: brain lower grade glioma) samples were 
# downloaded from GDC Data Portal using TCGAbiolink and processed as above for the TF analysis. 
# For gastric cancer we used 372 tumor samples and 30 normal samples belonging to all 
# gastrointestinal cancers. For esophagus cancer, we used 88 esophagus adenocarcinoma 
# tumor samples and 30 normal samples belonging to all gastrointestinal cancers. Esophagus 
# squamous cell carcinoma samples were not included for analysis. For gliomas, we used 102 
# tumor samples (GBM, LGG) and 400 solid tissue normal samples across all TCGA cancers as normal reference.  
# For esophagus and glioma cancer, the CIMP-H and L (intermediate methylation) subtypes were combined 
# into a high methylation group (CIMP-H’) and subtypes L3 and L4 were combined into a low 
# methylation group (CIMP-L’) due to low sample number availability.

#######################################################################################################


