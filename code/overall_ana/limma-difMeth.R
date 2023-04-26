#### differential Methylation Analysis of the Tumor Vs Normal sampels
#### Using limma package

#updating tibble
# Installing
#install.packages("tibble")
#install.packages("VennDiagram")
library(limma)
library("tibble")


#### load the Elmer MAE object with all probes  ####
####################################################

load("relative_path_to_code/mae_allprobes.Rda")
meth.tcga = assays(mae)[[1]]
# exp.tcga = assays(mae)[[2]]
dim(meth.tcga)
pheno.elmer = colData(mae)
table(pheno.elmer$TN)
normalTCGA.samples = rownames(pheno.elmer[pheno.elmer$TN == "Normal", ])
length(normalTCGA.samples)

## load the set of autosomal probes and filter for autosomal probes
load("relative_path_to_code/data/autosomal_450k_Probes.Rda")
meth.tcga.auto <- meth.tcga[rownames(meth.tcga) %in% autosomalProbes,] 
dim(meth.tcga); dim(meth.tcga.auto)
#fix the colnames starting with X
colnames(meth.tcga.auto) = substr(colnames(meth.tcga.auto),1,16)

## load cCRE probes identified from the ENCODE project
## Expanded encyclopaedias of DNA elements in the
## human and mouse genomes. ENCODE Project Consortium, Nature 583, 699–710 (2020).
load("relative_path_to_code/data/filtered_CRE_Probes.Rda")
length(dELS_CGI_probes) ; length(dELS_nonCGI_probes)
length(pELS_CGI_probes); length(pELS_nonCGI_probes)
length(PLS_CGI_probes); length(PLS_nonCGI_probes)

## select one (dELS, pELS, PLS probes)
cre.probes = unique(na.omit(c(dELS_CGI_probes, dELS_nonCGI_probes)))
length(cre.probes)

meth.tcga.sort <- meth.tcga.auto[cre.probes[cre.probes %in% rownames(meth.tcga.auto)],]
dim(meth.tcga.sort)
meth.tcga.sort.unmethNormal <- meth.tcga.sort[apply(meth.tcga.sort[ ,normalTCGA.samples], 1, mean) <= 0.2, ]
meth.tcga.sort.methNormal <- meth.tcga.sort[apply(meth.tcga.sort[ ,normalTCGA.samples], 1, mean) >= 0.6, ]
dim(meth.tcga.sort);  dim(meth.tcga.sort.unmethNormal); dim(meth.tcga.sort.methNormal) 

## 1.Idendify the hypermathylated probes using CpGs unmethylated in Normals
lapply(list(meth.tcga.sort.unmethNormal), FUN=function(x, deltaBeta=0.2)
{

design <- model.matrix(~0+factor(pheno.elmer$TN))
colnames(design) = levels(factor(pheno.elmer$TN))
contrast.matrix = makeContrasts(Tumor - Normal, levels =design)
## check the design matrix to see if samples are read correctly acc. to pheno.elmer$TN
head(design); dim(design); head(pheno.elmer$TN); head(contrast.matrix)
fit <- lmFit(x, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# All comparison labels:
allComparisons <-  gsub(" - ",".vs.",colnames(contrast.matrix))
# Create list to hold all diff methylated probes
diffMethProbesAllComparisons=vector(mode="list", length=length(allComparisons))
names(diffMethProbesAllComparisons) <- allComparisons

# Create list to hold all diff methylated genes
diffMethGenesAllComparisons=vector(mode="list", length=length(allComparisons))
names(diffMethGenesAllComparisons) <- allComparisons

# Create list to hold toptables
topTablesAllComparisons <- vector(mode="list", length=length(allComparisons))
names(topTablesAllComparisons) <- allComparisons
# Get top differentially methylated probes between for each comparison:
for(i in 1:length(allComparisons)){
  diffProbes <- topTable(fit2, coef=i, adjust="BH", number=nrow(x))
  hyperMethProbes <- (diffProbes[diffProbes$logFC >= deltaBeta & diffProbes$adj.P.Val <= 0.05, "ID"])
  hypoMethProbes <- (diffProbes[diffProbes$logFC <= -deltaBeta & diffProbes$adj.P.Val <= 0.05, "ID"])
  
  diffMethProbesAllComparisons[[i]]["HyperMethProbes"] <- list(unique(hyperMethProbes))
  diffMethProbesAllComparisons[[i]]["HypoMethProbes"] <- list(unique(hypoMethProbes))
  topTablesAllComparisons[[i]] <- diffProbes
}

# Decide tests
decideTestsResults <- decideTests(fit2)

# Save output from this function to the workspace (global environment)
{
  assign("diffMethProbesAllComparisons", diffMethProbesAllComparisons, envir=.GlobalEnv)
  assign("topTablesAllComparisons", topTablesAllComparisons, envir=.GlobalEnv)
  assign("decideTestsResults", decideTestsResults, envir=.GlobalEnv)
}  
gc()
}
)
save(diffMethProbesAllComparisons, topTablesAllComparisons, decideTestsResults,
     file = "relative_path_to_code/diff_meth/unmethNormal_dELS_Probes.Rda")

## 2.Idendify the hypomethylated probes using CpGs methylated in Normals
lapply(list(meth.tcga.sort.methNormal), FUN=function(x, deltaBeta=0.2)
{
  design <- model.matrix(~0+factor(pheno.elmer$TN))
  colnames(design) = levels(factor(pheno.elmer$TN))
  contrast.matrix = makeContrasts(Tumor - Normal, levels =design)
  ## check the design matrix to see if samples are read correctly acc. to pheno.elmer$TN
  head(design); dim(design); head(pheno.elmer$TN); head(contrast.matrix)
  fit <- lmFit(x, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # All comparison labels:
  allComparisons <-  gsub(" - ",".vs.",colnames(contrast.matrix))
  # Create list to hold all diff methylated probes
  diffMethProbesAllComparisons=vector(mode="list", length=length(allComparisons))
  names(diffMethProbesAllComparisons) <- allComparisons
  
  # Create list to hold all diff methylated genes
  diffMethGenesAllComparisons=vector(mode="list", length=length(allComparisons))
  names(diffMethGenesAllComparisons) <- allComparisons
  
  # Create list to hold toptables
  topTablesAllComparisons <- vector(mode="list", length=length(allComparisons))
  names(topTablesAllComparisons) <- allComparisons
  # Get top differentially methylated probes between for each comparison:
  for(i in 1:length(allComparisons)){
    diffProbes <- topTable(fit2, coef=i, adjust="BH", number=nrow(x))
    hyperMethProbes <- (diffProbes[diffProbes$logFC >= deltaBeta & diffProbes$adj.P.Val <= 0.05, "ID"])
    hypoMethProbes <- (diffProbes[diffProbes$logFC <= -deltaBeta & diffProbes$adj.P.Val <= 0.05, "ID"])

    diffMethProbesAllComparisons[[i]]["HyperMethProbes"] <- list(unique(hyperMethProbes))
    diffMethProbesAllComparisons[[i]]["HypoMethProbes"] <- list(unique(hypoMethProbes))
    topTablesAllComparisons[[i]] <- diffProbes
  }
  
  # Decide tests
  decideTestsResults <- decideTests(fit2)
  
  # Save output from this function to the workspace (global environment)
  {
    assign("diffMethProbesAllComparisons", diffMethProbesAllComparisons, envir=.GlobalEnv)
    assign("topTablesAllComparisons", topTablesAllComparisons, envir=.GlobalEnv)
    assign("decideTestsResults", decideTestsResults, envir=.GlobalEnv)
  }  
  gc()
}
)
save(diffMethProbesAllComparisons, topTablesAllComparisons, decideTestsResults,
     file = "relative_path_to_code/diff_meth/methNormal_dELS_Probes.Rda")



