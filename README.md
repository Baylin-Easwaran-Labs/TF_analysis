# TF_analysis
This folder contains scripts to generate the transcription factors and K-means clustering as presented in the manuscript. 
TCGA CRC methylation and expression data is downloaded using the commands in the script itself. 
Additional data needed for running the script including CIMP labels, autosomal probes, ENCODE cCRE probes, 
and normalized expression data for the combined CRC datasets (TCGA + GSE39582) for clustering analysis are provided 
as supplementary data to the manuscript.

Below are two major analysis pipelines utilized in our study.

A.	Transcription Factor (TF) Analysis using ELMER V2:

  1)	Refer to the directory “overall_ana” for identifying TFs related to overall analysis of TCGA CRC vs normal samples. 
      Refer to the directory “cimp_ana” for identifying TFs related to CIMP specific analysis of TCGA CRC samples
  2)	You will need to sequentially run three scripts for the TF calculation. Change to either the “overall_ana” 
      directory or the “cimp_ana” directory.
      
    a)	Run the script “get_Data.R” in order to download and prepare the methylation and expression data
    b)	Determine differentially methylated probes using the script “limma-difMeth.R”
    c)	Identify TFs using the script “get_TFs_ELMERv2.R”
    d)	Filter the TFs based on the adjusted p-value using the script “filter_TFs.R”
 
B.	Carry out K-means clustering on the expression data:

  Run the script “kmeans_clust.R”. This will carry out the K-means clustering as well as generate a PCA plot of the gene expression, 
  along with the K-means cluster, CIMP and CMS labels.
  
