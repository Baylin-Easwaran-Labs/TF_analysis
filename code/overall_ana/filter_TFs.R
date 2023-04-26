#Select top Transcription Factors (TFs) from the list using adjusted p value filter 
###################################################################################


## specify your working directory
wdir = "relative_path_to_code"
## specify which cis regulatory regions the analysis corresponds to: options("PLS", "pELS", "dELS")
cre.type = "dELS"  
setwd(paste(wdir, "overall_ana", cre.type, sep = "/"))
      
#load the enriched motif data from elmer
hyper.path = "hyper/getTF.hyper.TFs.with.motif.pvalue.rda"
hypo.path =  "hypo/getTF.hypo.TFs.with.motif.pvalue.rda"
load(hyper.path)
meth.cor.hyper = TF.meth.cor
dim(meth.cor.hyper)
remove(TF.meth.cor)
load(hypo.path)
meth.cor.hypo = TF.meth.cor
dim(meth.cor.hypo)
meth.cor.hyper[1:5, 1:4]
meth.cor.hypo[1:5, 1:4]
# identical(meth.cor.hyper, meth.cor.hypo)

## Function to determine the locations of n smallest elements of a numeric vector.
which.minn <- function(x,n=1){
  if (n==1)
    which.min(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=FALSE)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}

## Get the top TFs, for each TF carry out a p.adj across each enriched motifs (for each column)
## ..................
get.top.tx = function(cor.matrix){
  n.col = ncol(cor.matrix)
  #initialize a dataframe by pre-allocating space (quickest way to fill a dataframe)
  #top.tx <- data.frame(txf = character(n.row), minPadj = numeric(n.row), stringsAsFactors = FALSE)
  top.tf = character()
  for(col in 1:n.col){
    col.adjusted = p.adjust(cor.matrix[, col], method = "BH") ## fdr = BH
    col.adjusted.filter = col.adjusted[col.adjusted < 0.05]
    col.adjusted.tfs = rownames(cor.matrix)[col.adjusted < 0.05]
    ## get the tfs with minimum fdr
    col.adjusted.tfs.final = col.adjusted.tfs[which.minn(col.adjusted.filter, 100)]
    top.tf = c(top.tf, col.adjusted.tfs.final)
  }
  return(top.tf)
}

top.tx.hyper = unique(get.top.tx(meth.cor.hyper))
top.tx.hypo = unique(get.top.tx(meth.cor.hypo))
length(top.tx.hyper); length(top.tx.hypo)
intersect(top.tx.hyper, top.tx.hypo)

# tf.hyper.dELS = top.tx.hyper
# tf.hypo.dELS = top.tx.hypo

## Save the Filtered TFs  
save(top.tx.hyper, top.tx.hypo,
     file = paste0(wdir, "/overall_ana/", cre.type, "_topTFs.Rda"))
