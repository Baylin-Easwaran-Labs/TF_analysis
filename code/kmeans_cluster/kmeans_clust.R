## Carry out K-means clustering to the CRC expression datasets (TCGA + GSE39582)
## And plot the PCA components along with the K-means, CIMP and CMS classes 

#### Load Libraries  ###############
# install.packages("factoextra")
library(factoextra)
library(ggfortify)

## specify working directory ####
wdir = "relative_path_to_code"

#### Load the expression data related to CRC (TCGA + GSE39582) 
#### annotated with CIMP and CMS subtypes, and methylation associated TFs
load(paste(wdir, "data/exp_TCGA_GSE39582_cimpTFs.Rda", sep = "/"))
dim(exp.df)

####  Carry out Kmeans clustering  ##############
# set.seed(500)
## n.features is the number of CIMP specific cCRE-associated TFs
n.features = 476
print(n.features)
cl = kmeans(exp.df[1:n.features], centers = 3, iter.max = 10,  nstart = 10)
# cl$cluster
exp.df$KmeanClust = as.factor(cl$cluster)

dim(exp.df)
tail(exp.df)
res.pca <- prcomp(exp.df[1:n.features], scale = F)
fviz_eig(res.pca) ## needs factoextra lib

## Results for individuals
ind.contrib = as.data.frame(res.pca$x)
plot(ind.contrib$PC1, ind.contrib$PC2)
dim(ind.contrib); dim(exp.df)

####   Plot the clustering results  ############
# Graph settings
require(ggplot2)
require(RColorBrewer)
require(reshape)
require(plyr)

cols <- brewer.pal(8,"Set1")
cols1 <- c(cols[2], cols[1], cols[3:8]) # inverts red/blue order
# cols2 <- c("#33CC00", "#0072B2","chocolate1")
cols2 <- c("#0072B2","chocolate1","#33CC00")
col.cimp = c("cyan","#33CC00", "chocolate1", "gray")
col.cms =  c("#3399CC","#009900", "#99CCFF", "#CC9900")

# Convex hulls.
pc.df = ind.contrib; dim(pc.df)
identical(rownames(pc.df), rownames(exp.df))
pc.df$kmean = exp.df$KmeanClust
pc.df$cimp = exp.df$cimp
pc.df$cms = exp.df$cms_label
find_hull <- function(issp06) issp06[chull(issp06$PC2, issp06$PC1), ]
hulls <- ddply(pc.df, "kmean", find_hull)
dim(hulls)

fig1 <- ggplot(data=pc.df, aes(PC1, PC2, colour=kmean, fill=kmean)) + 
  geom_point() + 
  geom_hline(yintercept=0, colour="darkgrey") + 
  geom_vline(xintercept=0, colour="darkgrey") +
  scale_colour_manual("kmeans", values = cols2) +
  scale_fill_manual("kmeans", values = cols2) +
  labs(x = "PC1", y = "PC2") +
  theme_classic(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")); 
fig1b <- fig1 + geom_polygon(data=hulls, alpha=.1, color = "grey"); fig1b

fig2 <- ggplot(data=pc.df, aes(PC1, PC2, colour=cimp, fill=kmean)) + 
  geom_point(size = 2.5) + 
  geom_hline(yintercept=0, colour="darkgrey") + 
  geom_vline(xintercept=0, colour="darkgrey") +
  scale_fill_manual("kmeans", values = cols2) +
  scale_colour_manual("CIMP", values = col.cimp) +
  labs(x = "PC1", y = "PC2") +
  theme_classic(base_size = 30) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")); 
fig2b <- fig2 + geom_polygon(data=hulls, alpha=.1, color = "grey"); fig2b

fig3 <- ggplot(data=pc.df, aes(PC1, PC2, colour=cms, fill=kmean)) + 
  geom_point(size = 2.5) + 
  geom_hline(yintercept=0, colour="darkgrey") + 
  geom_vline(xintercept=0, colour="darkgrey") +
  scale_colour_manual("CMS", values = col.cms) +
  scale_fill_manual("kmeans", values = cols2) +
  labs(x = "PC1", y = "PC2")+
  theme_classic(base_size = 30) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"));
# Convex hulls
fig3b <- fig3 + geom_polygon(data=hulls, alpha=.1, color = "grey"); fig3b

#### Save Plots  #######
ggsave(fig1b,width = 7, height = 6, filename=paste0(wdir, "/kmeans_cluster/pca_kmeans.pdf"))
ggsave(fig2b,width = 9, height = 6,filename=paste0(wdir, "/kmeans_cluster/pca_cimp.pdf"))
ggsave(fig3b,width = 9, height = 6,filename=paste0(wdir, "/kmeans_cluster/pca_cms.pdf"))



