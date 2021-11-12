###### Codes for estimating the number of clusters in scRNA-seq data
rm(list = ls())

## R codes for computing consensus matrix
library(Seurat)
library(dplyr)

setwd("/Users/suoqinjin/Google Drive/projects/project_Chen/e135_higher_depth/brachial_LMC")
# load Seurat object
load("e135_brachial_LMC_cLMCl.RData")
w10x <- w10x_bLMC_cLMCl

N <- ncol(w10x@data)
C <- matrix(0, N, N)
resR <- seq(0.1,3,by = 0.05) # a range of resolution for clustering
for (res in resR) {
  w10x <- FindClusters(object = w10x, reduction.type = "pca", dims.use = 1:numPC, resolution = res, algorithm = 1,save.SNN = TRUE,print.output = 0,force.recalc = T)
  print(length(unique(w10x@ident)))
  clusIndex <- as.numeric(as.character(w10x@ident))
  adjMat <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
  C <- C + adjMat
}
CM <- C/length(resR)
# save the consensus matrix for running the Matlab code 'determineNumClusters.m'
write.table(CM,file = "consensusMatrix_brachial_LMC_cLMCl.txt",sep = '\t', row.names = F, col.names = F)

## Matlab codes for determining the number of clusters
# Please do this in Matlab
fileName = 'brachial_LMC_cLMCl'; # the file name of the figure to save
CM = importdata(['consensusMatrix_',fileName,'.txt']);
[numCluster, numCluster0, eigenvalues] = determineNumClusters(CM,fileName);
# Output:
 #  numCluster: Number of inferred clusters
 #  numCluster0: the minimum number of inferred clusters based on the number of zero eigenvalues
 #  eigenvalues: eigenvalues





