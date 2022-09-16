library(Seurat)
library(readxl)
library(ggplot2)
library(tidyverse)
library(plyr)

#######################################################
## If the file is formatted in Excel
## I hope the first column in the excel file is the true
##  					XFP label and the second col is ObjectName
########################################################
prep_data = function(file_name = file){
	f = data.frame(read_excel(file_name),check.names=F)
	rownames(f) = paste(f$XFP,f$ObjectNumber,sep="_")
	return(f)
}

#############################################
### Function to visualize hyperspectral data
### Currently done using SEURAT application
#############################################
## Explanation of the parameters
##
## file: file object stored as a data frame
## col_start: Column number where the data actually begins
## meta: name of the meta data column that's used to label the plot
## num_neighbors: Number of nearest neighbors used in local approximations of manifold structure (UMAP hyperparameter)
## min_dist: This controls how tightly the embedding is allowed compress points together (UMAP hyperparameter)
## cluster: Boolean (True/False) whether clustering needs to be done
## plots: Print plots
##############################################
run_primary_analysis = function(dat,col_start=3,num_neighbors,min_dist,meta="orig.ident",cluster=F,plots=T){
	dat = dat[,col_start:ncol(dat)]
	dat = t(dat)
	out = split(c(410:696), sort(c(410:696)%%96)) 
	rownames(dat) = laply(1:length(out),function(i) paste(out[i],collapse="_"))
	obj = CreateSeuratObject(counts=dat)
	VariableFeatures(obj) = rownames(obj)
	obj = ScaleData(obj,verbose=F)
	obj = RunPCA(obj,verbose=F,npcs=20)
	obj = RunUMAP(obj,dims=1:10,umap.method="uwot",n.neighbors=num_neighbors,min.dist=min_dist,metric="cosine",verbose=F)
	if(plots) print(DimPlot(obj,group.by=meta)+labs(title=""))
	if(cluster){
		obj = FindNeighbors(obj, reduction = "pca", dims = 1:10,verbose=F)
		obj = FindClusters(obj,resolution=c(0.4,0.6,0.8),random.seed=123,verbose=F)
		print(DimPlot(obj)+labs(title=""))
	}
	return(obj)
}

## Run the analysis
file="96Bins_Mostly512x512_MeanObjectIntensity.xlsx"

f = prep_data(file_name=file)
run1 = run_primary_analysis(dat=f,col_start=3,num_neighbors=80L,min_dist=0.25)
head(run1)

## If clustering is necessary
run1 = FindNeighbors(run1, reduction = "pca", dims = 1:10,verbose=F) %>% FindClusters(resolution=c(0.4,0.6,0.8),random.seed=123,verbose=F)
DimPlot(run1,group.by="RNA_snn_res.0.6")
