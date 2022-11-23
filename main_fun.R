## Load the following libraries
library(Seurat)
library(readxl)
library(ggplot2)
library(tidyverse)
library(plyr)
library(dittoColors)
library(clustree)

###############################################################
### Load utility functions
### Note: Make sure that this R script is in the same directory
###############################################################
source("codes/util_funcs.R")

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

####################################################################
## Function to generate barplots of XFPs/Clusters
##
## Parameters:
## obj- The processed matrix object  (for example, Seurat object)
## filename- name of the file containing the original data
## save- Boolean instructing whether to save the file
####################################################################
drawBarplot = function(obj,filename,save=T){
	dat = data.frame(t(as.matrix(obj@assays$RNA@data)),check.names=F)
	dat$Clusters = Idents(run1)
	dat_melt = reshape2::melt(dat)
	dat_melt$variable = factor(dat_melt$variable,levels=rev(names(table(dat_melt$variable))))
	if(save){
		pdf(paste0("BARPLOT_",filename,".pdf"),,width=12,height=12)
		print(ggbarplot(dat_melt,x="variable",y="value",fill="black",lab.size=3)+facet_grid(~Clusters)+coord_flip()+geom_hline(yintercept=0,color="red"))
		dev.off()
	}else
	{
		print(ggbarplot(dat_melt,x="variable",y="value",fill="black",lab.size=3)+facet_grid(~Clusters)+coord_flip()+geom_hline(yintercept=0,color="red"))
	}
}

###########################################################################################
### Function to visualize hyperspectral data
### Currently done using SEURAT application
###########################################################################################
## Explanation of the parameters
##
## file: file object stored as a data frame
## col_start: Column number where the data actually begins
## meta: name of the meta data column that's used to label the plot
## num_neighbors: Number of nearest neighbors used in local approximations of manifold structure (UMAP hyperparameter)
## min_dist: This controls how tightly the embedding is allowed compress points together (UMAP hyperparameter)
## cluster: Boolean (True/False) whether clustering needs to be done
## plots: Print plots
############################################################################################
run_primary_analysis = function(file,col_start=3,num_neighbors,min_dist,meta="orig.ident",cluster=F,plots=T){
	dat = prep_data(file_name=file)
	dat = dat[,col_start:ncol(dat)]
	dat = t(dat)
	out = split(c(410:696), sort(c(410:696)%%96))
	rownames(dat) = laply(1:length(out),function(i) paste(out[i],collapse="_"))
	datN = transformData(dat,cfLow=0, cfHigh=5)
	obj = CreateSeuratObject(counts=datN)
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
file="mTagBFP2_mTurq2_mUKG_mGreenLantern_mKOk_mScarlet_SuperNovaRed_mNeptune_mCardinal.xlsx"
run1 = run_primary_analysis(file,col_start=3,num_neighbors=50L,min_dist=0.01,plot=T)

## If clustering is necessary
run1 = FindNeighbors(run1, reduction = "pca", dims = 1:10,verbose=F) %>% FindClusters(resolution=c(0.4,0.6,0.8),random.seed=123,verbose=F)

## Cluster Stability plots
clustree(run1,prefix="RNA_snn_res.",layout="sugiyama")
clustree(run1, prefix = "RNA_snn_res.", node_colour = "sc3_stability")
run1 = BuildClusterTree(run1)
PlotClusterTree(run1,type="phylogram",use.edge.length=F,edge.color="black",edge.width = 2)

## Plot UMAP scattered plots
Idents(run1)<-'RNA_snn_res.0.6"'
DimPlot(run1,combine=T,cols=dittoColors(),label=T,label.size=5)+NoLegend()

## Draw barplots
drawBarplot(obj=run1,save=T,filename=gsub(".xlsx","",file))

## Plotting
# xfp_colors=c("#1F1850","#FF5900","#770B00","#095A98","#31C331","#FF1800","#FEB300","#250FB0","#298297","#438EB3","#29A37B","#FF8B00")
# DimPlot(run1,group.by="orig.ident",cols=xfp_colors)
# Idents(run1)<-'RNA_snn_res.0.6"'
# DimPlot(run1,combine=T,cols=dittoColors(),label=T,label.size=5)+NoLegend()
# FeaturePlot(run1, features = rownames(run1)[c(1,10,30,40)], min.cutoff = "q05", max.cutoff = "q95", cols = c("lightblue", "magenta"), pt.size = 1,order=T,label=F,ncol=2)
# DoHeatmap(run1,features=rownames(run1),size=3,group.by="orig.ident",raster=F,group.colors =xfp_colors)
# ALLmarkers = FindAllMarkers(object = run1,only.pos = TRUE, min.pct = 0.15,logfc.threshold=0.01)
# ALLmarkers = ALLmarkers[ALLmarkers$p_val_adj<=0.05,]
# ALLmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(run1,features=top10$gene,size=3,raster=F)
