#!/share/bin/Rscript
rm(list=ls())
options(warn=-1)
options(stringsAsFactors = FALSE)
################################################################################################################################################################
##
## Authors : Dr Stephen J Newhouse
## Email :   stephen..j.newhouse@gmail.com
## Date :    15-Nov-2013
## Title :   Microarry pre-processing for Illumina BeadArray data.
## Verison:  1.00
## Source:   Genomestudio
##
## Usage:    $ Rscrpit gx_illumina_pipeline_preProcessing.1.0.R <project.config>
##
## Essential Files/Data
##
## 1) Final reports from Genomestudio: a) Sample, b) Group Probe, c) Control Probe, d) probe annotations
## 2) Phenotype file:  must contain the following columns:- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID
## 3) Technical information file with batch processing and lab realted data. must contain the following columns:- Sample.ID followed by any and all processing data information, RIN, RNA yeilds and concetrations
##
#################################################################################################################################################################
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat("########################################################################################","\r","\n")
cat("## Authors : Dr Stephen J Newhouse","\r","\n")
cat("## Email :   stephen..j.newhouse@gmail.com","\r","\n")
cat("## Date :    15-Nov-2013","\r","\n")
cat("## Title :   Microarry pre-processing for Illumina BeadArray data.","\r","\n")
cat("## Verison:  1.00","\r","\n")
cat("## Source:   Genomestudio","\r","\n")
cat("## Usage:    $ Rscrpit gx_illumina_pipeline_preProcessing.1.0.R <project.config>","\r","\n")
cat("########################################################################################","\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" Loading Libraries ","\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)
library(gplots)
library(limma)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)
source("./pre_process_gx.R")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

###################
## some functins ##
##################
## subset eset
removeSamples_eset_lumi <- function(eset, sampleRemove) {
sample <- sampleNames(eset)
samples_to_remove <- sampleRemove
samples_to_keep <- (sample %in% samples_to_remove)==FALSE
sel_samp_names <- sampleNames(eset)[samples_to_keep]
eset <- eset[,samples_to_keep]
ControlData <- getControlData(eset)
ControlData  <- ControlData[,c("controlType","ProbeID",sel_samp_names)]
eset <- addControlData2lumi(ControlData , eset)
return(eset)
}
## shuffle_cols
shuffle_cols <- function(data) {
data_var <- names(data)
data_var_shuffled <- sample(data_var,size=length(data_var),replace=FALSE)
return(data_var_shuffled)
}
## shuffle_rows
shuffle_rows <- function(data) {
data_var <- rownames(data)
data_var_shuffled <- sample(data_var,size=length(data_var),replace=FALSE)
return(data_var_shuffled)
}
#####################################################################################################
############################
## negBeadOutlierRepMean
#############################
negBeadOutlierRepMean <- function(x) {
		z_out_samp <- abs(  as.numeric( scale(x) )  ) > 2
		mean_pop <- mean(x[z_out_samp==FALSE])
		sd_pop <- sd(x[z_out_samp==FALSE])
		new_x <- ifelse( abs(  as.numeric( scale(x) )  ) > 2, mean_pop, x )
		return(new_x)
}
##############
## quantfun ##
##############
quantfun <- function(x) {
as.integer(cut(x, quantile(x, probs=c(seq(0,1,0.20)  ) ), include.lowest=TRUE))
}
##############
## var_gene ##
##############
zero_var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
zero_var <- var_gx==0
return(zero_var)
}
mean_probe <- function(gx_matrix) {
gx <- gx_matrix
mean_gx <- apply(gx,1,mean)
return(mean_gx)
}
sd_probe <- function(gx_matrix) {
gx <- gx_matrix
sd_gx <- apply(gx,1,sd)
return(sd_gx)
}
var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
return(var_gx)
}
max_probe <- function(gx_matrix) {
gx <- gx_matrix
max_gx <- apply(gx,1,max)
return(max_gx)
}
min_probe <- function(gx_matrix) {
gx <- gx_matrix
min_gx <- apply(gx,1,min)
return(min_gx)
}
has_var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
has_var <- var_gx!=0
return(has_var)
}
has_var_probe2 <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
has_var <- var_gx >= quantile(var_gx, probs="0.1" );
return(has_var)
}

##############
## qc_plots ##
##############
gx_qc_plots_lumi <- function(eset, outfile ,do_pca=TRUE ) {
outpdf <- outfile;
##  eset <- eset_raw;
## get pheno data
pheno <- pData(eset)
## basic colours
chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
group_col <- labels2colors( as.character(pData(eset)$GROUPS))
pheno_col <- labels2colors( as.character(pData(eset)$PHENOTYPE))
gender_col <- labels2colors( as.character(pData(eset)$SEX))
tissue_col <- labels2colors( as.character(pData(eset)$TISSUE))
## batch pheno data
sel_tech <- grep("tech",names(pheno))
batch_pheno <- pheno[,sel_tech]
## id what is char/fav versus numerical
sel_batch <- sapply(batch_pheno ,class) %in% c("character","factor")
sel_num_batch <- sapply(batch_pheno ,class) %in% c("numeric")
## get names
batch_var_names <- names( batch_pheno[sel_batch])
batch_var_names_numeric <- names(batch_pheno[sel_num_batch])
## quantfun numeric data
numeric_tech <- as.data.frame(batch_pheno[,sel_num_batch])
colnames(numeric_tech) <- batch_var_names_numeric
quant_numeric <- apply( numeric_tech,2,quantfun)
## colours
batch_col <- apply(batch_pheno[,sel_batch],2,labels2colors) ## colours
quant_numeric_col <- apply(quant_numeric,2,numbers2colors)
datColors <- cbind(chip_col,group_col,pheno_col,gender_col,tissue_col,batch_col,quant_numeric_col)
## expression matrix and IAC
gx <- t(exprs(eset));
datExprs <- exprs(eset)
IAC <- cor(datExprs)
IAC_d <- 1-IAC
## flashClust
dist_exprs <- dist( t(datExprs), method="e" )
sampleTree <- flashClust( dist_exprs, method = "average");
## PLOTS
if(do_pca==TRUE) {
pca_raw <- prcomp(gx)$x;
pdf(outpdf, width=16.5,height=11.7)
## Standard plots
	plot(eset, what='boxplot', col=chip_col )
	plot(eset, what='density' )
	plot(eset, what='cv'  )
## PCA plots
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color="black", pch=21,bg=chip_col)
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Group ",color="black", pch=21,bg=group_col)
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Phenotype ",color="black", pch=21,bg=pheno_col)
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Gender ",color="black", pch=21,bg=gender_col)
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Gender ",color="black", pch=21,bg=tissue_col)
## loop PCA and plot by tech var
	for(tech_var in batch_var_names) {
		tech_var_name <- paste(tech_var,sep="")
		tech_var_col <- labels2colors( as.character(batch_pheno[,tech_var_name]) )
		scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main=paste("3D Scatterplot coloured by ",tech_var_name,sep=""),color="black",pch=21,bg=tech_var_col)
		}
			for(tech_var in batch_var_names_numeric) {
		tech_var_name <- paste(tech_var,sep="")
		tech_var_col <- numbers2colors( quant_numeric[,tech_var_name])
		scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main=paste("3D Scatterplot coloured by ",tech_var_name,sep=""),color="black",bg=tech_var_col,pch=21)
		}
## other plots
	tmp_iac <- outlierSamples(datExprs,thresh=iac_sd_thrs, showplots=TRUE)
	plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,dendroLabels = pData(eset)$Sample.ID, main="Sample dendrogram and trait heatmap")
	heatmap.2(IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )
dev.off()
} else {
pdf(outpdf, width=16.5,height=11.7)
	plot(eset, what='boxplot', col=chip_col )
	plot(eset, what='density' )
	plot(eset, what='cv'  )
	tmp_iac <- outlierSamples(datExprs,thresh=iac_sd_thrs, showplots=TRUE)
	plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,dendroLabels = pData(eset)$Sample.ID, main="Sample dendrogram and trait heatmap")
	heatmap.2(IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )

dev.off()
  }
}

############################
## write_expression_files ##
############################

write_expression_files <- function(eset, outfile) {

cat(" Writing probe exprs matrix [", paste(outfile,".exprs_matrix.txt",sep="")  ,"]","\r","\n")
gx <- exprs(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".exprs_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe se.exprs matrix [", paste(outfile,".se.exprs_matrix.txt",sep="")  ,"]","\r","\n")
gx <- se.exprs(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".se.exprs_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe detection matrix [", paste(outfile,".detection_matrix.txt",sep="")  ,"]","\r","\n")
gx <- detection(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".detection_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe beadNum matrix [", paste(outfile,".beadNum_matrix.txt",sep="")  ,"]","\r","\n")
gx <- beadNum(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".beadNum_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe PCA matrix [", paste(outfile,".pca_matrix.txt",sep="")  ,"]","\r","\n")
gx <- exprs(eset)
gx <- gx[zero_var_probe(gx)==FALSE,]
pca_gx <- prcomp(t(gx))$x
pca_gx <- pca_gx[,1:20]
pca_gx <- cbind(rownames(pca_gx),pca_gx)
grep_pc <- grep("PC",colnames(pca_gx))
colnames(pca_gx) <-   c("Sample.ID",colnames(pca_gx)[grep_pc])
write.table(pca_gx  , file=paste(outfile,".pca_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing pData slot of eset and adding PCA data to [", paste(outfile,".pData.txt",sep="")  ,"]","\r","\n")
pgx <- pData(eset)
pgx <- merge(pgx, pca_gx, by.x="Sample.ID",by.y="Sample.ID",sort=FALSE,all.x=TRUE)
write.table(pgx  , file=paste(outfile,".pData.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing fData slot of eset [", paste(outfile,".fData.txt",sep="")  ,"]","\r","\n")
fgx <- fData(eset)
sel_probe <- paste(fgx$nuID,sep="")
gx <- exprs(eset)
gx <- gx[sel_probe,]
fgx$mean_probe <- mean_probe(gx);
fgx$sd_probe <- sd_probe(gx);
fgx$var_probe <- var_probe(gx);
fgx$min_probe <- min_probe(gx);
fgx$max_probe <- max_probe(gx);
write.table(fgx  , file=paste(outfile,".fData.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
rm("gx","pca_gx","pgx","fgx")
}

###################
## sampleNetwork ##
###################
## eset <- eset_raw  ## testing
## col_by_chip <- 1
## enforce all then group then all on good samples

basic_sampleNetwork <- function(eset,col_by_chip=c("1","0"),groups=c("all","byGroup"),outfile,sd_thrs=2) {

if(groups=="all") {

mgroup <- "all"
n_groups <- length(mgroup)
res_names <- c("Group","mean_IAC","n_outliers","min_Z.K","rho","rho_pvalue","Z.K_outliers")
res <- matrix(ncol=length(res_names),nrow=n_groups, dimnames=list(mgroup,res_names))

cat(" doing group [",mgroup,"]","\r","\n")
gx <- exprs(eset)
sample_info <- pData(eset)
samle_names <- sampleNames(eset)
groups <- as.character(sample_info$GROUPS)
pheno <- as.character(sample_info$PHENOTYPE)
gpcolors <- labels2colors(groups)
## IAC, ADJACENCY & fundamentalNetworkConcepts, Z.K, Z.C, Z.MAR
cat(" Calculating fundamentalNetworkConcepts Metrics ","\r","\n")
IAC=cor(gx,method="p",use="p")
diag(IAC)=0
A.IAC=((1+IAC)/2)^2  ## ADJACENCY MATRIX
FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
K2=FNC$ScaledConnectivity
Z.K=(K2-mean(K2))/sd(K2)
Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
Z.MAR=(FNC$MAR-mean(FNC$MAR))/sd(FNC$MAR)
rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
## OUTLIERS
Z.K_outliers <- Z.K < -sd_thrs
Z.K_outliers <- names(Z.K_outliers[Z.K_outliers==TRUE])
n_outliers <- length(Z.K_outliers)
mean_IAC <- mean(IAC[upper.tri(IAC)])
min_Z.K <- min(Z.K)
cat(" Number of Z.K outliers [", n_outliers,"]","\r","\n")
cat(" mean_IAC [", mean_IAC,"]","\r","\n")
## Data frame of .SampleNetwork_Stats.txt
cat(" Making Data fram of fundamentalNetworkConcepts Metrics ","\r","\n")
FNC_DF <- as.data.frame(FNC)
FNC_DF$Z.K=(K2-mean(K2))/sd(K2)
FNC_DF$Z.C=(FNC_DF$ClusterCoef-mean(FNC_DF$ClusterCoef))/sd(FNC_DF$ClusterCoef)
FNC_DF$Z.MAR=(FNC_DF$MAR-mean(FNC_DF$MAR))/sd(FNC_DF$MAR)
FNC_DF$mean_IAC <- mean(IAC[upper.tri(IAC)])
FNC_DF$Mean_Connectivity <- mean(FNC$Connectivity)
FNC_DF$Mean_ScaledConnectivity <- mean(FNC$ScaledConnectivity)
FNC_DF$Mean_ClusterCoef <- mean(FNC$ClusterCoef)
FNC_DF$Mean_MAR <- mean(FNC$MAR)
FNC_DF$Decentralization <- 1-FNC_DF$Centralization
FNC_DF$Homogeneity <- 1-FNC_DF$Heterogeneity
FNC_DF$rho <- signif(cor.test(FNC_DF$Z.K,FNC_DF$Z.C,method="s")$estimate,2)
FNC_DF$rho_pvalue <- signif(cor.test(FNC_DF$Z.K,FNC_DF$Z.C,method="s")$p.value,2)
FNC_DF <- cbind(rownames(FNC_DF),FNC_DF)
colnames(FNC_DF) <- c("Sample.ID",names(FNC_DF[-1]))
## write data
cat(" Saving Data fram of fundamentalNetworkConcepts Metrics [",paste(outfile,".SampleNetwork_Stats.txt",sep=""),"]","\r","\n")
write.table(FNC_DF,file=paste(outfile,".group.",mgroup,".SampleNetwork_Stats.txt",sep=""),sep="\t",row.name=FALSE,quote=FALSE)
## write data for IAC
cat(" Saving Data fram of fundamentalNetworkConcepts Z.K outliers [",paste(outfile,".SampleNetwork_Stats_Z.K_outliers.txt",sep=""),"]","\r","\n")
write.table(FNC_DF[Z.K_outliers,],file=paste(outfile,".group.",mgroup,".SampleNetwork_Stats_Z.K_outliers.txt",sep=""),sep="\t",row.name=FALSE,quote=FALSE)
## set colours by chip of GROUPS
if(col_by_chip == 1) {   colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode)) }
if(col_by_chip == 0) {   colorvec <- gpcolors }
## plots of fundamentalNetworkConcepts
local({
colLab <<- function(n,treeorder) {
if(is.leaf(n)) {
a <- attributes(n)
i <<- i+1
attr(n, "nodePar") <-   c(a$nodePar, list(lab.col = colorvec[treeorder][i], lab.font = i%%3))
}
n
}
i <- 0
})
## Cluster for pics
cat(" Plotting SampleNetwork Metrics [",paste(outfile,".group.",mgroup,".SampleNetwork.qc.pdf",sep=""),"]","\r","\n")
meanIAC <- mean(IAC[upper.tri(IAC)],na.rm=T)
cluster1 <- hclust(as.dist(1-A.IAC),method="average")
cluster1order <- cluster1$order
cluster2 <- as.dendrogram(cluster1,hang=0.1)
cluster3 <- dendrapply(cluster2,colLab,cluster1order)
## PLOTS
## cluster IAC
pdf(file=paste(outfile,".group.",mgroup,".SampleNetwork.qc.pdf",sep=""),width=11,height=8)
par(mfrow=c(2,2))
par(mar=c(5,6,4,2))
plot(cluster3,nodePar=list(lab.cex=1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
## Connectivity
par(mar=c(5,5,4,2))
plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
abline(h=-2)
abline(h=-3)
## ClusterCoef
par(mar=c(5,5,4,2))
plot(Z.C,main="ClusterCoef", ylab="Z.C",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4,type="n")
text(Z.C,labels=samle_names,cex=0.8,col=colorvec)
abline(h=-2)
abline(h=-3)
## Connectivity vs ClusterCoef
par(mar=c(5,5,4,2))
plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=gpcolors,cex.main=1.8,cex.lab=1.4)
abline(lm(Z.C~Z.K),col="black",lwd=2)
mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
abline(v=-2,lty=2,col="grey")
abline(h=-2,lty=2,col="grey")
dev.off()
##
Z.K_outliers <- paste(Z.K_outliers,collapse=";")
res[mgroup,] <- c(mgroup,mean_IAC,n_outliers,min_Z.K,rho,rho_pvalue,Z.K_outliers)
res_out <- as.data.frame(res)
return(res_out)
}

################# if byGroup ####################

if(groups=="byGroup") {

group_list <- unique(pData(eset)$GROUPS)

n_groups <- length(group_list)

res_names <- c("Group","mean_IAC","n_outliers","min_Z.K","rho","rho_pvalue","Z.K_outliers")

res <- matrix(ncol=length(res_names),nrow=1000, dimnames=list(1:1000,res_names))

cat(" Groups [",group_list,"]","\r","\n")

	row <- 1

	for(mgroup in group_list) {

	cat(" doing group [",mgroup,"]","\r","\n")
	gx <- exprs(eset[,pData(eset)$GROUPS==mgroup])
	sample_info <- pData(eset[,pData(eset)$GROUPS==mgroup])
	samle_names <- sampleNames(eset[,pData(eset)$GROUPS==mgroup])
	groups <- as.character(sample_info$GROUPS)
	pheno <- as.character(sample_info$PHENOTYPE)
	gpcolors <- labels2colors(groups)
	## IAC, ADJACENCY & fundamentalNetworkConcepts, Z.K, Z.C, Z.MAR
	cat(" Calculating fundamentalNetworkConcepts Metrics ","\r","\n")
	IAC=cor(gx,method="p",use="p")
	diag(IAC)=0
	A.IAC=((1+IAC)/2)^2  ## ADJACENCY MATRIX
	FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
	K2=FNC$ScaledConnectivity
	Z.K=(K2-mean(K2))/sd(K2)
	Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
	Z.MAR=(FNC$MAR-mean(FNC$MAR))/sd(FNC$MAR)
	rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
	rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
	## OUTLIERS
	Z.K_outliers <- Z.K < -sd_thrs
	Z.K_outliers <- names(Z.K_outliers[Z.K_outliers==TRUE])
	n_outliers <- length(Z.K_outliers)
	mean_IAC <- mean(IAC[upper.tri(IAC)])
	min_Z.K <- min(Z.K)
	cat(" Number of Z.K outliers [", n_outliers,"]","\r","\n")
	cat(" mean_IAC [", mean_IAC,"]","\r","\n")
	## Data frame of .SampleNetwork_Stats.txt
	cat(" Making Data fram of fundamentalNetworkConcepts Metrics ","\r","\n")
	FNC_DF <- as.data.frame(FNC)
	FNC_DF$Z.K=(K2-mean(K2))/sd(K2)
	FNC_DF$Z.C=(FNC_DF$ClusterCoef-mean(FNC_DF$ClusterCoef))/sd(FNC_DF$ClusterCoef)
	FNC_DF$Z.MAR=(FNC_DF$MAR-mean(FNC_DF$MAR))/sd(FNC_DF$MAR)
	FNC_DF$mean_IAC <- mean(IAC[upper.tri(IAC)])
	FNC_DF$Mean_Connectivity <- mean(FNC$Connectivity)
	FNC_DF$Mean_ScaledConnectivity <- mean(FNC$ScaledConnectivity)
	FNC_DF$Mean_ClusterCoef <- mean(FNC$ClusterCoef)
	FNC_DF$Mean_MAR <- mean(FNC$MAR)
	FNC_DF$Decentralization <- 1-FNC_DF$Centralization
	FNC_DF$Homogeneity <- 1-FNC_DF$Heterogeneity
	FNC_DF$rho <- signif(cor.test(FNC_DF$Z.K,FNC_DF$Z.C,method="s")$estimate,2)
	FNC_DF$rho_pvalue <- signif(cor.test(FNC_DF$Z.K,FNC_DF$Z.C,method="s")$p.value,2)
	FNC_DF <- cbind(rownames(FNC_DF),FNC_DF)
	colnames(FNC_DF) <- c("Sample.ID",names(FNC_DF[-1]))
	## write data
	cat(" Saving Data fram of fundamentalNetworkConcepts Metrics [",paste(outfile,".SampleNetwork_Stats.txt",sep=""),"]","\r","\n")
	write.table(FNC_DF,file=paste(outfile,".group.",mgroup,".SampleNetwork_Stats.txt",sep=""),sep="\t",row.name=FALSE,quote=FALSE)
	## write data for IAC
	cat(" Saving Data fram of fundamentalNetworkConcepts Z.K outliers [",paste(outfile,".SampleNetwork_Stats_Z.K_outliers.txt",sep=""),"]","\r","\n")
	write.table(FNC_DF[Z.K_outliers,],file=paste(outfile,".group.",mgroup,".SampleNetwork_Stats_Z.K_outliers.txt",sep=""),sep="\t",row.name=FALSE,quote=FALSE)
	## set colours by chip of GROUPS
	colorvec <- labels2colors(as.character(pData(eset[,pData(eset)$GROUPS==mgroup])$Sentrix.Barcode))
	##if(col_by_chip == 1) {   colorvec <- labels2colors(as.character(pData(eset[,pData(eset)$GROUPS==mgroup])$Sentrix.Barcode)) }
	##if(col_by_chip == 0) {   colorvec <- gpcolors }
	## plots of fundamentalNetworkConcepts
	local({
	colLab <<- function(n,treeorder) {
	if(is.leaf(n)) {
	a <- attributes(n)
	i <<- i+1
	attr(n, "nodePar") <-   c(a$nodePar, list(lab.col = colorvec[treeorder][i], lab.font = i%%3))
	}
	n
	}
	i <- 0
	})
	## Cluster for pics
	cat(" Plotting SampleNetwork Metrics [",paste(outfile,".group.",mgroup,".SampleNetwork.qc.pdf",sep=""),"]","\r","\n")
	meanIAC <- mean(IAC[upper.tri(IAC)],na.rm=T)
	cluster1 <- hclust(as.dist(1-A.IAC),method="average")
	cluster1order <- cluster1$order
	cluster2 <- as.dendrogram(cluster1,hang=0.1)
	cluster3 <- dendrapply(cluster2,colLab,cluster1order)
	## PLOTS
	## cluster IAC
	pdf(file=paste(outfile,".group.",mgroup,".SampleNetwork.qc.pdf",sep=""),width=11,height=8)
	par(mfrow=c(2,2))
	par(mar=c(5,6,4,2))
	plot(cluster3,nodePar=list(lab.cex=1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
	mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
	## Connectivity
	par(mar=c(5,5,4,2))
	plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
	text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
	abline(h=-2);	abline(h=-3)
	## ClusterCoef
	par(mar=c(5,5,4,2))
	plot(Z.C,main="ClusterCoef", ylab="Z.C",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4,type="n")
	text(Z.C,labels=samle_names,cex=0.8,col=colorvec)
	abline(h=-2);	abline(h=-3)
	## Connectivity vs ClusterCoef
	par(mar=c(5,5,4,2))
	plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=gpcolors,cex.main=1.8,cex.lab=1.4)
	abline(lm(Z.C~Z.K),col="black",lwd=2)
	mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
	abline(v=-2,lty=2,col="grey");abline(h=-2,lty=2,col="grey")
	dev.off()

	## RESULTS

	Z.K_outliers <- paste(Z.K_outliers, collapse=";")

	res[row,] <- c(mgroup,mean_IAC,n_outliers,min_Z.K,rho,rho_pvalue,Z.K_outliers)

	row <- row + 1
	}
	      ## FINAL RESULTS
	res_out <- as.data.frame(res)
	res_out <- subset(res_out, res_out$mean_IAC!="NA")

	save(res_out,file=paste(outfile,".group.",mgroup,".SampleNetwork.qc.RData",sep=""))

	return(res_out)
	      }
}

## test
##debug( basic_sampleNetwork)
## test <- basic_sampleNetwork(eset=eset_raw,col_by_chip="0",groups="byGroup",outfile=paste(out_dir,"/",project_name,".XXXXeset_raw",sep=""),sd_thrs=2)

#################################
##  basic_sampleNetworkIterate ##
#################################

basic_sampleNetworkIterate <- function(eset, col_by_chip=c("1","0"), groups=c("all","byGroup") , outfile, IACthresh=0.95, sd_thrs=2) {
## sampleNetwork round 1
if(groups=="all") {
mgroup <- "all"
outlier_running_count <- 0;
iteration <- 1;
sd_thrs <- sd_thrs
iac_outlier_samples <- c();
n_groups <- 1
basic_sampleNetworkIterate_summary_names <- c("Group","round","nSamp","nOutlier","mean_IAC","min_Z.K","rho","rho_pvalue","Z.K_outliers")
basic_sampleNetworkIterate_summary <- matrix(ncol=length(basic_sampleNetworkIterate_summary_names),nrow=1000, dimnames=list(1:1000,basic_sampleNetworkIterate_summary_names))
out <- basic_sampleNetwork(eset,col_by_chip,groups,outfile=paste(outfile,".group.",mgroup,".round.",iteration,sep="" ))
n_samp <-  length(sampleNames(eset))
mgroup <- out$Group
mean_IAC <- out$mean_IAC;
min_Z.K <- out$min_Z.K;
rho_pvalue <- out$rho_pvalue
rho <- out$rho
Z.K_outliers <- out$Z.K_outliers_list
Z.K_outliers_list <- paste(out$Z.K_outliers)
Z.K_outliers_list <- strsplit(Z.K_outliers_list,";")[[1]]
outlier_running_count <- outlier_running_count + length(Z.K_outliers_list);
iac_outlier_samples <- c(Z.K_outliers_list, iac_outlier_samples)
res <- list(Group=mgroup,round=iteration,nSamp=n_samp,nOutlier=length(Z.K_outliers_list),mean_IAC=mean_IAC,min_Z.K=min_Z.K,KvC_rho=rho,KvC_rho_pvalue=rho_pvalue,Z.K_outliers=Z.K_outliers)
basic_sampleNetworkIterate_summary <- rbind(basic_sampleNetworkIterate_summary,res)
## if outlier samples then remove them
if(outlier_running_count >= 1) {
eset <- removeSamples_eset_lumi(eset,iac_outlier_samples)
} else {
eset <- eset;
cat(" You have no outliers!","\r","\n") ;
write.table(basic_sampleNetworkIterate_summary,file=paste(outfile,".basic_sampleNetworkIterate_summary.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
}
##
n_samp_left <- length(sampleNames(eset))
cat(" Number of outliers after round [",iteration,"] = [",outlier_running_count,"].  Percentage [",round(outlier_running_count/n_samp,3),"]. Mean IAC [",mean_IAC,"]. Min Z.K [" ,min_Z.K,"]. KvC [",rho,"] [",rho_pvalue ,"] N SAMPLE LEFT [",n_samp_left,"]","\r","\n")
iteration <- iteration + 1;
## keep going until Z.K > -2
while(iteration >= 1 &  round( as.numeric( min_Z.K ),3 )  < -sd_thrs) {
n_samp_start <- length(sampleNames(eset))
out <- basic_sampleNetwork(eset,col_by_chip,groups,outfile=paste(outfile,".round.",iteration,sep="" ))
mean_IAC <- out$mean_IAC;
min_Z.K <- out$min_Z.K;
Z.K_outliers <- out$Z.K_outliers
rho_pvalue <- out$rho_pvalue
rho <- out$rho
Z.K_outliers_list <- paste(out$Z.K_outliers)
Z.K_outliers_list <- strsplit(Z.K_outliers_list,";")[[1]]
outlier_running_count <- outlier_running_count + length(Z.K_outliers_list);
iac_outlier_samples <- c(Z.K_outliers_list, iac_outlier_samples)
eset <- removeSamples_eset_lumi(eset,iac_outlier_samples)
n_samp_left <- length(sampleNames(eset))
res <- list(Group=mgroup,round=iteration,nSamp=n_samp,nOutlier=length(Z.K_outliers_list),mean_IAC=mean_IAC,min_Z.K=min_Z.K,KvC_rho=rho,KvC_rho_pvalue=rho_pvalue,Z.K_outliers=Z.K_outliers)
basic_sampleNetworkIterate_summary <- rbind(basic_sampleNetworkIterate_summary,res)
cat(" Number of outliers after round [",iteration,"] = [",outlier_running_count,"].  Percentage [",round(outlier_running_count/n_samp,3),"]. Mean IAC [",mean_IAC,"]. Min Z.K [" ,min_Z.K,"]. KvC [",rho,"] [",rho_pvalue ,"] N SAMPLE LEFT [",n_samp_left,"]","\r","\n")
iteration <- iteration + 1;
};
final_basic_sampleNetworkIterate_summary <- as.data.frame(basic_sampleNetworkIterate_summary)
final_basic_sampleNetworkIterate_summary <- subset(final_basic_sampleNetworkIterate_summary, final_basic_sampleNetworkIterate_summary$mean_IAC!="NA")
write.table(final_basic_sampleNetworkIterate_summary,file=paste(outfile,".group.basic_sampleNetworkIterate_summary.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
outlier_samples <- iac_outlier_samples
iac_outlier_samples_df <- as.data.frame(outlier_samples)
colnames(iac_outlier_samples_df) <- c("Sample.ID")
write.table(iac_outlier_samples_df,file=paste(outfile,".group.",mgroup,".iac_outlier_samples.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
outlier_samples <- list(iac_outlier_samples=outlier_samples)
return(outlier_samples)
}

## if byGroup ##

if(groups=="byGroup") {

## make results matrix basic_sampleNetworkIterate_summary
basic_sampleNetworkIterate_summary_names <- c("Group","round","nSamp","nOutlier","mean_IAC","min_Z.K","rho","rho_pvalue","Z.K_outliers")

basic_sampleNetworkIterate_summary <- matrix(ncol=length(basic_sampleNetworkIterate_summary_names),nrow=1000, dimnames=list(1:1000,basic_sampleNetworkIterate_summary_names))

## get group names
group_list <- unique(pData(eset)$GROUPS); cat(" Number of Groups [",length(group_list),"]\nGroup Names [",group_list,"]","\r","\n")

## set outlier running count to 0
outlier_running_count <- 0;

## iteration count
iteration <- 1;

## sd threshold
sd_thrs <- sd_thrs

## outlier count
iac_outlier_samples <- c();

## n groups
n_groups <- length(group_list)

## loop through groups and run basic_sampleNetwork

 for(mgroup in group_list)  {

## subset eset to group samples
group_samples <- sampleNames(eset[,pData(eset)$GROUPS==mgroup]) ## sample names in group

group_eset <- eset[,pData(eset)$GROUPS==mgroup] ## get eset for group

cat(" doing group [",mgroup,"]","\r","\n")

## run basic_sampleNetwork
out <- basic_sampleNetwork(group_eset,col_by_chip=0,groups="byGroup",outfile=paste(outfile,".group.",mgroup,".round.",iteration,sep="" ))

## get stats from basic_sampleNetwork
n_samp <-  length(sampleNames(group_eset))
mgroup <- out$Group
mean_IAC <- out$mean_IAC;
min_Z.K <- out$min_Z.K;
rho_pvalue <- out$rho_pvalue
rho <- out$rho

## Outlier samples
Z.K_outliers <- out$Z.K_outliers

Z.K_outliers_list <- paste(out$Z.K_outliers)

Z.K_outliers_list <- strsplit(Z.K_outliers_list,";")[[1]]

## update outlier_running_count
outlier_running_count <- outlier_running_count + length(Z.K_outliers_list);

## make list of outlier samples
iac_outlier_samples <- c(Z.K_outliers_list, iac_outlier_samples)

## store in basic_sampleNetworkIterate_summary results matrix
res <- c(mgroup,iteration,n_samp,length(Z.K_outliers_list),mean_IAC,min_Z.K,rho,rho_pvalue,Z.K_outliers)

basic_sampleNetworkIterate_summary[iteration,] <- res

## if outlier samples then remove them
if(outlier_running_count >= 1) {
group_eset <- removeSamples_eset_lumi(group_eset,iac_outlier_samples)
} else {
group_eset <- group_eset;
cat(" You have no outliers!","\r","\n") ;
write.table(basic_sampleNetworkIterate_summary,file=paste(outfile,".basic_sampleNetworkIterate_summary.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
}

## number of samples remaining
n_samp_left <- length(sampleNames(group_eset))

cat(" Number of outliers after round [",iteration,"] = [",outlier_running_count,"].  Percentage [",round(outlier_running_count/n_samp,3),"]. Mean IAC [",mean_IAC,"]. Min Z.K [" ,min_Z.K,"]. KvC [",rho,"] [",rho_pvalue ,"] N SAMPLE LEFT [",n_samp_left,"]","\r","\n")

iteration <- iteration + 1; ## increment iteration number

## keep going until Z.K > -2 ##
while(iteration >= 1 &  round( as.numeric( min_Z.K ),3 )  < -sd_thrs) {

n_samp_start <- length(sampleNames(eset))

out <- basic_sampleNetwork(group_eset,col_by_chip=0,groups="byGroup",outfile=paste(outfile,".group.",mgroup,".round.",iteration,sep="" ))

mean_IAC <- out$mean_IAC;
min_Z.K <- out$min_Z.K;
Z.K_outliers <- out$Z.K_outliers
rho_pvalue <- out$rho_pvalue
rho <- out$rho

Z.K_outliers_list <- paste(out$Z.K_outliers)

Z.K_outliers_list <- strsplit(Z.K_outliers_list,";")[[1]]

outlier_running_count <- outlier_running_count + length(Z.K_outliers_list);

iac_outlier_samples <- c(Z.K_outliers_list, iac_outlier_samples)

group_eset <- removeSamples_eset_lumi(group_eset,iac_outlier_samples)

n_samp_left <- length(sampleNames(group_eset))

res <- c(mgroup,iteration,n_samp_left,length(Z.K_outliers_list),mean_IAC,min_Z.K,rho,rho_pvalue,Z.K_outliers)

basic_sampleNetworkIterate_summary[iteration,] <- res

cat(" Number of outliers after round [",iteration,"] = [",outlier_running_count,"].  Percentage [",round(outlier_running_count/n_samp,3),"]. Mean IAC [",mean_IAC,"]. Min Z.K [" ,min_Z.K,"]. KvC [",rho,"] [",rho_pvalue ,"] N SAMPLE LEFT [",n_samp_left,"]","\r","\n")

iteration <- iteration + 1;

};
##
}
## DONE LOOPING THROUGH GROUPS ##

final_basic_sampleNetworkIterate_summary <- as.data.frame(basic_sampleNetworkIterate_summary)

final_basic_sampleNetworkIterate_summary <- subset(final_basic_sampleNetworkIterate_summary, final_basic_sampleNetworkIterate_summary$mean_IAC!="NA")

write.table(final_basic_sampleNetworkIterate_summary,file=paste(outfile,".group.basic_sampleNetworkIterate_summary.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)

outlier_samples <- iac_outlier_samples

iac_outlier_samples_df <- as.data.frame(outlier_samples)

colnames(iac_outlier_samples_df) <- c("Sample.ID")

write.table(iac_outlier_samples_df,file=paste(outfile,".group.iac_outlier_samples.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

outlier_samples <- list(iac_outlier_samples=outlier_samples)

return(outlier_samples)

}

}




##################
##lumi.N <- lumiExpresso(eset_raw, variance.stabilize=TRUE,  varianceStabilize.param=list(method='vst'),  normalize.param=list(method='rsn'), bg.correct=FALSE)
##outsn <- basic_sampleNetworkIterate(eset=lumi.N,col_by_chip=0,groups="byGroup",outfile=paste(out_dir,"/",project_name,".lumi.N",sep=""),IACthresh=0.95, sd_thrs=2 );



##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

## eset_raw
## eset_bgcor
## eset_norm <- lumiExpresso(eset_bgcor, variance.stabilize=TRUE,  varianceStabilize.param=list(method='vst'),  normalize.param=list(method='rsn'), bg.correct=FALSE)
## eset_final <- eset_norm[good_probes, good_samples]

## 1 RAW qc plots, write dat, sample networks, pc regressions
## 2 bg correct qc plots, write dat, sample networks, pc regressions
## 3 Probe Detection 2SD > mean neg probe
## 4 VST qc plots, write dat, sample networks, pc regressions
## 5 RSN qc plots, write dat, sample networks, pc regressions
## 6 Batch Adjustment qc plots, write dat, sample networks, pc regressions
## 7 Final Data good probes, good samples qc plots, write dat, sample networks, pc regressions

##args <- commandArgs(trailingOnly=TRUE);

##config_file <- args[1];

config_file <- "epilepsy.config"

project_settings <- read.table(config_file, head=TRUE, sep="\t") ## as.is=T, fill=TRUE)

t(project_settings)

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

cat("\n"," START PRE_PROCESSING!" ,"\r","\n")

cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")

#########################
## 1: Project settings ##
#########################

cat(" Getting Project settings and file names","\r","\n")

project_dir <- project_settings$project_dir ## /scratch/project/pipelines/DATA/expression/GAP_Expression
gs_report <- project_settings$gs_report ## "Sample_and_Control_Probe_Profile_FinalReport.txt"
gs_probe <- project_settings$gs_probe ## "Group_Probe_Profile_Final_Report.txt" # genomestudio report
gs_sample <- project_settings$gs_sample ## "sample_table_Final_Report.txt" # genomestudio report
gs_control <- project_settings$gs_control ## "control_probe_profile_Final_Report.txt" # genomestudio report

anno_table <- project_settings$anno_table ## annotation table for project
pheno_file <- project_settings$pheno_file ## FILE NAME must contain : Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID
project_name <- project_settings$project_name ## "GAP" # GAP
out_dir <- paste(project_name,"_lumi_processing" ,sep="") # output dir for lumi process and plots

probe_det <- project_settings$probe_det ## 80 # 50, 80, 90,100
sample_det <- project_settings$sample_det ## 80 # 50,80,90,100

sex_check <- project_settings$sex_check ## 1 # 1 or 0
iac_check <- project_settings$iac_check ## 1 # 1 or 0
iac_sd_thrs <- project_settings$iac_sd_thrs ## 2 #

norm_method <- project_settings$norm_method ## "rsn" # quantile, rsn, or both
transform_method <- project_settings$transform_method ## "vst" # log2, vst or both
mbcb_method <- project_settings$mbcb_method ## NP or MLE

tech_pheno_file <- project_settings$tech_pheno_file ## "sample_tech_info.txt"

cat(" Data Dir ",project_dir,"\r","\n")
cat(" Set Working Dir to [", project_dir,"]","\r","\n")

setwd(project_dir)

cat(" Project Name:- [", project_name,"]" ,"\r","\n")
cat(" Merged Final Report:- [", gs_report,"]" ,"\r","\n")
cat(" Probe Final Report:- [", gs_probe,"]" ,"\r","\n")
cat(" Sample Table Final Report:- [", gs_sample,"]" ,"\r","\n")
cat(" Control Probe Final Report:- [", gs_control,"]" ,"\r","\n")
cat(" Probe Annotations Final Report:- [", anno_table,"]" ,"\r","\n")
cat(" Phenotype File:- [", pheno_file,"]" ,"\r","\n")
cat(" Batch & Technical Phenotype File:- [", tech_pheno_file,"]" ,"\r","\n")
cat(" Probe Detection Call Rate:- [", probe_det,"]" ,"\r","\n")
cat(" Sample Detection Call Rate:- [", sample_det,"]" ,"\r","\n")
cat(" Sex Check:- [",sex_check,"]","\n","\r")
cat(" IAC Check:- [",iac_check,"]","\n","\r")
cat(" IAC Check SD Threshold:- [",iac_sd_thrs,"]","\n","\r")
cat(" Transformation Method:- [",transform_method,"]","\r","\n")
cat(" Normalisation Method:- [",norm_method,"]","\r","\n")

cat(" Making processing directory:- ",paste(" mkdir ./",out_dir,sep=""),"\r","\n")

make_dir_command <- paste(" if [ ! -e ./",out_dir," ]; then mkdir ./",out_dir,"; fi",sep="")

system( make_dir_command )

cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

####################
## summary matrix ##
####################
eset_summary <- matrix(ncol=4,nrow=6, dimnames=list(c("eset_raw","eset_raw","eset_bkcor","eset_lumiT","eset_lumiN","eset_final"),c("eset_step","n_samples","n_probes","IAC")) )
eset_summary[,1] <- rownames(eset_summary )

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

########################
### Read sample table ##
########################

cat(" Reading Sample Table Final Report generated in Genomestudio [", gs_sample ,"]","\r","\n")

raw_sample_data <- read.table(paste(gs_sample) ,skip=8,as.is=T,fill=T,head=T,sep="\t")

rownames(raw_sample_data) <- raw_sample_data$Sample.ID

raw_sample_data <- raw_sample_data[,names(raw_sample_data)!="X"]

raw_n_samples <- dim(raw_sample_data)[1]

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

###########################################################
## if phenotype file is provided add this to sample data ##
###########################################################

if(is.na(pheno_file))  stop(" WARNING!: YOU HAVENT PROVIDED ANY PHENOTYPE INFORMATION!!!" )

cat(" Phenotype Data Provided.[",pheno_file,"] Mergeing with Sample Table Final Report generated in Genomestudio [", gs_sample,"]" ,"\r","\n")

pheno_dat <- read.table(paste(pheno_file), as.is=T,fill=T,head=T,sep="\t")

has_pheno_cols <- c("Sample.ID","SEX","GROUPS","TISSUE","PHENOTYPE","Study_ID") %in% names(pheno_dat);

missing_pheno_cols <- "FALSE" %in% has_pheno_cols

if(missing_pheno_cols == "TRUE") stop(" WARNING!: YOU ARE MISSING ESSENTIAL SAMPLE INFORMATION! MAKE SURE YOUR PHENO_FILE HAS:- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID !!!")

raw_sample_data <- merge(raw_sample_data , pheno_dat , by.x="Sample.ID",by.y="Sample.ID",sort=FALSE,all.x=TRUE)

rownames(raw_sample_data) <- raw_sample_data$Sample.ID

write.table(raw_sample_data,file=paste(out_dir,"/",project_name,".sample.pheno.txt",sep=""),sep="\t",quote=F,row.names=F )

cat(" Number of samples = [",raw_n_samples,"]","\r","\n","\n")

###########################################################
## if tech Batch file is provided add this to sample data ##
###########################################################

if(is.na(tech_pheno_file))  stop(" WARNING!: YOU HAVENT PROVIDED ANY BATCH INFORMATION!!!" )

cat(" Technical/Batch Data Provided.[",tech_pheno_file,"] Mergeing with Sample Table Final Report generated in Genomestudio [", gs_sample,"]" ,"\r","\n")

tech_pheno <- read.table(paste(tech_pheno_file),head=T,sep="\t")

tech_pheno$Sentrix.Barcode <- as.character(tech_pheno$Sentrix.Barcode)

colnames(tech_pheno) <- paste("tech.",names(tech_pheno),sep="")

raw_sample_data <- merge(raw_sample_data , tech_pheno , by.x="Sample.ID",by.y="tech.Sample.ID",sort=FALSE,all.x=TRUE)

rownames(raw_sample_data) <- raw_sample_data$Sample.ID

write.table(raw_sample_data,file=paste(out_dir,"/",project_name,".eset_raw.sample.pheno.txt",sep=""),sep="\t",quote=F,row.names=F )

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

###############################################################
## Readging Genomestudio Final Report file & Calling lumiR() ##
###############################################################

cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" Reading Genomestudio Final Report file :[",gs_report,"]","\r","\n")
cat(" Calling lumiR() on [", gs_report,"]","\r","\n")
cat(" Saved as object [ eset_raw ] ","\r","\n")

eset_raw <- lumiR(paste(gs_report), lib.mapping="lumiHumanIDMapping",annotationColumn = c('PROBE_ID','CHROMOSOME','SYMBOL','DEFINITION','ACCESSION','ENTREZ_GENE_ID','PROBE_TYPE','PROBE_START','PROBE_SEQUENCE','PROBE_CHR_ORIENTATION','PROBE_COORDINATES','CHROMOSOME','TRANSCRIPT','ILMN_GENE','REFSEQ_ID','UNIGENE_ID','SYMBOL','PROTEIN_PRODUCT'))

cat(" Getting Chip Info ","\r","\n","\n")

chip <- getChipInfo(eset_raw)

cat(" chip version: [", chip$chipVersion,"]", "\r","\n")
cat(" chip species: [", chip$species,"]", "\r","\n")
cat(" chip IDType: [", chip$IDType,"]", "\r","\n")
cat(" chip chipProbeNumber: [", chip$chipProbeNumber,"]", "\r","\n")
cat(" chip inputProbeNumber: [", chip$inputProbeNumber,"]", "\r","\n")
cat(" chip matchedProbeNumber: [", chip$matchedProbeNumber,"]", "\r","\n")
cat(" Number of samples: [",raw_n_samples,"]","\r","\n")

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

##########################################
## adding Sample Report to pData() slot ##
##########################################

cat(" adding Sample Report to pData() slot of [eset_raw]","\r","\n")

pData(eset_raw) <- raw_sample_data;

################################################
## getting probe annotaions from fData() slot ##
################################################

if( !is.na(anno_table) ) {

cat(" getting probe annotations from Final Report ","\r","\n")

system( paste("  cut -f 1-24 ", anno_table," | awk 'NR > 8 ' > ",anno_table,".tmp",sep="" ) ) ## cut -f 1-24 probe_annotation_Final_Report.txt | awk 'NR > 8' | less -S

tmp_anno <- paste(anno_table,".tmp",sep="")

probe_annotations <- read.table(tmp_anno ,as.is=T,fill=T,head=T,sep="\t")

fData(eset_raw)$nuID <- rownames(fData(eset_raw))
fData(eset_raw)$forder <- 1:dim(fData(eset_raw)[1])

n_probe_anno <- dim(probe_annotations)[1]
n_probe_fdata <- dim(fData(eset_raw))[1]

cat(" N probes read [",n_probe_anno,"]","\r","\n")
cat(" N probes fData [",n_probe_fdata,"]","\r","\n")


probe_annotations$nuID <- probeID2nuID( probe_annotations$ProbeID , lib.mapping = "lumiHumanIDMapping")[,"nuID"]

probe_annotations$order <- 1:dim(probe_annotations)[1]

probe_annotations <- merge(probe_annotations, fData(eset_raw)[,c("nuID","forder")] ,by.x="nuID",by.y="nuID", sort=FALSE, all.y=TRUE )

probe_annotations <- probe_annotations[order(probe_annotations$order),names(probe_annotations)!="forder"]

rownames(probe_annotations) <- probe_annotations$nuID

fData(eset_raw) <- probe_annotations

n_probe_fdata <- dim(fData(eset_raw))[1]
cat(" N probes new fData [",n_probe_fdata,"]","\r","\n")

system(paste(" rm ",anno_table,".tmp",sep="") )

} else {

cat(" getting probe annotations from fData() slot of [ eset_raw ]","\r","\n")

probe_annotations <- fData(eset_raw)

probe_annotations$nuID <- rownames(probe_annotations)

fData(eset_raw) <- probe_annotations

}

cat(" saving probe annotaions as [",paste(out_dir,"/",project_name,".eset_raw.probe_annotations.RData",sep=""),"]","\r","\n")

save(probe_annotations,file=paste(out_dir,"/",project_name,".eset_raw.probe_annotations.RData",sep="") )

write.table(probe_annotations,file=paste(out_dir,"/",project_name,".eset_raw.probe_annotations.txt",sep=""),sep="\t",quote=F,row.names=F )

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

##########################
## Save eset_raw #########
##########################
cat(" Saving eset_raw as [", paste(out_dir,"/",project_name,".eset_raw.RData",sep=""),"]","\r","\n")

save(eset_raw , file=paste(out_dir,"/",project_name,".eset_raw.RData",sep="") )## This is raw data pre-qc

#########################################
## beadNum, detection, exprs, se.exprs ##
#########################################
cat(" Writing eset_raw [beadNum, detection, exprs, se.exprs] data to file ", paste(out_dir,"/",project_name,".eset_raw.[beadNum, detection, exprs, se.exprs].txt",sep=""), "\r","\n")

write_expression_files(eset=eset_raw,outfile= paste(out_dir,"/",project_name,".eset_raw",sep="") )

####################
## PLOTS RAW DATA ##
####################
gx_qc_plots_lumi(eset=eset_raw, outfile= paste(out_dir,"/",project_name,".eset_raw.qc.plots.pdf",sep="") ,do_pca=TRUE );

#9249896080_C 9257622046_J
#   -3.768176    -8.802562
##############################################################
## background correct using mbcb                            ##
##############################################################

cat(" Start background correction ", "\r","\n")

sig <- paste(out_dir,"/",project_name,".eset_raw_exprsSignal.txt",sep="")
neg <- paste(out_dir,"/",project_name,".eset_raw_negative.control.txt",sep="")

#############################################################
## get negative control bead data
#############################################################
cat(" get negative control bead data","\r","\n")

negativeControl <- getControlData(eset_raw)
negativeControl <- subset(negativeControl, negativeControl$controlType=="NEGATIVE")
negativeControl_orig <- negativeControl[,c(3:dim(negativeControl)[2])]

#############################################################
## removing outlier neg beads and replace with mean
#############################################################
negativeControl <- apply(negativeControl_orig ,2, negBeadOutlierRepMean)

neg_max <- apply(negativeControl,2,max)
neg_sd <- apply(negativeControl,2,sd)
neg_mean <- apply(negativeControl,2,mean)

#############################################################
## save negativeControl bead expression data
#############################################################
write.table(negativeControl, file=neg, sep="\t");

#############################################################
## get expression data
#############################################################
cat(" get expression bead data ","\r","\n")

expressionSignal <- exprs(eset_raw)

#############################################################
## save expression data
#############################################################
write.table(expressionSignal, file=sig, sep="\t");

#############################################################
## set data for mbcb
#############################################################
cat(" set data for mbcb ","\r","\n")

data <- mbcb.parseFile(sig, neg);

signal <- data$sig;
negCon <- data$con;

#############################################################
## run mbcb
#############################################################
cat(" run background correct using mbcb  ","\r","\n")

gx_mbcb <- mbcb.correct(signal,negCon,npBool=TRUE,mleBool=TRUE, isRawBead=FALSE)

###############
## save mbcb ##
###############
cat(" saveing raw mbcb data to ",paste(out_dir,"/",project_name,".eset_raw.mbcb.RData",sep=""), "\r","\n")

save(gx_mbcb , file=paste(out_dir,"/",project_name,".eset_raw.mbcb.correct.RData",sep=""))

cat(" mbcb complete ","\r","\n")

################################
## select mbcb method results ##
################################
cat(" Selected mbcb method to save: ",mbcb_method,"\r","\n")

if(mbcb_method=="NP") { gx_mbcb <- gx_mbcb$NP }

if(mbcb_method=="MLE") { gx_mbcb <- gx_mbcb$MLE }


#############################################################
## replace names with original as R adds "X" to numbers #####
#############################################################

cat(" replace names with original sampleNames(eset_raw), as R adds X to numbers ","\r","\n")

colnames(gx_mbcb) <- sampleNames(eset_raw)

#############################################################
## make new eset for bk corrected data
#############################################################

cat(" Creating Background Corrected Data set: eset_bkcor  ","\r","\n")

eset_bkcor <- eset_raw

#############################################################
## replace old exprs data with new mbcb bk corrected data
#############################################################
cat(" replace old exprs data with new mbcb Background corrected data ","\r","\n")

exprs(eset_bkcor) <- as.matrix(gx_mbcb)

cat(" Saveing Background Corrected Data set: ",paste(out_dir,"/",project_name,".eset_bkcor.RData",sep=""),"\r","\n")

####################################
## SAVE BACKGROUND CORRECTED DATA ##
####################################
save(eset_bkcor, file=paste(out_dir,"/",project_name,".eset_bkcor.RData",sep="")  , compress=T)

########################################
## beadNum, detection, exprs, se.exprs #
########################################
cat(" Writing eset_bkcor [beadNum, detection, exprs, se.exprs] data to file ", paste(out_dir,"/",project_name,".eset_bkcor.[beadNum, detection, exprs, se.exprs].txt",sep=""), "\r","\n")

write_expression_files(eset=eset_bkcor,outfile= paste(out_dir,"/",project_name,".eset_bkcor",sep="") )

#############################
## PLOTS POST bkcor        ##
#############################

gx_qc_plots_lumi(eset=eset_bkcor, outfile= paste(out_dir,"/",project_name,".eset_bkcor.qc.plots.pdf",sep="") ,do_pca=TRUE );
##9249896080_C 9257622046_J
##   -3.768696    -8.802189
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

##################################
## Sex Check based on XIST GENE ##
##################################
cat(" Checking Gender based on XIST gene expression","\r","\n")

#####################################################
## get neg control data from eset_bkcor
#####################################################
negativeControl <- getControlData(eset_bkcor)
negativeControl <- subset(negativeControl, negativeControl$controlType=="NEGATIVE")
negativeControl <- negativeControl[,c(3:dim(negativeControl)[2])]

#####################################################
## get neg control info mean,sd,max etc
#####################################################
neg_max <- apply(negativeControl,2,max)
neg_sd <- apply(negativeControl,2,sd)
neg_mean <- apply(negativeControl,2,mean)
neg_2sd <- neg_mean + 2*neg_sd

#####################################################
## get XIST gene postion
#####################################################
xist_raw <- fData(eset_raw)$ILMN_GENE=="XIST";

xist_bkcor <- fData(eset_bkcor)$ILMN_GENE=="XIST";

#####################################################
## get XIST gene expression signal
#####################################################
xist_gx_raw  <- exprs(eset_raw[xist_raw, ]  )
xist_gx_raw <- as.data.frame( t(xist_gx_raw));

xist_gx_bkcor <- exprs(eset_bkcor[xist_bkcor , ]  )
xist_gx_bkcor <- as.data.frame( t(xist_gx_bkcor));

#####################################################
## cobine raw and bkCro gx data
#####################################################
xist_gx <- cbind(xist_gx_raw,xist_gx_bkcor)

colnames(xist_gx) <- c("raw_XIST","bkcor_XIST")

xist_gx$neg_2sd <- neg_2sd

xist_gx$neg_max <- neg_max

xist_gx$neg_mean <- neg_mean

xist_gx$neg_sd <- neg_sd

#####################################################
## gender based on XIST expression 1=female , 0 = male
#####################################################
xist_gx$XIST_Gender_max <- ifelse(xist_gx$bkcor_XIST > xist_gx$neg_max,1,0)

xist_gx$XIST_Gender_2sd <- ifelse(xist_gx$bkcor_XIST > xist_gx$neg_2sd,1,0)  ### THIS WORKS!! IF GENE IS > 2SD_NEG_BEAD THEN ITS EXPRESSED/DETECTED!

xist_gx$XIST_z <-  ( xist_gx$bkcor_XIST - xist_gx$neg_mean ) / xist_gx$neg_sd

xist_gx$Sample.ID <- rownames(xist_gx)

head(xist_gx)

xist_gx$xist_gender <- ifelse(xist_gx$bkcor_XIST > xist_gx$neg_2sd,"Female","Male")

head(xist_gx)

pdf(file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.pdf",sep=""),height=8,width=11)
plot(xist_gx$bkcor_XIST, type="p", col=2+xist_gx$XIST_Gender_2sd, pch=19 , cex=0.6);points(xist_gx$neg_2sd,col="blue", type="l")
dev.off()

save(xist_gx, file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.RData",sep=""))

write.table(xist_gx, file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.txt",sep=""),sep="\t",row.names=F)

#####################################################
## check sex
#####################################################


###pheno_update <- cbind(pData(eset_bkcor),xist_gx[,names(xist_gx)!="Sample.ID"]])

pheno_update <- merge(pData(eset_bkcor), xist_gx, by.x="Sample.ID", by.y="Sample.ID",all.x=TRUE,sort=FALSE)

pData(eset_bkcor) <- pheno_update



pData(eset_bkcor)$gender_missmatch <- ifelse( pData(eset_bkcor)$xist_gender == pData(eset_bkcor)$SEX, "PASS","GENDER_MISSMATCH" )

save(eset_bkcor, file=paste(out_dir,"/",project_name,".eset_bkcor.RData",sep=""))

n_gender_fails <- sum( pData(eset_bkcor)$gender_missmatch=="GENDER_MISSMATCH" )

if ( n_gender_fails > 0 ) {
    cat(" WARNING: Youn have GENDER_MISSMATCH samples!!!!!!! N=[",n_gender_fails,"]","\r","\n")
} else {
  cat(" Congratulations! \n All your Males are Male and Females are Female. \n You have NO GENDER_MISSMATCH samples!!!", "\r","\n")
}

#############################
## write file of sex fails ##
#############################
gender_missmatch_table <- subset(pData(eset_bkcor),pData(eset_bkcor)$gender_missmatch=="GENDER_MISSMATCH")

write.table(gender_missmatch_table, file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.gender_missmatch_table.txt",sep=""),sep="\t",row.names=F)

####################################
## SAVE BACKGROUND CORRECTED DATA sex checked ##
####################################
save(eset_bkcor, file=paste(out_dir,"/",project_name,".eset_bkcor.RData",sep="")  , compress=T)
cat(" Writing eset_bkcor [beadNum, detection, exprs, se.exprs] data to file ", paste(out_dir,"/",project_name,".eset_bkcor.[beadNum, detection, exprs, se.exprs].txt",sep=""), "\r","\n")
write_expression_files(eset=eset_bkcor,outfile= paste(out_dir,"/",project_name,".eset_bkcor",sep="") )


##--------------------------------------------------------------------------------------------------------------------------------------------------------#

##############################################
## PROBE DETECTED 2SD ABOVE MEAN BACKGROUND ##
##############################################

cat(" Calculating Probe Detection rates. \n Probe is seen as Detected if it has background corrected signal intensity >= 2SD of the mean intensity of the negative control beads","\r","\n")

#####################################################
## get expression matrix
#####################################################
gx <- exprs(eset_bkcor)

#####################################################
## get negative bead ranges mean or max or >2SD mean  of neg beads
#####################################################
neg_2sd <- neg_mean + 2*neg_sd

#######################################################################
## sweep through gx matrix to id probes > 2SD Mean of negative beads
#######################################################################
det <- sweep(gx, 1, neg_2sd,">=")

##########################################
## Writing Probe Detection Calls to file #
##########################################
det_tmp  <- as.data.frame(det)

det_tmp$nuID <- rownames(det_tmp)

probe_detected <- merge(fData(eset_bkcor), det_tmp , by.x="nuID",by.y="nuID",sort=FALSE)

cat(" Writing Probe Detection Calls to",paste(out_dir,"/",project_name,".eset_bkcor.probe_detected.txt",sep="") ,"\r","\n")

write.table(probe_detected , file=paste(out_dir,"/",project_name,".eset_bkcor.probe_detected.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

rm("det_tmp")

#####################################################
## probe_detection counts
#####################################################
probe_detection <- rowSums(det)

##########################################
## n samples
##########################################
n_samples <- dim(gx)[2]

#####################################################
## probe annotations
#####################################################

probes_not_detected_in_any_sample <- probe_detection==0

probes_detected_in_80_sample <- probe_detection>=n_samples*0.80

probes_detected_in_all_sample <- probe_detection==n_samples

probe_annotations_0_detected <- fData(eset_bkcor[probes_not_detected_in_any_sample,])
probe_annotations_80_detected  <- fData(eset_bkcor[probes_detected_in_80_sample ,])
probe_annotations_100_detected  <- fData(eset_bkcor[probes_detected_in_all_sample,])

cat(" Adding detetion call rate for all probes and samples to fData() slot for eset_bkcor","\r","\n")

fData(eset_bkcor)$n_detected <- probe_detection
fData(eset_bkcor)$n_detected_call_rate <- round( probe_detection/n_samples ,3)

#####################################################
## sample_detection counts
#####################################################
n_probes <- dim(eset_bkcor)[1]

sample_detection <- colSums(det)

pData(eset_bkcor)$n_probes_detected <- sample_detection

pData(eset_bkcor)$n_probes_detected_call_rate <- round( sample_detection/n_probes ,3)

##plot(pData(eset_bkcor)$n_probes_detected_call_rate)
##plot(pData(eset_bkcor)$n_probes_detected ,pData(eset_bkcor)$Detected.Genes..0.01.)

save(eset_bkcor, file=paste(out_dir,"/",project_name,".eset_bkcor.RData",sep=""))

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

#####################################################
## get group information from pData() slot
#####################################################
group_names <- unique(pData(eset_bkcor)$EPI_TYPE);

groups  <- pData(eset_bkcor)$EPI_TYPE

n_groups <- length(group_names)

###########################
## get expression matrix ##
###########################
gx <- exprs(eset_bkcor)

neg_2sd <- as.numeric(neg_2sd)

##########################################################################################
## loop through each group and id probes > 2SD mean neg beads in X% of samples/group
##########################################################################################
for(n in group_names ) {

cat(" Finding probes in ",probe_det," of sample group [",n,"] with signal intensity >= 2SD mean intensity of the negative control beads ","\r","\n")

group_label <- paste(n)

sel_samples <- pData(eset_bkcor)$EPI_TYPE==n;

n_samples_in_group <- dim(gx[,sel_samples])[2];

cat(" Number of samples in group [",n,"] = ",n_samples_in_group,"\r","\n")

detection_matrix <- sweep(gx[,sel_samples],1,neg_2sd[sel_samples],">=")

group_probe_detection <- rowSums(detection_matrix) >= probe_det*n_samples_in_group

group_probe_detection_nuID <- rownames( gx[group_probe_detection, ])

cat(" Number of probes in group [",n,"] with signal intensity >= 2SD mean intensity of the negative control beads = ", length(group_probe_detection_nuID) ,"\r","\n")

cat(" Writing probe list to ",paste(out_dir,"/",project_name,".GROUP.",group_label,".detected_probes_nuID.txt",sep=""), "\r","\n" )

det_probes <- as.data.frame(group_probe_detection_nuID)

colnames(det_probes) <- c("nuID")

write.table( det_probes ,file=paste(out_dir,"/",project_name,".GROUP.",group_label,".detected_probes_nuID.txt",sep=""),row.names=FALSE,quote=FALSE,col.names=FALSE)

}

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

#####################################################
## Y CHROM EXPRESSION IN XIST MALES
#####################################################
cat(" Y Chromosome probe detection based on XIST males","\r","\n")

xist_males <- pData(eset_bkcor)$xist_gender=="Male"

gx_y <- exprs(eset_bkcor[fData(eset_bkcor)$CHR=="Y",])

detection_matrix_y <- sweep( gx_y[,xist_males],1, neg_2sd[xist_males ],">=")

y_probe_detection <- rowSums(detection_matrix_y) >= probe_det * sum(xist_males==TRUE)

y_probe_detection_nuID <- rownames( gx_y[y_probe_detection, ])

y_det_probes <- as.data.frame(y_probe_detection_nuID)

colnames(y_det_probes) <- c("nuID")

write.table( y_det_probes ,file=paste(out_dir,"/",project_name,".GROUP.Y.detected_probes_nuID.txt",sep=""),row.names=FALSE,quote=FALSE,col.names=FALSE)

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

#####################################################
## writing final good probe list
######################################################

cat(" writing final good probe list to ",paste(out_dir,"/",project_name,".detected_probes_nuID_final.txt",sep=""),"\r","\n")

system(paste("cat ",out_dir,"/",project_name,"****.detected_probes_nuID.txt | sort | uniq >> ",out_dir,"/",project_name,".detected_probes_nuID_final.txt", sep="") )

good_probes <- read.table(file=paste(out_dir,"/",project_name,".detected_probes_nuID_final.txt",sep=""),head=FALSE)

good_probes <- paste(good_probes[,1])

n_good_probes <- length(good_probes)

cat(" Total number of good probes = ",n_good_probes,"\n","\r" )

good_probes_annotation <- fData(eset_bkcor[paste(good_probes,sep=""),])

good_probes_annotation$raw_mean <- apply( exprs(eset_bkcor[good_probes,]), 1,mean)
good_probes_annotation$raw_sd <- apply( exprs(eset_bkcor[good_probes,]), 1,sd)
good_probes_annotation$raw_var <- apply( exprs(eset_bkcor[good_probes,]), 1,var)
good_probes_annotation$raw_min <- apply( exprs(eset_bkcor[good_probes,]), 1,min)
good_probes_annotation$raw_max <- apply( exprs(eset_bkcor[good_probes,]), 1,max)

head(good_probes_annotation)

cat(" saving good probe annotations to ",paste(out_dir,"/",project_name,".detected_probes_nuID_final.***",sep=""),"\r","\n")

save(good_probes_annotation, file=paste(out_dir,"/",project_name,".detected_probes_nuID_final.RData",sep="") )

write.table(good_probes_annotation, file=paste(out_dir,"/",project_name,".detected_probes_nuID_final.txt",sep=""),quote=F,sep="\t",row.names=F )

#--------------------------------------------------------------------------------------------------------------------------------------------------------#



##################
## lumiExpresso ##
##################

eset_norm <- lumiExpresso(eset_bkcor, variance.stabilize=TRUE,  varianceStabilize.param=list(method='vst'),  normalize.param=list(method='rsn'), bg.correct=FALSE)
save(eset_norm, file=paste(out_dir,"/",project_name,".eset_norm.RData",sep="")  , compress=T)
write_expression_files(eset=eset_norm,outfile= paste(out_dir,"/",project_name,".eset_norm",sep="") )
gx_qc_plots_lumi(eset=eset_norm, outfile= paste(out_dir,"/",project_name,".eset_norm.qc.plots.pdf",sep="") ,do_pca=TRUE );



#########################
## basic_sampleNetwork ##
#########################

basic_sampleNetworkIterate_eset_norm <- basic_sampleNetworkIterate(eset=eset_norm,col_by_chip=0,groups="byGroup",outfile=paste(out_dir,"/",project_name,".eset_norm.",sep=""),IACthresh=0.95, sd_thrs=2 );




####################################################################
## Transform
####################################################################

cat(" TRANSFORM DATA. Method =",transform_method,"\r","\n")

if(transform_method=="vst") {
  eset <- eset_bkcor
  eset_lumiT  <- lumiT(eset, ifPlot = FALSE, method="vst" )
}
if(transform_method=="log2") {
  eset <-  eset_bkcor
  eset_lumiT  <- lumiT(eset, ifPlot = FALSE, method="log2" )
}

################################
## SAVE transform
################################
save(eset_lumiT , file=paste(out_dir,"/",project_name,".eset_lumiT.RData",sep="") )

########################################
## beadNum, detection, exprs, se.exprs #
########################################
cat(" Writing eset_lumiT [beadNum, detection, exprs, se.exprs] data to file ", paste(out_dir,"/",project_name,".eset_lumiT.[beadNum, detection, exprs, se.exprs].txt",sep=""), "\r","\n")

gx <- exprs(eset_lumiT)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiT),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiT.exprs.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- se.exprs(eset_lumiT)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiT),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiT.se.exprs.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- detection(eset_lumiT)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiT),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiT.detection.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- beadNum(eset_lumiT)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiT),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiT.beadNum.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

###############################
## PLOTS OF TRANSFORMED DATA ##
###############################
qc_plots_pdf <- paste(out_dir,"/",project_name,".eset_lumiT.qc.plots.pdf",sep="")
chip_col <- labels2colors( as.character(pData(eset_lumiT  )$Sentrix.Barcode))
cat(" Start basic QC plots of eset_lumiT. Saved to: ",qc_plots_pdf,"\r","\n")
datExprs0 <- exprs(eset_lumiT)
cat(" calculate pca for plots ")
pca_raw <- prcomp(t(datExprs0))$x
cat(" plotting... ","\r","\n")
pdf(qc_plots_pdf, width=16.5,height=11.7) ## A3
plot(eset_lumiT, what='boxplot', main=paste(project_name," eset_raw Boxplot",sep=""), col=chip_col )
plot(eset_lumiT, what='density', main=paste(project_name," eset_raw Density",sep="") )
plot(eset_lumiT, what='cv',main=paste(project_name," eset_raw density plot of coefficient of varience",sep="")  )
scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()
cat(" finished basic QC plots of eset_lumiT Saved to: ",qc_plots_pdf,"\r","\n")
rm("datExprs0","tmp_iac","pca_raw");
gc()
###############################
## IAC OF TRANSFORMED DATA ##
###############################
cat(" Checking IAC after eset_lumiT ","\n","\n")
rm("iac_outlierSamples")
datExprs0 <- exprs(eset_lumiT)
pdf(file=paste(out_dir,"/",project_name,".eset_lumiT.iac_outliers.pdf",sep=""), width=8,height=6 )
iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection
plot(iac_outlierSamples$samplemeanIAC)
dev.off()
eset_lumiT_IAC <- iac_outlierSamples$meanIAC
cat(" meanIAC post VST: ",  eset_lumiT_IAC,"\r","\n")
nsamp <- dim(exprs(eset_lumiT))[2]
nprobe <- dim(exprs(eset_lumiT))[1]
IAC <- round(iac_outlierSamples$meanIAC, 4)
eset_summary["eset_lumiT",2:4] <- c(nsamp,nprobe,IAC)

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

####################################################################
# normalise
####################################################################

cat(" NORMALISE DATA. Method =",norm_method,"\r","\n")

if(norm_method=="rsn") {
  eset_lumiN  <- lumiN(eset_lumiT, ifPlot = FALSE, method="rsn" )
}

if(norm_method=="quantile") {
  eset_lumiN  <- lumiN(eset_lumiT, ifPlot = FALSE, method="quantile" )
}

########################################
## SAVE normalise
########################################
save(eset_lumiN , file=paste(out_dir,"/",project_name,".eset_lumiN.RData",sep="") )

########################################
## beadNum, detection, exprs, se.exprs #
########################################
cat(" Writing eset_lumiT [beadNum, detection, exprs, se.exprs] data to file ", paste(out_dir,"/",project_name,".eset_lumiT.[beadNum, detection, exprs, se.exprs].txt",sep=""), "\r","\n")

gx <- exprs(eset_lumiN)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiN),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiN.exprs.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- se.exprs(eset_lumiN)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiN),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiN.se.exprs.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- detection(eset_lumiN)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiN),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiN.detection.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- beadNum(eset_lumiN)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_lumiN),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_lumiN.beadNum.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

###############################
## PLOTS OF normalised DATA ##
###############################
qc_plots_pdf <- paste(out_dir,"/",project_name,".eset_lumiN.qc.plots.pdf",sep="")
chip_col <- labels2colors( as.character(pData(eset_lumiN  )$Sentrix.Barcode))
cat(" Start basic QC plots of eset_lumiN. Saved to: ",qc_plots_pdf,"\r","\n")
datExprs0 <- exprs(eset_lumiN)
cat(" calculate pca for plots ")
pca_raw <- prcomp(t(datExprs0))$x
cat(" plotting... ","\r","\n")
pdf(qc_plots_pdf, width=16.5,height=11.7) ## A3
plot(eset_lumiN, what='boxplot', main=paste(project_name," eset_raw Boxplot",sep=""), col=chip_col )
plot(eset_lumiN, what='density', main=paste(project_name," eset_raw Density",sep="") )
plot(eset_lumiN, what='cv',main=paste(project_name," eset_raw density plot of coefficient of varience",sep="")  )
scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()
cat(" finished basic QC plots of eset_lumiN Saved to: ",qc_plots_pdf,"\r","\n")
rm("datExprs0","tmp_iac","pca_raw");
gc()
###############################
## IAC OF normalised DATA ##
###############################
cat(" Checking IAC after eset_lumiN ","\n","\n")
rm("iac_outlierSamples")
datExprs0 <- exprs(eset_lumiN)
pdf(file=paste(out_dir,"/",project_name,".eset_lumiN.iac_outliers.pdf",sep=""), width=8,height=6 )
iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection
plot(iac_outlierSamples$samplemeanIAC)
dev.off()
eset_lumiN_IAC <- iac_outlierSamples$meanIAC
cat(" meanIAC post normalise: ",  eset_lumiN_IAC,"\r","\n")
nsamp <- dim(exprs(eset_lumiN))[2]
nprobe <- dim(exprs(eset_lumiN))[1]
IAC <- round(iac_outlierSamples$meanIAC, 4)
eset_summary["eset_lumiN",2:4] <- c(nsamp,nprobe,IAC)
eset_summary

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

################################################
## make final data set : remove bad probes
################################################

eset_final <- eset_lumiN[good_probes,]

#####################
## SAVE FINAL DATA ##
#####################
save(eset_final , file=paste(out_dir,"/",project_name,".eset_final.RData",sep="") )

########################################
## beadNum, detection, exprs, se.exprs #
########################################
cat(" Writing eset_final [beadNum, detection, exprs, se.exprs] data to file ", paste(out_dir,"/",project_name,".eset_final.[beadNum, detection, exprs, se.exprs].txt",sep=""), "\r","\n")

gx <- exprs(eset_final)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_final),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_final.exprs.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- se.exprs(eset_final)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_final),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_final.se.exprs.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- detection(eset_final)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_final),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_final.detection.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

gx <- beadNum(eset_final)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset_final),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(out_dir,"/",project_name,".eset_final.beadNum.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

###############################
## PLOTS OF FINAL DATA ##
###############################
qc_plots_pdf <- paste(out_dir,"/",project_name,".eset_final.qc.plots.pdf",sep="")
chip_col <- labels2colors( as.character(pData(eset_final  )$Sentrix.Barcode))
cat(" Start basic QC plots of eset_lumiN. Saved to: ",qc_plots_pdf,"\r","\n")
datExprs0 <- exprs(eset_final)
cat(" calculate pca for plots ")
pca_raw <- prcomp(t(datExprs0))$x
cat(" plotting... ","\r","\n")
pdf(qc_plots_pdf, width=16.5,height=11.7) ## A3
plot(eset_final, what='boxplot', main=paste(project_name," eset_final Boxplot",sep=""), col=chip_col )
plot(eset_final, what='density', main=paste(project_name," eset_final Density",sep="") )
plot(eset_final, what='cv',main=paste(project_name," eset_raw density plot of coefficient of varience",sep="")  )
scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()
cat(" finished basic QC plots of eset_final Saved to: ",qc_plots_pdf,"\r","\n")
rm("datExprs0","tmp_iac","pca_raw");
gc()
###############################
## IAC OF FINAL DATA ##
###############################
cat(" Checking IAC after eset_lumiN ","\n","\n")
rm("iac_outlierSamples")
datExprs0 <- exprs(eset_final)
pdf(file=paste(out_dir,"/",project_name,".eset_final.iac_outliers.pdf",sep=""), width=8,height=6 )
iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection
plot(iac_outlierSamples$samplemeanIAC)
dev.off()
eset_final_IAC <- iac_outlierSamples$meanIAC
cat(" meanIAC post normalise and removal of bad probes: ",  eset_final_IAC,"\r","\n")
nsamp <- dim(exprs(eset_final))[2]
nprobe <- dim(exprs(eset_final))[1]
IAC <- round(iac_outlierSamples$meanIAC, 4)
eset_summary["eset_final",2:4] <- c(nsamp,nprobe,IAC)
##eset_summary
##iac_outlierSamples$samples_to_remove
##summary(iac_outlierSamples$samples_to_remove)
##[1] "dataClean"         "IAC"               "meanIAC"
##[4] "meanIACdiag"       "samplemeanIAC"     "numbersd"
##[7] "clust"             "samples_to_remove"

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

########################
## BATCH REGRESSIONS ##
#######################

## making this a requirement if [tech_pheno_file] then do this
## test for assoc of PHENO with PC1 and batches, test for association of PC1 with batches

## tech_pheno_file <- "sample_tech_info.txt"

if( !is.na(tech_pheno_file)==TRUE ) {

cat(" Starting batch versus PC1 and PHENOTYPE Rregressions ","\n","\r")

cat(" Getting Gene expression matrix ","\n","\r")

gx <- exprs(eset_lumiN)

gx <- t(gx)

cat(" Reading in technical information on eg [Sample.ID, RIN, RNA_YIELD, BATCH, CHIP, DATE_CHIP_RUN, DATE_RNA_EXTRACTED] ",tech_pheno_file ,"\n","\r")

tech_pheno <- read.table(paste(tech_pheno_file),head=TRUE,sep="\t") ## Sample.ID,RIN,RNA_YIELD,BATCH,CHIP, DATE_CHIP_RUN,DATE_RNA_EXTRACTED

tech_pheno$Sentrix.Barcode <- as.factor(tech_pheno$Sentrix.Barcode)

pdata <- pData(eset_lumiN)

pdata <- as.data.frame(pdata[,c("Sample.ID","Index")])

colnames(pdata) <- c("Sample.ID","Index")

tech_pheno <- merge( pdata, tech_pheno, by.x="Sample.ID", by.y="Sample.ID", sort=FALSE, all.x=TRUE)

tech_batch <- tech_pheno[,3:dim(tech_pheno)[2]]

##head(tech_batch)

######################
## get names of var ##
######################
batch_var_names <- names(tech_batch)

##########################
## which ones are dates ##
##########################
date_vars <- grep("Date",batch_var)

#############
## Run PCA ##
#############
cat(" Running PCA on t(exprs(eset)) ","\n","\r")

pca_gx <- prcomp(gx)$x

pca_gx <- pca_gx[,"PC1"]

################################
## PHENOTYPES FOR REGRESSIONS ##
################################

PC1 <- as.numeric(pca_gx)

PHENOTYPE <- as.numeric(as.factor(pData(eset_lumiN)$PHENOTYPE))

GROUPS <-    as.numeric(as.factor(pData(eset_lumiN)$GROUPS))

#################################################
## Test for association of batch vars with PC1 ##
#################################################

#############################
## multivariate full model ##
#############################

multivariate_model_terms <- paste(batch_var_names ,collapse="+")

for(pheno in c("PC1","PHENOTYPE","GROUPS") ) {

multivariate_model <- paste(pheno,"~",multivariate_model_terms,sep="") ## prot ~ c1+c2+c3...

multivariate_model <- as.formula(multivariate_model)

cat(" running full multivariate models ",pheno," ~ multivariate_model ", "\r","\n")

######################
## multivariate lm  ##
######################
lm_batch <- lm(multivariate_model, data=tech_batch)

#########################
## RSQUARED summary lm ##
#########################
lm_r2 <- round( summary(lm_batch)$adj.r.squared, 3)

###############
## summary lm #
###############
summary_lm_batch <- summary(lm_batch)$coef

summary_lm_batch <- as.data.frame(summary_lm_batch)

summary_lm_batch$terms <- rownames(summary_lm_batch)

summary_lm_batch$significant <- ifelse(summary_lm_batch$"Pr(>|t|)"<=0.05,1,0)

summary_lm_batch$model_rsq <- lm_r2

##########################
## save summary lm #######
##########################
write.table(summary_lm_batch , file=paste(out_dir,"/",project_name,".eset_final.",pheno,".multivariate_model_batch_variables.csv",sep=""),row.names=FALSE,quote=FALSE, sep=",")

###########################
## multivariate ANOVA    ##
###########################
anova_lm_batch <- anova(lm_batch)

anova_lm_data <- as.data.frame(anova_lm_batch)

anova_lm_data$terms <- rownames(anova_lm_data)

anova_lm_data <- subset(anova_lm_data, anova_lm_data$terms!="Residuals")

##################
## plot ANOVA P ##
##################
pdf(file=paste(out_dir,"/",project_name,".eset_final.",pheno,".ANOVA_multivariate_model_batch_variables.pdf",sep=""), width=8,height=6 )
par(mar=c(10,5,4,2))
barplot(-log10(anova_lm_data $"Pr(>F)"), las=3, names=c(anova_lm_data$terms) , ylab="ANOVA -log10(P)", main=paste("multivariate_model. R2=",lm_r2,sep=""),cex.names=0.8,cex.main=0.8,cex.lab=1)
abline(h=-log10(0.05),col="blue")
abline(h=-log10(0.01),col="red")
dev.off()

###############
## save ANOVA #
###############
write.table(anova_lm_data , file=paste(out_dir,"/",project_name,".eset_final.",pheno,".ANOVA_multivariate_model_batch_variables.csv",sep=""),row.names=FALSE,quote=FALSE, sep=",")

############################################################
## test for sig ANOVA terms ################################
############################################################

min_anova_p <- min(anova_lm_data$"Pr(>F)")

if(min_anova_p <= 0.05) {

sig_anova_lm_data <- subset(anova_lm_data, anova_lm_data$"Pr(>F)"<=0.05)
most_sig_term <- sig_anova_lm_data$terms[sig_anova_lm_data$"Pr(>F)"==min_anova_p]
sig_terms <- sig_anova_lm_data$terms
sig_terms <- paste(sig_terms ,collapse=" + ")

cat(" WARNING!: SIGNIFICANT ASSOCIATION BETWEEN [ ",pheno," ] ~ [",sig_terms," ]. R2=",lm_r2 ," [",most_sig_term,"] MIN_P=",min_anova_p,". YOU MAY WANT TO CORRECT FOR THIS BEFORE THE FINAL ANLYSIS ","\r","\n")

##################################################
## STEP find best terms ##########################
##################################################

cat(" Finding a best model by AIC in a Stepwise Algorithm ","\r","\n" )

step_lm_batch <- stepAIC(lm_batch,direction = "both")

#####################
## summary step lm ##
#####################
summary_step_lm_batch <- summary(step_lm_batch)$coef

summary_step_lm_batch  <- as.data.frame(summary_step_lm_batch)

summary_step_lm_batch$terms <- rownames(summary_step_lm_batch)

###########################
## save summary step lm ###
###########################
write.table(summary_step_lm_batch, file=paste(out_dir,"/",project_name,".eset_final.stepAIC_multivariate_model_batch_variables.csv",sep=""),row.names=FALSE,quote=FALSE, sep=",")

###############
## anova step #
###############
anova_step_lm_batch <- anova(step_lm_batch)

anova_data <- as.data.frame(anova_step_lm_batch)

anova_data$terms <- rownames(anova_data)

anova_data <- subset(anova_data, anova_data$terms!="Residuals")

best_model <- paste(anova_data$terms, collapse=" + ")

best_model <- as.formula( paste(pheno,"~",best_model,sep="") )

###################
## RSQ anova step #
###################
anova_step_r2 <- round( summary(step_lm_batch)$adj.r.squared, 3)

######################
## plot ANOVA step P #
######################
pdf(file=paste(out_dir,"/",project_name,".eset_final.",pheno,".stepANOVA_multivariate_model_batch_variables.pdf",sep=""), width=8,height=6 )
par(mar=c(10,5,4,2))
barplot(-log10(anova_data$"Pr(>F)"), las=3, names=c(anova_data$terms) , ylab="ANOVA -log10(P)", main=paste(step_lm_batch$call[2]," R2=",anova_step_r2 ,sep=""),cex.names=0.8,cex.main=0.8,cex.lab=1)
abline(h=-log10(0.05),col="blue")
abline(h=-log10(0.01),col="red")
dev.off()

## save stepANOVA
write.table(anova_data , file=paste(out_dir,"/",project_name,".eset_final.",pheno,".stepANOVA_multivariate_model_batch_variables.csv",sep=""),row.names=FALSE,quote=FALSE, sep=",")

cat(" BEST MODEL BASED ON stepAIC [",pheno,"] ~ [",best_model,"] ","\r","\n")

} else {

cat(" NICE!: NO SIGNIFICANT ASSOCIATION BETWEEN [ ",pheno," ] ~ [ BATCH PHENOTYPES ]","\r","\n")

}

#######################
## univariate models ##
#######################
cat(" running full univariate models ","\r","\n")

for(pheno in c("PC1","GROUPS","PHENOTYPE") ) {

univar_results <- data.frame()

	for(covar in batch_var_names) {

		pheno_name <- paste(pheno,sep="")

		model <- as.formula( paste(pheno,"~",covar) )

		univar_lm <- lm( model, data=tech_batch)

		summary_univar_lm <- summary(univar_lm)

		rsq <- summary_univar_lm$adj.r.squared

		anova_univar_lm <- anova(univar_lm)

		phenotype <- pheno_name
		covariate <- paste(covar,sep="")

		summary_univar_lm_data <- as.data.frame(coef(summary_univar_lm))
		summary_univar_lm_data$adj.r.squared  <- rsq
		summary_univar_lm_data$anova_F  <- anova_univar_lm$F[1]
		summary_univar_lm_data$anova_p <- anova_univar_lm$"Pr(>F)"[1]
		summary_univar_lm_data$terms <- rownames(summary_univar_lm_data)
		summary_univar_lm_data$phenotype <- pheno_name
		summary_univar_lm_data$covariate <- covariate

		univar_results <- rbind(univar_results , summary_univar_lm_data )

		if( !is.na(anova_univar_lm$"Pr(>F)"[1]) & anova_univar_lm$"Pr(>F)"[1]<=0.05) { cat(" WARNING!: [",pheno,"] ~ [", covar,"] ARE ASSOCIATED ","\r","\n" ) }

	}

write.table(univar_results, file=paste(out_dir,"/",project_name,".eset_final.",pheno,".univariate_model_batch_variables.csv",sep="")  ,sep=",",row.names=FALSE,quote=FALSE)

}

} else {

cat(" WARNING!: SKIPPING BATCH ASSOCIATION TESTING. NOT REALLY RECOMMENDED!","\r","\n")

}

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

######################
## SAVE ESET_SUMMARY #
######################

cat(" Writing summary information on pre-processing steps to ",paste(out_dir,"/",project_name,".eset_summary.txt",sep=""),"\r","\n")

save(eset_summary, file=paste(out_dir,"/",project_name,".eset_summary.RData",sep=""))

write.table(eset_summary, file=paste(out_dir,"/",project_name,".eset_summary.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

eset_summary <- as.data.frame(eset_summary)

print(eset_summary)

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

#########################
## END PRE_PROCESSING! ##
#########################

cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" END PRE_PROCESSING!" ,"\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

###########
## CHMOD ##
###########
system( paste(" chmod 776 ",out_dir,"/",project_name,"***",sep=""))

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

########################
## PRINT SESSION INFO ##
########################
sessionInfo()

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

#######################
## CLEAN UP AND QUIT ##
#######################
rm(list=ls())
gc()
q()
n

#--------------------------------------------------------------------------------------------------------------------------------------------------------#











