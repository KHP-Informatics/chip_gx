#!/share/bin/Rscript
rm(list=ls())
######################################################################################################################################################################################
## Authors : Stephen Newhouse, Amos Forlarin
## Email : stephen..j.newhouse@gmail.comp
## Date : 15 Nov 2013
## Verison: 1.00
######################################################################################################################################################################################

######################################################################################################################################################################################

cat(" Loading Libraries ","\r","\n")

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
source("./pre_process_gx.R")
source("./SampleNetwork_1.0.r")
source("./ModuleSampleNetwork_0.5.r")

######################################################################################################################################################################################

## some functins
negBeadOutlierRepMean <- function(x) { 
		z_out_samp <- abs(  as.numeric( scale(x) )  ) > 2
		mean_pop <- mean(x[z_out_samp==FALSE])
		sd_pop <- sd(x[z_out_samp==FALSE])
		new_x <- ifelse( abs(  as.numeric( scale(x) )  ) > 2, mean_pop, x ) 
		return(new_x)
}
######################################################################################################################################################################################

options(stringsAsFactors = FALSE)

##args <- commandArgs(trailingOnly=TRUE);
##  config_file <- args[1];

config_file <- "project.config"

project_settings <- read.table(config_file, head=T,sep="\t",as.is=T,fill=TRUE)

######################################################################################################################################################################################

#########################
## 1: Project settings ##
#########################

cat(" Getting Project settings and file names","\r","\n")

project_dir <- project_settings$project_dir ## /scratch/project/pipelines/DATA/expression/GAP_Expression
gs_report <- project_settings$gs_report ## "Sample_and_Control_Probe_Profile_FinalReport.txt"
gs_probe <- project_settings$gs_probe ## "Group_Probe_Profile_Final_Report.txt" # genomestudio report
gs_sample <- project_settings$gs_sample ## "sample_table_Final_Report.txt" # genomestudio report
gs_control <- project_settings$gs_control ## "control_probe_profile_Final_Report.txt" # genomestudio report
anno_table <- project_settings$anno_table ## "NULL" "annotation table for project
pheno_file <- project_settings$pheno_file ## "NULL" # FILE NAME OR "NULL"study specific information 
project_name <- project_settings$project_name ## "GAP" # GAP
out_dir <- paste(project_name,"_lumi_processing" ,sep="") # output dir for lumi process and plots
probe_det <- project_settings$probe_det ## 80 # 50, 80, 90,100
sample_det <- project_settings$sample_det ## 80 # 50,80,90,100
sex_check <- project_settings$sex_check ## 1 # 1 or 0
iac_check <- project_settings$iac_check ## 1 # 1 or 0
iac_sd_thrs <- project_settings$iac_sd_thrs ## 2 #
norm_method <- project_settings$norm_method ## "rsn" # quantile, rsn, or both
transform_method <- project_settings$transform_method ## "vst" # log2, vst or both

## mbcb_method=="NP"
mbcb_method=="MLE"


cat(" Data Dir ",project_dir,"\r","\n")
cat(" Set Working Dir to ", project_dir,"\r","\n")

setwd(project_dir)

cat(" Project Name:- ", project_name ,"\r","\n")
cat(" Merged Final Report:- ", gs_report ,"\r","\n")
cat(" Probe Final Report:- ", gs_probe ,"\r","\n")
cat(" Sample Table Final Report:- ", gs_sample ,"\r","\n")
cat(" Control Probe Final Report:- ", gs_control ,"\r","\n")
cat(" Probe Annotations Final Report:- ", anno_table ,"\r","\n")
cat(" Phenotype File:- ", pheno_file ,"\r","\n")
cat(" Probe Detection Call Rate:- ", probe_det ,"\r","\n")
cat(" Sample Detection Call Rate:- ", sample_det ,"\r","\n")
cat(" Sex Check:- ",sex_check,"\n","\r")
cat(" IAC Check:- ",iac_check,"\n","\r")
cat(" IAC Check SD Threshold:- ",iac_sd_thrs,"\n","\r")
cat(" Transformation Method:- ",transform_method,"\r","\n")
cat(" Normalisation Method:- ",norm_method,"\r","\n")

cat(" Making processing directory. mkdir ./",out_dir,"\r","\n")

system( paste(" mkdir ./",out_dir,sep="") )

## groups data ??


######################################################################################################################################################################################

########################
### Read sample table ##
########################

cat(" Reading Sample Table Final Report generated in Genomestudio ", gs_sample ,"\r","\n")

raw_sample_data <- read.table(gs_sample ,skip=8,as.is=T,fill=T,head=T,sep="\t")

rownames(raw_sample_data) <- raw_sample_data$Sample.ID

raw_n_samples <- dim(raw_sample_data)[1]

###########################################################
## if phenotype file is provided add this to sample data ##
###########################################################

if( !is.na(pheno_file) ) {  

cat(" Phenotype Data Provided. Mergeing with Sample Table Final Report generated in Genomestudio ", gs_sample ,"\r","\n")

pheno_dat <- read.table(pheno_file, as.is=T,fill=T,head=T,sep="\t")

raw_sample_data <- merge(raw_sample_data , pheno_dat , by.x="Sample.ID",by.y="Sample.ID",sort=FALSE,all.x=TRUE)  

rownames(raw_sample_data) <- raw_sample_data$Sample.ID

write.table(raw_sample_data,file=paste(out_dir,"/",project_name,".sample.pheno.txt",sep=""),sep="\t",quote=F,row.names=F )

} else {

cat(" No Phenotype Data Provided! Moving on :( .....","\r","\n")

}

cat(" Number of samples = ",raw_n_samples,"\r","\n")

######################################################################################################################################################################################

############################################################
## Create Raw LumiBatch object                            ##
############################################################

###############################################################
## Readging Genomestudio Final Report file & Calling lumiR() ##
###############################################################

cat(" Readging Genomestudio Final Report file ","\r","\n")
cat(" Calling lumiR() on ", gs_report,"\r","\n")

raw_lumiBatch <- lumiR(gs_report, lib.mapping="lumiHumanIDMapping",
annotationColumn = c('PROBE_ID','CHROMOSOME','SYMBOL','DEFINITION','ACCESSION','ENTREZ_GENE_ID','PROBE_TYPE','PROBE_START','PROBE_SEQUENCE','PROBE_CHR_ORIENTATION','PROBE_COORDINATES',
'CHROMOSOME','TRANSCRIPT','ILMN_GENE','REFSEQ_ID','UNIGENE_ID','SYMBOL','PROTEIN_PRODUCT'))

cat( print(raw_lumiBatch) )

cat(" Getting Chip Info ","\r","\n")

chip <- getChipInfo(raw_lumiBatch)

cat(" chip version: ", chip$chipVersion, "\r","\n")
cat(" chip species: ", chip$species, "\r","\n")
cat(" chip IDType: ", chip$IDType, "\r","\n")
cat(" chip chipProbeNumber: ", chip$chipProbeNumber, "\r","\n")
cat(" chip inputProbeNumber: ", chip$inputProbeNumber, "\r","\n")
cat(" chip matchedProbeNumber: ", chip$matchedProbeNumber, "\r","\n")
cat(" Number of samples = ",raw_n_samples,"\r","\n")

## add sample data/information from genome studio report or external
## this should be the sample report from genomestudio plus any useful demographic data and phenotype data ie age,sex,case/control status, tissue type, phenotype etc etc 

##########################################
## adding Sample Report to pData() slot ##
##########################################

cat(" adding Sample Report to pData() slot","\r","\n")

pData(raw_lumiBatch) <- raw_sample_data;

################################################
## getting probe annotaions from fData() slot ##
################################################

if( !is.na(anno_table) ) { 

cat(" getting probe annotations from Final Report ","\r","\n")

system( paste("  cut -f 1-24 ", anno_table," | awk 'NR > 8 ' > ",anno_table,".tmp",sep="" ) ) ## cut -f 1-24 probe_annotation_Final_Report.txt | awk 'NR > 8' | less -S

tmp_anno <- paste(anno_table,".tmp",sep="")

probe_annotations <- read.table(tmp_anno ,as.is=T,fill=T,head=T,sep="\t") 

probe_annotations$nuID <- probeID2nuID( probe_annotations$ProbeID , lib.mapping = "lumiHumanIDMapping")[,"nuID"]

rownames(probe_annotations) <- probe_annotations$nuID

fData(raw_lumiBatch) <- probe_annotations

system(paste(" rm ",anno_table,".tmp",sep="") )

} else { 

cat(" getting probe annotations from fData() slot ","\r","\n")

probe_annotations <- fData(raw_lumiBatch) 

probe_annotations$nuID <- rownames(probe_annotations)

}

cat(" saving probe annotaions as",paste(out_dir,"/",project_name,".raw_lumiBatch.probe_annotations.RData",sep=""),"\r","\n")

save(probe_annotations,file=paste(out_dir,"/",project_name,".raw_lumiBatch.probe_annotations.RData",sep="") )

write.table(probe_annotations,file=paste(out_dir,"/",project_name,".raw_lumiBatch.probe_annotations.txt",sep=""),sep="\t",quote=F,row.names=F )

##########################
## Save raw_lumiBatch ####
##########################

cat(" Saving raw_lumiBatch as ", paste(out_dir,"/",project_name,".raw_lumiBatch.RData",sep=""), "\r","\n")

save(raw_lumiBatch , file=paste(out_dir,"/",project_name,".raw_lumiBatch.RData",sep="") )## This is raw data pre-qc

######################################################################################################################################################################################

####################
## PLOTS RAW DATA ##
####################

library(scatterplot3d)

raw_qc_plots_pdf <- paste(out_dir,"/",project_name,".raw_lumiBatch.qc.plots.pdf",sep="")

chip_col <- labels2colors( as.character(pData(raw_lumiBatch)$Sentrix.Barcode))

cat(" Start basic QC plots of raw_lumiBatch. Saved to: ",raw_qc_plots_pdf,"\r","\n")

datExprs0 <- exprs(raw_lumiBatch) 

cat(" calculate pca for plots ")

pca_raw <- prcomp(t(datExprs0))$x

cat(" plotting... ","\r","\n")

pdf(raw_qc_plots_pdf, width=16.5,height=11.7) ## A3 
	plot(raw_lumiBatch, what='boxplot', main=paste(project_name," raw_lumiBatch Boxplot",sep=""), col=chip_col )
	plot(raw_lumiBatch, what='density', main=paste(project_name," raw_lumiBatch Density",sep="") )
	plot(raw_lumiBatch, what='cv',main=paste(project_name," raw_lumiBatch density plot of coefficient of varience",sep="")  )
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
	tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()

cat(" finished basic QC plots of raw_lumiBatch. Saved to: ",raw_qc_plots_pdf,"\r","\n")

rm("datExprs0","tmp_iac","pca_raw");

gc()

######################################################################
## IAC QC                                                           ##
## ## optional :- but worth one round of outlier detection          ##
######################################################################

if(iac_check==1) {

cat(" Starting inter-array correlation quality control step ","\n","\r")

datExprs0 <- exprs(raw_lumiBatch) 

pdf(file=paste(out_dir,"/",project_name,".raw_lumiBatch.iac_outliers.pdf",sep=""), width=8,height=6 )

iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection

plot(iac_outlierSamples$samplemeanIAC)

dev.off()

## meanIAC 
cat(" meanIAC: ", iac_outlierSamples$meanIAC , "\r","\n")

## add iac sd to pData
cat(" addingd IAC sd to pData() slot of raw_lumiBatch ", "\r","\n")

pData(raw_lumiBatch)$iac_sd <- iac_outlierSamples$numbersd

## iac outliers
iac_outliers <- iac_outlierSamples$samples_to_remove

## n iac out 
n_iac_outliers <- length(iac_outlierSamples$samples_to_remove)

cat(" number of iac outliers: ", n_iac_outliers , "\r","\n")

## iac sample information
poorChips <- iac_outlierSamples$samples_to_remove

poorChips <- data.frame(Sample.ID=names(iac_outliers),iac_sd=poorChips,samplemeanIAC=iac_outlierSamples$samplemeanIAC[names(iac_outliers)],QC=rep("iac_outliers", n_iac_outliers)  )

save(poorChips, file=paste(out_dir,"/",project_name,".IAC.outliers.RData",sep="") )

cat(" saving list of iac outliers to: ",paste(out_dir,"/",project_name,".IAC.outliers.txt",sep=""), "\r","\n")

write.table(poorChips,file=paste(out_dir,"/",project_name,".IAC.outliers.txt",sep=""),quote=F,row.names=F,sep="\t")

##  remove iac samples
sel_samp <- (sampleNames(raw_lumiBatch) %in% names(iac_outliers) )==FALSE

sel_samp_names <- sampleNames(raw_lumiBatch)[sel_samp]

cat(" removing iac outliers from raw_lumiBatch","\r","\n")

raw_lumiBatch <- raw_lumiBatch[,sel_samp]

cat(" update Control Probe data slot..subsetting to remaining samples ", "\r","\n")

## update control data slot
update_controlData <- getControlData(raw_lumiBatch)

update_controlData  <- update_controlData[,c("controlType","ProbeID",sel_samp_names)]

raw_lumiBatch <- addControlData2lumi(update_controlData , raw_lumiBatch)

cat(" saving updated raw_lumiBatch post IAC to: ",paste(out_dir,"/",project_name,".raw_lumiBatch_postIAC.RData",sep=""), "\r","\n" )

save(raw_lumiBatch , file=paste(out_dir,"/",project_name,".raw_lumiBatch_postIAC.RData",sep="") )

cat(" Number of samples after iac outlier removal: ", length(sampleNames(raw_lumiBatch)), "\r","\n" )

} else {

cat(" Skipping inter-array correlation quality control step ","\n","\r")

}

#############################
## PLOTS POST IAC RAW DATA ##
#############################

raw_qc_plots_pdf <- paste(out_dir,"/",project_name,".raw_lumiBatch.post.iac.qc.plots.pdf",sep="")

chip_col <- labels2colors( as.character(pData(raw_lumiBatch)$Sentrix.Barcode))

cat(" Start basic QC plots of post IAC raw_lumiBatch. Saved to: ",raw_qc_plots_pdf,"\r","\n")

datExprs0 <- exprs(raw_lumiBatch) 

cat(" calculate pca for plots post IAC ")

pca_raw <- prcomp(t(datExprs0))$x

cat(" plotting post IAC... ","\r","\n")

pdf(raw_qc_plots_pdf, width=16.5,height=11.7) ## A3 
	plot(raw_lumiBatch, what='boxplot', main=paste(project_name," raw_lumiBatch Boxplot",sep=""), col=chip_col )
	plot(raw_lumiBatch, what='density', main=paste(project_name," raw_lumiBatch Density",sep="") )
	plot(raw_lumiBatch, what='cv',main=paste(project_name," raw_lumiBatch density plot of coefficient of varience",sep="")  )
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
	tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()

cat(" finished basic QC plots of post IAC raw_lumiBatch. Saved to: ",raw_qc_plots_pdf,"\r","\n")

rm("datExprs0","tmp_iac","pca_raw")
gc()


#########################################
## Check IAC after 1st round of iac qc ##
#########################################

cat(" Checking IAC after 1st round of IAC QC ","\n","\n")

preIAC <- iac_outlierSamples$meanIAC

rm("iac_outlierSamples")

datExprs0 <- exprs(raw_lumiBatch) 

pdf(file=paste(out_dir,"/",project_name,".raw_lumiBatch.post.iac_outliers.pdf",sep=""), width=8,height=6 )

iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection

plot(iac_outlierSamples$samplemeanIAC)

dev.off()

postIAC <- iac_outlierSamples$meanIAC

cat(" meanIAC pre  IAC QC: ", preIAC, "\r","\n") 
cat(" meanIAC post IAC QC: ", postIAC,"\r","\n")

##############################################################
## background correct using mbcb                            ##
##############################################################

cat(" Start background correction ", "\r","\n")

sig <- paste(out_dir,"/",project_name,".raw_lumiBatch_exprsSignal.txt",sep="")
neg <- paste(out_dir,"/",project_name,".raw_lumiBatch_negative.control.txt",sep="")

## get negative control bead data
cat(" get negative control bead data","\r","\n")

negativeControl <- getControlData(raw_lumiBatch)
negativeControl <- subset(negativeControl, negativeControl$controlType=="NEGATIVE")
negativeControl_orig <- negativeControl[,c(3:dim(negativeControl)[2])]

## removing outlier neg beads and replace with mean
negativeControl <- apply(negativeControl_orig ,2, negBeadOutlierRepMean)

neg_max <- apply(negativeControl,2,max)
neg_sd <- apply(negativeControl,2,sd)
neg_mean <- apply(negativeControl,2,mean)


## save negativeControl bead expression data
write.table(negativeControl, file=neg, sep="\t");

## get expression data
cat(" get expression bead data ","\r","\n")

expressionSignal <- exprs(raw_lumiBatch)

## save expression data
write.table(expressionSignal, file=sig, sep="\t");

## set data for mbcb
cat(" set data for mbcb ","\r","\n")

data <- mbcb.parseFile(sig, neg);

signal <- data$sig;
negCon <- data$con;
     
## run mbcb
cat(" run background correct using mbcb  ","\r","\n")

gx_mbcb <- mbcb.correct(signal,negCon,npBool=TRUE,mleBool=TRUE, isRawBead=FALSE)

cat(" saveing raw mbcb data to ",paste(out_dir,"/",project_name,".raw_lumiBatch.mbcb.RData",sep=""), "\r","\n")

save(gx_mbcb , file=paste(out_dir,"/",project_name,".raw_lumiBatch.mbcb.correct.RData",sep=""))

if(mbcb_method=="NP") { gx_mbcb <- gx_mbcb$NP }

if(mbcb_method=="MLE") { gx_mbcb <- gx_mbcb$MLE }

cat(" mbcb complete ","\r","\n")

## replace names with original as R adds "X" to numbers

cat(" replace names with original sampleNames(raw_lumiBatch) as R adds X to numbers ","\r","\n")

colnames(gx_mbcb) <- sampleNames(raw_lumiBatch)

## make new eset for bk corrected data

cat(" Creating Background Corrected Data set: bgCor_lumiBatch ","\r","\n")

bkCor_lumiBatch <- raw_lumiBatch

## replace old exprs data with new mbcb bk corrected data
cat(" replace old exprs data with new mbcb Background corrected data ","\r","\n")

exprs(bkCor_lumiBatch) <- as.matrix(gx_mbcb)

cat(" Saveing Background Corrected Data set: ",paste(out_dir,"/",project_name,".bkCor_lumiBatch.RData",sep=""),"\r","\n")

save(bkCor_lumiBatch, file=paste(out_dir,"/",project_name,".bkCor_lumiBatch.RData",sep="")  , compress=T)

#############################
## PLOTS POST bkCor        ##
#############################

bkCor_qc_plots_pdf <- paste(out_dir,"/",project_name,".bkCor_lumiBatch.post.bkCor.qc.plots.pdf",sep="")

chip_col <- labels2colors( as.character(pData(bkCor_lumiBatch)$Sentrix.Barcode))

cat(" Start basic QC plots of post IAC raw_lumiBatch. Saved to: ",bkCor_qc_plots_pdf,"\r","\n")

datExprs0 <- exprs(bkCor_lumiBatch) 

cat(" calculate pca for plots post bkCor ")

pca_raw <- prcomp(t(datExprs0))$x

cat(" plotting post IAC... ","\r","\n")

pdf(bkCor_qc_plots_pdf, width=16.5,height=11.7) ## A3

	plot(bkCor_lumiBatch, what='boxplot', main=paste(project_name," bkCor_lumiBatch Boxplot",sep=""), col=chip_col )
	plot(bkCor_lumiBatch, what='density', main=paste(project_name," bkCor_lumiBatch Density",sep="") )
	plot(bkCor_lumiBatch, what='cv',main=paste(project_name," bkCor_lumiBatch density plot of coefficient of varience",sep="")  )
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
	tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()

cat(" finished basic QC plots of post bkCor on raw_lumiBatch. Saved to: ",bkCor_qc_plots_pdf,"\r","\n")

rm("datExprs0","tmp_iac","pca_raw")
gc()


#########################################
## Check IAC after 1st round of iac qc ##
#########################################

cat(" Checking IAC after bkCor ","\n","\n")

rm("iac_outlierSamples")

datExprs0 <- exprs(bkCor_lumiBatch) 

pdf(file=paste(out_dir,"/",project_name,".bkCor_lumiBatch.iac_outliers.pdf",sep=""), width=8,height=6 )

iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection

plot(iac_outlierSamples$samplemeanIAC)

dev.off()

postbkCorIAC <- iac_outlierSamples$meanIAC

cat(" meanIAC pre  IAC QC: ", preIAC, "\r","\n") 
cat(" meanIAC post IAC QC: ", postIAC,"\r","\n")
cat(" meanIAC post bkCor: ",  postbkCorIAC,"\r","\n")


##################################
## Sex Check based on XIST GENE ##
##################################

neg_max <- apply(negativeControl,2,max)

neg_sd <- apply(negativeControl,2,sd)

neg_mean <- apply(negativeControl,2,mean)

neg_P95 <- apply(negativeControl,2, quantile(probs=0.95) )

xist_raw <- fData(raw_lumiBatch)$ILMN_GENE=="XIST";
xist_bkCor <- fData(bkCor_lumiBatch)$ILMN_GENE=="XIST";


xist_gx_raw  <- exprs(raw_lumiBatch[xist_raw, ]  ) 
xist_gx_raw <- as.data.frame( t(xist_gx_raw)); 

xist_gx_bkCor <- exprs(bkCor_lumiBatch[xist_bkCor , ]  ) 
xist_gx_bkCor <- as.data.frame( t(xist_gx_bkCor)); 

## z diff ##
xist_gx_bkCor$z_diff <- (xist_gx_bkCor - neg_mean)/neg_sd 
xist_gx_bkCor$neg_mean <- neg_mean
xist_gx_bkCor$neg_max <- neg_max

xist_gx_raw$z_diff <- (xist_gx_raw - neg_max)/neg_sd 


xist_gx$Sample.ID <- rownames(xist_gx)

colnames(xist_gx) <- c("Raw_XIST","Sample.ID")

raw_sample_data$XIST <- xist_gx$Raw_XIST

## detecton call based on illumina p value

## addinng infot to raw_sample_data
## XIST
raw_sample_data$bgXIST <- xist_gx$XIST

## NEG BEAD MAX,MEDIAN,MEAN,SD 
neg_max <- apply(negativeControl,2,max)
neg_sd <- apply(negativeControl,2,sd)
neg_mean <- apply(negativeControl,2,mean)

(xist - neg_mean)/neg_sd 

#raw_sample_data$neg_max <- apply(negativeControl,2,max)
#raw_sample_data$neg_sd <- apply(negativeControl,2,sd)
#raw_sample_data$neg_med <- apply(negativeControl,2,median)

raw_sample_data$xist_gx <- ifelse(raw_sample_data$bgXIST > raw_sample_data$neg_max,1,0) ### THIS WORKS!! IF GENE IS > MAX_NEG_BEAD THEN ITS EXPRESSED/DETECTED!

raw_sample_data$xist_gender <- ifelse(raw_sample_data$bgXIST > raw_sample_data$neg_max,"Female","Male")

pdf(file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.pdf",sep=""))
plot( raw_sample_data$neg_max, raw_sample_data$XIST, col=(raw_sample_data$xist_gx+2), ylab="XIST",xlab="MAX_NEG_BEAD", main="Green=Expressed above background")
plot( log2(raw_sample_data$XIST), col=(raw_sample_data$xist_gx+2), ylab="XIST",xlab="Index", main="Green=Expressed above background")
#points(1:raw_n_samples, log2(raw_sample_data$neg_max), col="black")
dev.off()

pData(bkCor_lumiBatch) <- raw_sample_data




############################
## PROBE ABOVE BACKGROUND ##
############################

gx <- exprs(bgCor_lumiBatch)

neg <- raw_sample_data$neg_max

det <- sweep(gx, 1, neg,">")

probe_detection <- rowSums(det)
sample_detection <- colSums(det)
sample_detection

##-----------------------------------------------------------# 
## 4. Basic plots & QC on raw data 
##-----------------------------------------------------------# 



#-----------------------------------------------------------# 
# 7 probe & sample detection levels p<0.05 by pop, group
#-----------------------------------------------------------# 

###############################################
## probe
###############################################

probe_detection <- detectionCall(bgCor_lumiBatch, Th = 0.05, type = c('probe') )

probes_not_detected_in_any_sample <- probe_detection_mbcbNP==0 ## remove from all downstream analysis

n_samples <- dim( pData(eset_mbcbNP) )[1]

## add option for probe detection level

probes_detected_in_50percent_sample <- probe_detection_mbcbNP>=0.50*n_samples
probes_detected_in_80percent_sample <- probe_detection_mbcbNP>=0.80*n_samples
probes_detected_in_100percet_sample <- probe_detection_mbcbNP==n_samples

good_probes <- probes_detected_in_80percent_sample

##############################################################
## sample detection rate based on detected/expressed probes ## 
##############################################################

sample_detection_mbcbNPVSTQN <- detectionCall(eset_mbcbNP[probes_detected_in_80percent_sample , ], Th = 0.05, type = c('sample') )

n_probes <- sum(probes_detected_in_80percent_sample==TRUE)

good_samples <- sample_detection_mbcbNP>=0.80*n_probes

## subset to good samples and good probes

eset_mbcbNP_good <- eset_mbcbNP[ good_probes , good_samples ]


#-----------------------------------------------------------# 
# 8. vst transform removing poorChips
#-----------------------------------------------------------# 

## add option for log2 or vst transform

eset_mbcbNPVST  <- lumiT(eset_mbcbNP_good, ifPlot = FALSE, method="vst" )
eset_mbcbNPLOG2 <- lumiT(eset_mbcbNP_good, ifPlot = FALSE, method="log2" ) ###method = c("vst", 'log2', 'cubicRoot'))


#-----------------------------------------------------------# 
# 9. quantile normalise
#-----------------------------------------------------------# 

# add option for qn or normalise

eset_mbcbNPLOG2TQN <- lumiN(eset_mbcbNPLOG2, method="quantile")
eset_mbcnNPLOG2RSN <- lumiN(eset_mbcbNPLOG2, method="rsn")

##
eset_mbcbNPVSTQN <- lumiN(eset_mbcbNPVST, method="quantile")
eset_mbcnNPVSTRSN <- lumiN(eset_mbcbNPVST, method="rsn")



#  method: five different between chips normalization methods
#          ("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant")

system( paste(" chmod 776 ",out_dir,"/",project_name),sep="""

sessionInfo()
#-----------------------------------------------------------# 
# 10. basic plots 
#-----------------------------------------------------------# 

#-----------------------------------------------------------# 
# 11. Analysis ready data
#-----------------------------------------------------------# 
 
LIMMA/SAM/basic stats
WGCNA
ML
RF














