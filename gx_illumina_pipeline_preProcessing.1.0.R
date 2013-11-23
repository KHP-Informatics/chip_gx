#!/share/bin/Rscript

rm(list=ls())

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
##
## OPTIONAL BUT HIGHLY RECOMMENDED
## 1) Technical information file with batch processing and lab realted data eg Sample.ID,RIN,RNA_YIELD,DATE_CHIP_RUN,DATE_RNA_EXTRACTED
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

sessionInfo()

cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
## some functins
negBeadOutlierRepMean <- function(x) { 
		z_out_samp <- abs(  as.numeric( scale(x) )  ) > 2
		mean_pop <- mean(x[z_out_samp==FALSE])
		sd_pop <- sd(x[z_out_samp==FALSE])
		new_x <- ifelse( abs(  as.numeric( scale(x) )  ) > 2, mean_pop, x ) 
		return(new_x)
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

options(stringsAsFactors = FALSE)

##  args <- commandArgs(trailingOnly=TRUE);
##  config_file <- args[1];

config_file <- "project.config"

project_settings <- read.table(config_file, head=T,sep="\t",as.is=T,fill=TRUE)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

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

cat(" Phenotype Data Provided. Mergeing with Sample Table Final Report generated in Genomestudio ", gs_sample ,"\r","\n")

pheno_dat <- read.table(pheno_file, as.is=T,fill=T,head=T,sep="\t")

has_pheno_cols <- c("Sample.ID","SEX","GROUPS","TISSUE","PHENOTYPE","Study_ID") %in% names(pheno_dat);

missing_pheno_cols <- "FALSE" %in% has_pheno_cols 

if(missing_pheno_cols == "TRUE") stop(" YOU ARE MISSING ESSENTIAL SAMPLE INFORMATION! MAKE SURE YOUR PHENO_FILE HAS:- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID !!!")

raw_sample_data <- merge(raw_sample_data , pheno_dat , by.x="Sample.ID",by.y="Sample.ID",sort=FALSE,all.x=TRUE)  

rownames(raw_sample_data) <- raw_sample_data$Sample.ID

write.table(raw_sample_data,file=paste(out_dir,"/",project_name,".sample.pheno.txt",sep=""),sep="\t",quote=F,row.names=F )

cat(" Number of samples = ",raw_n_samples,"\r","\n")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

############################################################
## Create Raw LumiBatch object                            ##
############################################################

###############################################################
## Readging Genomestudio Final Report file & Calling lumiR() ##
###############################################################

 ## load(paste(out_dir,"/",project_name,".eset_raw.probe_annotations.RData",sep=""))

cat(" Reading Genomestudio Final Report file :",gs_report,"\r","\n")
cat(" Calling lumiR() on ", gs_report,"\r","\n")

eset_raw <- lumiR(gs_report, lib.mapping="lumiHumanIDMapping",
annotationColumn = c('PROBE_ID','CHROMOSOME','SYMBOL','DEFINITION','ACCESSION','ENTREZ_GENE_ID','PROBE_TYPE','PROBE_START','PROBE_SEQUENCE','PROBE_CHR_ORIENTATION','PROBE_COORDINATES',
'CHROMOSOME','TRANSCRIPT','ILMN_GENE','REFSEQ_ID','UNIGENE_ID','SYMBOL','PROTEIN_PRODUCT'))

## 9-12.4GB RAM!!!!

cat( print(eset_raw) )

cat(" Getting Chip Info ","\r","\n")

chip <- getChipInfo(eset_raw)

cat(" chip version: ", chip$chipVersion, "\r","\n")
cat(" chip species: ", chip$species, "\r","\n")
cat(" chip IDType: ", chip$IDType, "\r","\n")
cat(" chip chipProbeNumber: ", chip$chipProbeNumber, "\r","\n")
cat(" chip inputProbeNumber: ", chip$inputProbeNumber, "\r","\n")
cat(" chip matchedProbeNumber: ", chip$matchedProbeNumber, "\r","\n")
cat(" Number of samples = ",raw_n_samples,"\r","\n")

##########################################
## adding Sample Report to pData() slot ##
##########################################

cat(" adding Sample Report to pData() slot","\r","\n")

pData(eset_raw) <- raw_sample_data;

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

fData(eset_raw) <- probe_annotations

system(paste(" rm ",anno_table,".tmp",sep="") )

} else { 

cat(" getting probe annotations from fData() slot ","\r","\n")

probe_annotations <- fData(eset_raw) 

probe_annotations$nuID <- rownames(probe_annotations)

}

cat(" saving probe annotaions as",paste(out_dir,"/",project_name,".eset_raw.probe_annotations.RData",sep=""),"\r","\n")

save(probe_annotations,file=paste(out_dir,"/",project_name,".eset_raw.probe_annotations.RData",sep="") )

write.table(probe_annotations,file=paste(out_dir,"/",project_name,".eset_raw.probe_annotations.txt",sep=""),sep="\t",quote=F,row.names=F )

##########################
## Save eset_raw #########
##########################

cat(" Saving eset_raw as ", paste(out_dir,"/",project_name,".eset_raw.RData",sep=""), "\r","\n")

save(eset_raw , file=paste(out_dir,"/",project_name,".eset_raw.RData",sep="") )## This is raw data pre-qc

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## ADD TECH INFO OPTION OF PROCESSING DATES, RIN AND YEILDS. If true do pca col be these factors - 10,20,30,40,50...percentiles for quant data

####################
## PLOTS RAW DATA ##
####################

eset <- eset_raw

qc_plots_pdf <- paste(out_dir,"/",project_name,".eset_raw.qc.plots.pdf",sep="")

chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))

cat(" Start basic QC plots of eset_raw. Saved to: ",qc_plots_pdf,"\r","\n")

datExprs0 <- exprs(eset) 

cat(" calculate pca for plots ")

pca_raw <- prcomp(t(datExprs0))$x

cat(" plotting... ","\r","\n")

pdf(raw_qc_plots_pdf, width=16.5,height=11.7) ## A3 
	plot(eset, what='boxplot', main=paste(project_name," eset_raw Boxplot",sep=""), col=chip_col )
	plot(eset, what='density', main=paste(project_name," eset_raw Density",sep="") )
	plot(eset, what='cv',main=paste(project_name," eset_raw density plot of coefficient of varience",sep="")  )
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
	tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()

cat(" finished basic QC plots of eset_raw. Saved to: ",raw_qc_plots_pdf,"\r","\n")

rm("datExprs0","tmp_iac","pca_raw");

gc()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

######################################################################
## IAC QC                                                           ##
## ## optional :- but worth one round of outlier detection          ##
######################################################################

if(iac_check==1) {

cat(" Starting inter-array correlation quality control step ","\n","\r")

datExprs0 <- exprs(eset_raw)

pdf(file=paste(out_dir,"/",project_name,".eset_raw.iac_outliers.pdf",sep=""), width=8,height=6 )

iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection

plot(iac_outlierSamples$samplemeanIAC)

dev.off()

## meanIAC 
cat(" meanIAC: ", iac_outlierSamples$meanIAC , "\r","\n")

## add iac sd to pData
cat(" addingd IAC sd to pData() slot of eset_raw ", "\r","\n")

pData(eset_raw)$iac_sd <- iac_outlierSamples$numbersd

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

#########################
##  remove iac samples ##
#########################
sel_samp <- (sampleNames(eset_raw) %in% names(iac_outliers) )==FALSE

sel_samp_names <- sampleNames(eset_raw)[sel_samp]

cat(" removing iac outliers from eset_raw","\r","\n")

eset_raw_iac <- eset_raw[,sel_samp]

cat(" update Control Probe data slot..subsetting to remaining samples ", "\r","\n")

## update control data slot
update_controlData <- getControlData(eset_raw_iac)

update_controlData  <- update_controlData[,c("controlType","ProbeID",sel_samp_names)]

eset_raw_iac <- addControlData2lumi(update_controlData , eset_raw_iac)

cat(" saving updated eset_raw post IAC to: ",paste(out_dir,"/",project_name,".eset_raw_postIAC.RData",sep=""), "\r","\n" )

save(eset_raw_iac , file=paste(out_dir,"/",project_name,".eset_raw_postIAC.RData",sep="") )

cat(" Number of samples after iac outlier removal: ", length(sampleNames(eset_raw_iac)), "\r","\n" )

} else {

cat(" Skipping inter-array correlation quality control step ","\n","\r")

}

#############################
## PLOTS POST IAC RAW DATA ##
#############################

eset <- eset_raw_iac

qc_plots_pdf <- paste(out_dir,"/",project_name,".eset_raw.post.iac.qc.plots.pdf",sep="")

chip_col <- labels2colors( as.character(pData(eset_raw_iac)$Sentrix.Barcode))

cat(" Start basic QC plots of post IAC eset_raw_iac. Saved to: ",qc_plots_pdf,"\r","\n")

datExprs0 <- exprs(eset)

cat(" calculate pca for plots post IAC ")

pca_raw <- prcomp(t(datExprs0))$x

cat(" plotting post IAC... ","\r","\n")

pdf(raw_qc_plots_pdf, width=16.5,height=11.7) ## A3 
	plot(eset_raw, what='boxplot', main=paste(project_name," eset_raw_iac Boxplot",sep=""), col=chip_col )
	plot(eset_raw, what='density', main=paste(project_name," eset_raw_iac Density",sep="") )
	plot(eset_raw, what='cv',main=paste(project_name," eset_raw_iac density plot of coefficient of varience",sep="")  )
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
	tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()

cat(" finished basic QC plots of post IAC eset_raw_iac. Saved to: ",raw_qc_plots_pdf,"\r","\n")

rm("datExprs0","tmp_iac","pca_raw")
gc()

#########################################
## Check IAC after 1st round of iac qc ##
#########################################

cat(" Checking IAC after 1st round of IAC QC ","\n","\n")

preIAC <- iac_outlierSamples$meanIAC

rm("iac_outlierSamples")

datExprs0 <- exprs(eset_raw_iac) 

pdf(file=paste(out_dir,"/",project_name,".eset_raw.post.iac_outliers.pdf",sep=""), width=8,height=6 )

iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection

plot(iac_outlierSamples$samplemeanIAC)

dev.off()

postIAC <- iac_outlierSamples$meanIAC

cat(" meanIAC pre  IAC QC: ", preIAC, "\r","\n") 
cat(" meanIAC post IAC QC: ", postIAC,"\r","\n")

rm("datExprs0")
gc()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

##############################################################
## background correct using mbcb                            ##
##############################################################

cat(" Start background correction ", "\r","\n")

sig <- paste(out_dir,"/",project_name,".eset_raw_exprsSignal.txt",sep="")
neg <- paste(out_dir,"/",project_name,".eset_raw_negative.control.txt",sep="")

## get negative control bead data
cat(" get negative control bead data","\r","\n")

negativeControl <- getControlData(eset_raw_iac)
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

expressionSignal <- exprs(eset_raw_iac)

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

cat(" saveing raw mbcb data to ",paste(out_dir,"/",project_name,".eset_raw.mbcb.RData",sep=""), "\r","\n")

save(gx_mbcb , file=paste(out_dir,"/",project_name,".eset_raw.mbcb.correct.RData",sep=""))

if(mbcb_method=="NP") { gx_mbcb <- gx_mbcb$NP }

if(mbcb_method=="MLE") { gx_mbcb <- gx_mbcb$MLE }

cat(" mbcb complete ","\r","\n")

## replace names with original as R adds "X" to numbers

cat(" replace names with original sampleNames(eset_raw) as R adds X to numbers ","\r","\n")

colnames(gx_mbcb) <- sampleNames(eset_raw_iac)

## make new eset for bk corrected data

cat(" Creating Background Corrected Data set: bgCor_lumiBatch ","\r","\n")

eset_bkCor <- eset_raw_iac

## replace old exprs data with new mbcb bk corrected data
cat(" replace old exprs data with new mbcb Background corrected data ","\r","\n")

exprs(eset_bkCor) <- as.matrix(gx_mbcb)

cat(" Saveing Background Corrected Data set: ",paste(out_dir,"/",project_name,".eset_bkCor.RData",sep=""),"\r","\n")

save(eset_bkCor, file=paste(out_dir,"/",project_name,".eset_bkCor.RData",sep="")  , compress=T)

#############################
## PLOTS POST bkCor        ##
#############################

bkCor_qc_plots_pdf <- paste(out_dir,"/",project_name,".eset_bkCor.post.bkCor.qc.plots.pdf",sep="")

chip_col <- labels2colors( as.character(pData(eset_bkCor)$Sentrix.Barcode))

cat(" Start basic QC plots of post IAC eset_raw. Saved to: ",bkCor_qc_plots_pdf,"\r","\n")

datExprs0 <- exprs(eset_bkCor) 

cat(" calculate pca for plots post bkCor ")

pca_raw <- prcomp(t(datExprs0))$x

cat(" plotting post IAC... ","\r","\n")

pdf(bkCor_qc_plots_pdf, width=16.5,height=11.7) ## A3
	plot(eset_bkCor, what='boxplot', main=paste(project_name," eset_bkCor Boxplot",sep=""), col=chip_col )
	plot(eset_bkCor, what='density', main=paste(project_name," eset_bkCor Density",sep="") )
	plot(eset_bkCor, what='cv',main=paste(project_name," eset_bkCor density plot of coefficient of varience",sep="")  )
	scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color=chip_col)
	tmp_iac <- outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=TRUE)
dev.off()

cat(" finished basic QC plots of post bkCor on eset_raw. Saved to: ",bkCor_qc_plots_pdf,"\r","\n")

rm("datExprs0","tmp_iac","pca_raw")
gc()


#########################################
## Check IAC after 1st round of iac qc ##
#########################################

cat(" Checking IAC after bkCor ","\n","\n")

rm("iac_outlierSamples")

datExprs0 <- exprs(eset_bkCor) 

pdf(file=paste(out_dir,"/",project_name,".eset_bkCor.iac_outliers.pdf",sep=""), width=8,height=6 )

iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection

plot(iac_outlierSamples$samplemeanIAC)

dev.off()

postbkCorIAC <- iac_outlierSamples$meanIAC

cat(" meanIAC pre  IAC QC: ", preIAC, "\r","\n") 
cat(" meanIAC post IAC QC: ", postIAC,"\r","\n")
cat(" meanIAC post bkCor: ",  postbkCorIAC,"\r","\n")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
##################################
## Sex Check based on XIST GENE ##
##################################

cat(" Checking Gender based on XIST gene expression","\r","\n")

##load(paste(out_dir,"/",project_name,".eset_raw_postIAC.RData",sep=""))
##load(paste(out_dir,"/",project_name,".eset_bkCor.RData",sep=""))

## get neg control data from eset_bkCor
negativeControl <- getControlData(eset_bkCor)
negativeControl <- subset(negativeControl, negativeControl$controlType=="NEGATIVE")
negativeControl <- negativeControl[,c(3:dim(negativeControl)[2])]

## get neg control info mean,sd,max etc
neg_max <- apply(negativeControl,2,max)
neg_sd <- apply(negativeControl,2,sd)
neg_mean <- apply(negativeControl,2,mean)
neg_2sd <- neg_mean + 2*neg_sd

## get XIST gene postion
xist_raw <- fData(eset_raw_iac)$ILMN_GENE=="XIST";

xist_bkCor <- fData(eset_bkCor)$ILMN_GENE=="XIST";

## get XIST gene expression signal
xist_gx_raw  <- exprs(eset_raw_iac[xist_raw, ]  ) 
xist_gx_raw <- as.data.frame( t(xist_gx_raw)); 

xist_gx_bkCor <- exprs(eset_bkCor[xist_bkCor , ]  ) 
xist_gx_bkCor <- as.data.frame( t(xist_gx_bkCor)); 

## cobine raw and bkCro gx data
xist_gx <- cbind(xist_gx_raw,xist_gx_bkCor)

colnames(xist_gx) <- c("raw_XIST","bkCor_XIST")

xist_gx$neg_2sd <- neg_2sd

xist_gx$neg_max <- neg_max

xist_gx$neg_mean <- neg_mean

xist_gx$neg_sd <- neg_sd

## gender based on XIST expression 1=female , 0 = male
xist_gx$XIST_Gender_max <- ifelse(xist_gx$bkCor_XIST > xist_gx$neg_max,1,0)

xist_gx$XIST_Gender_2sd <- ifelse(xist_gx$bkCor_XIST > xist_gx$neg_2sd,1,0)  ### THIS WORKS!! IF GENE IS > 2SD_NEG_BEAD THEN ITS EXPRESSED/DETECTED!

xist_gx$XIST_z <-  ( xist_gx$bkCor_XIST - xist_gx$neg_mean ) / xist_gx$neg_sd

xist_gx$Sample.ID <- rownames(xist_gx)

head(xist_gx)

xist_gx$xist_gender <- ifelse(xist_gx$bkCor_XIST > xist_gx$neg_2sd,"Female","Male")

head(xist_gx)

pdf(file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.pdf",sep=""),height=8,width-11)
plot(xist_gx$bkCor_XIST, type="p", col=2+xist_gx$XIST_Gender_2sd, pch=19 , cex=0.6);points(xist_gx$neg_2sd,col="blue", type="l")
dev.off()

pData(eset_bkCor)$xist_gender <- xist_gx$xist_gender 

save(xist_gx, file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.RData",sep=""))

write.table(xist_gx, file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.XIST.csv",sep=""),sep=",",row.names=F)

## check sex 

pData(eset_bkCor)$gender_missmatch <- ifelse( pData(eset_bkCor)$xist_gender == pData(eset_bkCor)$SEX, "PASS","GENDER_MISSMATCH" ) 

save(eset_bkCor, file=paste(out_dir,"/",project_name,".eset_bkCor.RData",sep=""))

n_gender_fails <- sum( pData(eset_bkCor)$gender_missmatch=="GENDER_MISSMATCH" )

if ( n_gender_fails > 0 ) {
    cat(" WARNING: Youn have GENDER_MISSMATCH samples!!!!!!! N=[",n_gender_fails,"]","\r","\n")
} else {
  cat(" Congratulations! \n All your Males are Male and Females are Female. \n You have NO GENDER_MISSMATCH samples!!!", "\r","\n")
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------#
############################
## PROBE ABOVE BACKGROUND ##
############################

cat(" Calculating Probde Detection rates. \n Probe is seen as Detected if it has signal intensity >= 2SD of the mean intensity of the negative control beads","\r","\n")

## get expression matrix
gx <- exprs(eset_bkCor)

## get negative bead ranges mean or max or >2SD mean  of neg beads

neg_2sd <- neg_mean + 2*neg_sd

## sweep through gx matrix to id probes > 2SD Mean of 
det <- sweep(gx, 1, neg_2sd,">=")

## probe_detection counts
probe_detection <- rowSums(det)

## probes_not_detected_in_any_sample
n_samples <- dim(gx)[2]

## probe annotations 
probes_not_detected_in_any_sample <- probe_detection==0
probes_detected_in_80_sample <- probe_detection>=n_samples*0.80
probes_detected_in_all_sample <- probe_detection==n_samples 

probe_annotations_0_detected <- fData(eset_bkCor[probes_not_detected_in_any_sample,])
probe_annotations_80_detected  <- fData(eset_bkCor[probes_detected_in_80_sample ,])
probe_annotations_100_detected  <- fData(eset_bkCor[probes_detected_in_all_sample,])

fData(eset_bkCor)$n_detected <- probe_detection 
fData(eset_bkCor)$n_detected_call_rate <- round( probe_detection/n_samples ,3)


## sample_detection counts
n_probes <- dim(eset_bkCor)[1]

sample_detection <- colSums(det)

pData(eset_bkCor)$n_probes_detected <- sample_detection

pData(eset_bkCor)$n_probes_detected_call_rate <- round( sample_detection/n_probes ,3)

plot(pData(eset_bkCor)$n_probes_detected_call_rate)

plot(pData(eset_bkCor)$n_probes_detected ,pData(eset_bkCor)$Detected.Genes..0.01.)

save(eset_bkCor, file=paste(out_dir,"/",project_name,".eset_bkCor.RData",sep=""))

#--------------------------------------------------------------------------------------------------------------------------------------------------------#
## get group information from pData() slot 

group_names <- unique(pData(eset_bkCor)$GROUPS);

groups  <- pData(eset_bkCor)$GROUPS

n_groups <- length(group_names)

##
gx <- exprs(eset_bkCor)

## loop through each group and id probes > 2SD mean neg beads in X% of samples/group 

for(n in group_names ) {

cat(" Finding probes in ",probe_det," of sample group [",n,"] with signal intensity >= 2SD mean intensity of the negative control beads ","\r","\n")

group_label <- paste(n)

sel_samples <- pData(eset_bkCor)$GROUPS==n;

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
## Y CHROM EXPRESSION IN XIST MALES

cat(" Y Chromosome probe detection based on XIST males","\r","\n")

xist_males <- pData(eset_bkCor)$xist_gender=="Male"

gx_y <- exprs(eset_bkCor[fData(eset_bkCor)$CHR=="Y",])

detection_matrix_y <- sweep( gx_y[,xist_males],1, neg_2sd[xist_males ],">=")

y_probe_detection <- rowSums(detection_matrix_y) >= probe_det * sum(xist_males==TRUE)

y_probe_detection_nuID <- rownames( gx_y[y_probe_detection, ])

y_det_probes <- as.data.frame(y_probe_detection_nuID)

colnames(y_det_probes) <- c("nuID")

write.table( y_det_probes ,file=paste(out_dir,"/",project_name,".GROUP.Y.detected_probes_nuID.txt",sep=""),row.names=FALSE,quote=FALSE,col.names=FALSE)

#--------------------------------------------------------------------------------------------------------------------------------------------------------#
#######################

cat(" writing final good probe list to ",paste(out_dir,"/",project_name,".detected_probes_nuID_final.txt",sep=""),"\r","\n")

system(paste("cat ",out_dir,"/",project_name,"****.detected_probes_nuID.txt | sort | uniq >> ",out_dir,"/",project_name,".detected_probes_nuID_final.txt", sep="") )

good_probes <- read.table(file=paste(out_dir,"/",project_name,".detected_probes_nuID_final.txt",sep=""),head=FALSE)

good_probes <- paste(good_probes[,1])

n_good_probes <- length(good_probes)

cat(" Total number of good probes = ",n_good_probes,"\n","\r" )

good_probes_annotation <- fData(eset_bkCor[good_probes,])

good_probes_annotation$raw_mean <- apply( exprs(eset_bkCor[good_probes,]), 1,mean) 
good_probes_annotation$raw_sd <- apply( exprs(eset_bkCor[good_probes,]), 1,sd) 
good_probes_annotation$raw_var <- apply( exprs(eset_bkCor[good_probes,]), 1,var) 
good_probes_annotation$raw_min <- apply( exprs(eset_bkCor[good_probes,]), 1,min) 
good_probes_annotation$raw_max <- apply( exprs(eset_bkCor[good_probes,]), 1,max) 

head(good_probes_annotation)

cat(" saving good probe annotations to ",paste(out_dir,"/",project_name,".detected_probes_nuID_final.***",sep=""),"\r","\n")

save(good_probes_annotation, file=paste(out_dir,"/",project_name,".detected_probes_nuID_final.RData",sep="") )

write.table(good_probes_annotation, file=paste(out_dir,"/",project_name,".detected_probes_nuID_final.tsv",sep=""),quote=F,sep="\t",row.names=F )
#--------------------------------------------------------------------------------------------------------------------------------------------------------#


####################################################################
# 8. vst transform 
####################################################################

cat(" TRANSFORM DATA. Method =",transform_method,"\r","\n")

if(transform_method=="vst") {
  eset <- eset_bkCor
  eset_lumiT  <- lumiT(eset, ifPlot = FALSE, method="vst" )
}
if(transform_method=="log2") {
  eset <-  eset_bkCor
  eset_lumiT  <- lumiT(eset, ifPlot = FALSE, method="log2" )
}

## SAVE
save(eset_lumiT , file=paste(out_dir,"/",project_name,".eset_lumiT.RData",sep="") )

## PLOTS ##
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
cat(" finished basic QC plots of eset_lumiT Saved to: ",raw_qc_plots_pdf,"\r","\n")
rm("datExprs0","tmp_iac","pca_raw");
gc()
## iac
cat(" Checking IAC after eset_lumiT ","\n","\n")
rm("iac_outlierSamples")
datExprs0 <- exprs(eset_lumiT) 
pdf(file=paste(out_dir,"/",project_name,".eset_lumiT.iac_outliers.pdf",sep=""), width=8,height=6 )
iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection
plot(iac_outlierSamples$samplemeanIAC)
dev.off()
eset_lumiT_IAC <- iac_outlierSamples$meanIAC
cat(" meanIAC post VST: ",  eset_lumiT_IAC,"\r","\n")




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
## SAVE
save(eset_lumiN , file=paste(out_dir,"/",project_name,".eset_lumiN.RData",sep="") )

## PLOTS ##
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
cat(" finished basic QC plots of eset_lumiN Saved to: ",raw_qc_plots_pdf,"\r","\n")
rm("datExprs0","tmp_iac","pca_raw");
gc()
## iac
cat(" Checking IAC after eset_lumiN ","\n","\n")
rm("iac_outlierSamples")
datExprs0 <- exprs(eset_lumiN) 
pdf(file=paste(out_dir,"/",project_name,".eset_lumiN.iac_outliers.pdf",sep=""), width=8,height=6 )
iac_outlierSamples <-  outlierSamples(datExprs0,thresh=iac_sd_thrs, showplots=T) ## optional :- but worth one round of outlier detection
plot(iac_outlierSamples$samplemeanIAC)
dev.off()
eset_lumiN_IAC <- iac_outlierSamples$meanIAC
cat(" meanIAC post normalise: ",  eset_lumiN_IAC,"\r","\n")

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

## IAC on eset_lumiN

eset_final <- eset_lumiN[good_probes,]

save(eset_final , file=paste(out_dir,"/",project_name,".eset_final.RData",sep="") )



cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
cat(" END PRE_PROCESSING!" ,"\r","\n")
cat(" #--------------------------------------------------------------------------------------------------------------------------------------------------------# " ,"\r","\n")
##
system( paste(" chmod 776 ",out_dir,"/",project_name,"***",sep=""))

sessionInfo()

rm(list=ls())













