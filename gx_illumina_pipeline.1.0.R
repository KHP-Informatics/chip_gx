######################################################################################################################################################################################
## Authors : Stephen Newhouse, Amos Forlarin
## Email : stephen..j.newhouse@gmail.comp
## Date : 15 Nov 2013
## Verison: 1.00
######################################################################################################################################################################################
rm(list=ls())
data_dir <- getwd()
cat(" Data Dir ",data_dir,"\r","\n")
cat(" Set Working Dir to ", data_dir,"\r","\n")
setwd(data_dir)
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



#########################
## 1: Project settings ##
#########################

cat(" Getting Project settings and file names","\r","\n")

gs_report <- "Sample_and_Control_Probe_Profile_FinalReport.txt"
gs_probe <- "Group_Probe_Profile_Final_Report.txt" # genomestudio report
gs_sample <- "sample_table_Final_Report.txt" # genomestudio report
gs_control <- "control_probe_profile_Final_Report.txt" # genomestudio report

pheno_file <- "NULL" # FILE NAME OR "NULL"study specific information 
project_name <- "GAP" # GAP
out_dir <- paste(project_name,"_lumi_processing" ,sep="") # output dir for lumi process and plots

probe_det <- 80 # 50, 80, 90,100
sample_det <- 80 # 50,80,90,100

qc_plots <- 1 # 1 or 0
int_qc <- 1 # 1 or 0
tech_qc_sd_thrs <- 2 #
sex_check <- 1 # 1 or 0
iac_check <- 1 # 1 or 0
iac_sd_thrs <- 2 #

norm_method <- "rsn" # quantile, rsn
transform_method <- "vst" # log2, vst

cat(" Project Name ", project_name ,"\r","\n")
cat(" Merged Final Report ", project_name ,"\r","\n")
cat(" Probe Final Report ", gs_probe ,"\r","\n")
cat(" Sample Table Final Report ", gs_sample ,"\r","\n")
cat(" Control Probe Final Report ", gs_control ,"\r","\n")

cat(" Probe Detection Call Rate ", probe_det ,"\r","\n")
cat(" Sample Detection Call Rate ", sample_det ,"\r","\n")

cat(" Transformation Method ",transform_method,"\r","\n")
cat(" Normalisation Method ",norm_method,"\r","\n")

## groups data

cat(" Making processing directory. mkdir ./",out_dir,"\r","\n")

system( paste(" mkdir ./",out_dir,sep="") )

######################################################################################################################################################################################

cat(" Reading Sample Table Final Report generated in Genomestudio ", gs_sample ,"\r","\n")

raw_sample_data <- read.table(gs_sample ,skip=8,as.is=T,fill=T,head=T,sep="\t")

raw_n_samples <- dim(raw_sample_data)[1]

cat(" Number of samples = ",raw_n_samples,"\r","\n")


######################################################################################################################################################################################

############################################################
## Create Raw LumiBatch object                            ##
############################################################

cat(" Readging Genomestudio Final Report file ","\r","\n")
cat(" Calling lumiR() on ", gs_report,"\r","\n")

raw_lumiBatch <- lumiR(gs_report, lib.mapping="lumiHumanIDMapping",
annotationColumn = c('PROBE_ID','CHROMOSOME','SYMBOL','DEFINITION','ACCESSION','ENTREZ_GENE_ID','PROBE_TYPE','PROBE_START','PROBE_SEQUENCE','PROBE_CHR_ORIENTATION','PROBE_COORDINATES',
'CHROMOSOME','TRANSCRIPT','ILMN_GENE','REFSEQ_ID','UNIGENE_ID','SYMBOL','PROTEIN_PRODUCT'))

cat( print(raw_lumiBatch) )

## add sample data/information from genome studio report or external
## this should be the sample report from genomestudio plus any useful demographic data and phenotype data ie age,sex,case/control status, tissue type, phenotype etc etc 

cat(" adding Sample Report to pData() slot","\r","\n")

pData(raw_lumiBatch) <- raw_sample_data;

cat(" getting probe annotaions from fData() slot","\r","\n")

probe_annotations <- fData(raw_lumiBatch)

cat(" saving probe annotaions as",paste(out_dir,"/",project_name,".raw_lumiBatch.probe_annotations.RData",sep=""),"\r","\n")

save(probe_annotations,file=paste(out_dir,"/",project_name,".raw_lumiBatch.probe_annotations.RData",sep="") )

write.table(probe_annotations,file=paste(out_dir,"/",project_name,".raw_lumiBatch.probe_annotations.txt",sep=""),sep="\t",quote=F,row.names=F )

## Save raw_lumiBatch

cat(" Saving raw_lumiBatch as ", paste(out_dir,"/",project_name,".raw_lumiBatch.RData",sep=""), "\r","\n")

save(raw_lumiBatch , file=paste(out_dir,"/",project_name,".raw_lumiBatch.RData",sep="") )

## This is raw data pre-qc

######################################################################################################################################################################################



##-----------------------------------------------------------# 
## 4. Basic plots & QC on raw data 
##-----------------------------------------------------------# 

## plots on raw data pre mbcb ##

# boxplots

# cluster

# iac plots

# pca plots 

# xist plots - determine sex

# sample intensity plot 

## pca???


## QC in raw data pre ##



#-----------------------------------------------------------# 
# 5. background correct  basic plots and qc
#-----------------------------------------------------------# 

cat(" Start background correction ", "\r","\n")

sig <-  paste(out_dir,"/",project_name,".raw_lumiBatch_exprsSignal.txt",sep="")
neg <- paste(out_dir,"/",project_name,".raw_lumiBatch_negative.control.txt",sep="")

## get negative control bead data
cat(" get negative control bead data","\r","\n")

negativeControl <- getControlData(raw_lumiBatch)
negativeControl <- subset(negativeControl, negativeControl$controlType=="NEGATIVE")
negativeControl <- negativeControl[,c(3:dim(negativeControl)[2])]
negativeControl <- apply(negativeControl,2, negBeadOutlierRepMean)

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
cat(" run mbcb ","\r","\n")

gx_mbcb_NP <- mbcb.correct(signal,negCon,npBool=TRUE, isRawBead=FALSE)

##names(gx_mbcb_NP) # [1] "RMA"      "NP"       "MLE"      "Bayes"    "GMLE"     "AvgNP" "AvgRMA"   "AvgMLE"   "AvgBayes" "AvgGMLE" 

gx_mbcb <- gx_mbcb_NP$NP

cat(" mbcb complete ","\r","\n")

## replace names with original as R adds "X" to numbers

cat(" replace names with original sampleNames(raw_lumiBatch) as R adds X to numbers ","\r","\n")

colnames(gx_mbcb) <- sampleNames(raw_lumiBatch)

##c heck gx_mbcb[1:10,1:10]

## make new eset for bk corrected data

cat(" Creating Background Corrected Data set: bgCor_lumiBatch ","\r","\n")

bgCor_lumiBatch <- raw_lumiBatch

## replace old exprs data with new mbcb bk corrected data
cat(" replace old exprs data with new mbcb Background corrected data ","\r","\n")

exprs(bgCor_lumiBatch) <- as.matrix(gx_mbcb)

cat(" Saveing Background Corrected Data set: ",paste(out_dir,"/",project_name,".bgCor_lumiBatch.RData",sep=""),"\r","\n")

save(bgCor_lumiBatch, file=paste(out_dir,"/",project_name,".bgCor_lumiBatch.RData",sep="")  , compress=T)

##

probeDection <- function(gx,neg) {


}

max_neg <- max(neg)
probe_det <- ifelse( > max_neg,1,0)

probe_expressed <- apply(  )



negativeControl


#-----------------------------------------------------------# 
# 6. Sample qc - low_int_chips, iac outliers
#-----------------------------------------------------------# 


low_int_chips <- scale(pData(eset_mbcbNP)$Signal.Average)<(-2)  ## good measure of bad chips, can use other tech metircs from Genomestudio. can skip this if already done in genomestudio

datExprs0 <- exprs(eset_mbcbNP) 

iac_outliers_00 <-  outlierSamples(datExprs0) ## optional :- but worth one round of outlier detection

poorChips <- c(names(iac_outliers_00$samples_to_remove), sampleNames(eset_mbcbN)[low_int_chips] )

save(poorChips, file="poorChips.RData")

## remove poor samples/chips
eset_old <- eset_mbcbNP

eset_mbcbNP <- eset_old[  , (sampleNames(eset_old) %in% poorChips)==FALSE ]

rm("eset_old")

save(eset_mbcbNP, file="eset_mbcbNP.RData", compress=T)

#-----------------------------------------------------------# 
# 7 probe & sample detection levels p<0.05 by pop, group
#-----------------------------------------------------------# 

###############################################
## probe
###############################################

probe_detection_mbcbNP <- detectionCall(eset_mbcbNP, Th = 0.05, type = c('probe') )

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


#-----------------------------------------------------------# 
# 10. basic plots 
#-----------------------------------------------------------# 

#-----------------------------------------------------------# 
# 11. Analysis ready data
#-----------------------------------------------------------# 




#-----------------------------------------------------------# 
# 12. Gene Expression Module for (limma type diff ex)
#-----------------------------------------------------------# 
#if(gene_ex_check_run)
{
   # source(./gx_diff-ex_mod.R)
    # In progress Steve et al., basic pairwise factorial is done, but not tested
    # will add a option to supply custom design matrix
}



LIMMA/SAM/basic stats
WGCNA
ML
RF














