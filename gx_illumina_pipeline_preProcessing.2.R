#!/share/bin/Rscript
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
# begin
rm(list=ls())
options(stringsAsFactors = FALSE)

# setwd
setwd()

# load libs
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)
library(gplots)
library(limma)
library(vsn)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)

# load processing functions
path_to_scripts <- "" # path to gene expression processing scripts

source(paste(path_to_scripts,"/pre_process_gx.R",sep=""))

# read config file
# User to set a config file that stores the analysis settings and options
config_file <- "project.config"
project_settings <- read.table(config_file, head=TRUE, sep="\t")
# check project settings
t(project_settings)

# set project settings and I/O
# User is askes to manually provide options

# project directory
project_dir <- project_settings$project_dir ## /scratch/project/pipelines/DATA/expression/GAP_Expression
setwd(project_dir)

# project name
project_name <- project_settings$project_name ## "GAP

# output directory for lumi process and plots
out_dir <- paste(project_name,"_lumi_processing" ,sep="")

# make project pre-processing directory
make_dir_command <- paste(" if [ ! -e ./",out_dir," ]; then mkdir ./",out_dir,"; fi",sep="")
system( make_dir_command )

# genomestudio reports
gs_report <- project_settings$gs_report ## "Sample_and_Control_Probe_Profile_FinalReport.txt"
gs_probe <- project_settings$gs_probe ## "Group_Probe_Profile_Final_Report.txt" # genomestudio report
gs_sample <- project_settings$gs_sample ## "sample_table_Final_Report.txt" # genomestudio report
gs_control <- project_settings$gs_control ## "control_probe_profile_Final_Report.txt" # genomestudio report
anno_table <- project_settings$anno_table ## annotation table for project

# sample information
pheno_file <- project_settings$pheno_file ## FILE NAME must contain : Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID
tech_pheno_file <- project_settings$tech_pheno_file ## "sample_tech_info.txt"

#
probe_det <- 80
sample_det <- 80

sex_check <- 1
iac_check <- 1
iac_sd_thrs <- 2

# pre-processing methods

# Model based background correction method (MLE as default)
mbcb_method <- "MLE"

# Transformation
transform_method <- project_settings$transform_method ## "vst" # log2, vst or both

# Normalisation
norm_method <- project_settings$norm_method ## "rsn" # quantile, rsn, or both

# BEGIN PRE-PROCESSING

# read raw gene expression data



##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
## lumiB(x.lumi, method = c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy'), verbose = TRUE, ...)
## lumiN(x.lumi, method = c("quantile", "rsn", "vsn"), verbose = TRUE, ...)
## lumiT(x.lumi, method = c("vst", 'log2', 'cubicRoot'), ifPlot = FALSE, simpleOutput = TRUE, verbose = TRUE, ...)
## lumiExpresso(lumiBatch,bg.correct=TRUE, bgcorrect.param=list(method='bgAdjust'), variance.stabilize=TRUE, varianceStabilize.param = list(), normalize=TRUE, normalize.param=list(), QC.evaluation = TRUE, QC.param = list(), verbose = TRUE)
## bgcor_log_quantile <- lumiExpresso(bgcor, variance.stabilize=TRUE, varianceStabilize.param=list(method='log2'), normalize.param=list(method='quantile'), verbose=FALSE)
## bgcor_log_rsn      <- lumiExpresso(bgcor, variance.stabilize=TRUE, varianceStabilize.param=list(method='log2'), normalize.param=list(method='rsn'),      verbose=FALSE)
## bgcor_vst_quantile <- lumiExpresso(bgcor, variance.stabilize=TRUE, varianceStabilize.param=list(method='vst'),  normalize.param=list(method='quantile'), verbose=FALSE)
## bgcor_vst_rsn      <- lumiExpresso(bgcor, variance.stabilize=TRUE, varianceStabilize.param=list(method='vst'),  normalize.param=list(method='rsn'),      verbose=FALSE)


    Options for processing

	bg_rma
	bg_mb

	bg_vst_rsn*
	bg_vst_quantile*

	bg_log_rsn
	bg_log_quantile

	bg_vsn




## 1 RAW qc plots, write dat, sample networks, pc regressions
## 2 bg correct qc plots, write dat, sample networks, pc regressions
## 3 Probe Detection 2SD > mean neg probe
## 4 VST qc plots, write dat, sample networks, pc regressions
## 5 RSN qc plots, write dat, sample networks, pc regressions
## 6 Batch Adjustment qc plots, write dat, sample networks, pc regressions
## 7 Final Data good probes, good samples qc plots, write dat, sample networks, pc regressions


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











