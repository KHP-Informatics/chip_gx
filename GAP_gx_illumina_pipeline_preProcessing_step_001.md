


Microarry pre-processing workflow for Illumina BeadArray data
==============================================================
- Stephen Newhouse
- stephen.newhouse@kcl.ac.uk
Example workflow based on real data (GAP). 

# Project Directory and Files
Data files and directories used in this workflow

```r

# Data Directory
"/media/D/expression/GAP_Expression"

# Data Files Exported from GenomeStudio
control_probe_profile_Final_Report.txt
Group_Probe_Profile_Final_Report.txt
probe_annotation_Final_Report.txt
sample_table_Final_Report.txt

# Input for lumiR
Sample_and_Control_Probe_Profile_FinalReport.txt 

# User Provided Data Files
# Use NA to record missing values in batch_info.txt
# Use Unknown to record missing values in pheno_info.txt 

# REQUIRED COLUMNS (MATCH FORMAT SHOWN HERE ie SEX not Sex or Gender etc!):- "Sample.ID","SEX","GROUPS","TISSUE","PHENOTYPE","Study_ID"
pheno_info.txt

# REQUIRED COLUMNS:- "Sample.ID","RIN","RNA_YIELD" and any other related batch info eg dates or processing. 
batch_info.txt


# Naming Convensions
* All UPPERCASE
* SEX = MALE, FEMALE or UNKNOWN
* Missing data = NA for all data. The exceptioins are: SEX,GROUPS,PHENOTYPE, TISSUE. Use UNKNOWN

TO DO: ADD LINK TO Genomestudio SOP and scripts for making lumiR input

```


****

## setwd
This is the full path to the directory where all the raw genomestudio output is stored


```r
setwd("/media/D/expression/GAP_Expression")
```


## load libs

```r
# Load libraries
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)
```

```
## ==========================================================================
## *
## *  Package WGCNA 1.34 loaded.
## *
## *    Important note: It appears that your system supports multi-threading,
## *    but it is not enabled within WGCNA in R. 
## *    To allow multi-threading within WGCNA with all available cores, use 
## *
## *          allowWGCNAThreads()
## *
## *    within R. Use disableWGCNAThreads() to disable threading if necessary.
## *    Alternatively, set the following environment variable on your system:
## *
## *          ALLOW_WGCNA_THREADS=<number_of_processors>
## *
## *    for example 
## *
## *          ALLOW_WGCNA_THREADS=8
## *
## *    To set the environment variable in linux bash shell, type 
## *
## *           export ALLOW_WGCNA_THREADS=8
## *
## *     before running R. Other operating systems or shells will
## *     have a similar command to achieve the same aim.
## *
## ==========================================================================
```

```r
library(gplots)
library(limma)
library(vsn)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)
library(plyr)
library(ggplot2)
```


## load source file with some processing functions
email stephen.newhouse@kcl.ac.uk for code. This will all be on git soon

```r
# path to gene expression processing scripts
path_to_scripts <- "/media/D/sjnewhouse/GENE_EXPRESSION"
# load 'em
source(paste(path_to_scripts, "/sjnewhouse_misc_R.R", sep = ""))
ls()
```

```
##  [1] "basic_sampleNetwork"        "basic_sampleNetworkIterate"
##  [3] "bgcor_mbcb"                 "cv"                        
##  [5] "data_summary_plots"         "gx_qc_plots_lumi"          
##  [7] "gx_qc_plots_lumi_2"         "has_var_probe"             
##  [9] "has_var_probe2"             "max_probe"                 
## [11] "mean_probe"                 "min_probe"                 
## [13] "negBeadOutlierRepMean"      "outlierSamples"            
## [15] "outlierSamplesIterate"      "path_to_scripts"           
## [17] "preProcess"                 "quantfun"                  
## [19] "removeOutlierProbes"        "removeOutlierProbesIterate"
## [21] "removeSamples_eset_lumi"    "removeTooManyNAs"          
## [23] "runVSN"                     "sd_probe"                  
## [25] "shuffle_cols"               "shuffle_rows"              
## [27] "var_probe"                  "write_expression_files"    
## [29] "zero_var_probe"
```


## set project settings and I/O
User is asked to manually provide these options.  
This sets the working directoty, prjoect name, input and output files, along with qc options for transformation and normalisation methods.  
This project configuration data is written to *.csv file in your project directory.


```r
# project directory
project_dir <- "/media/D/expression/GAP_Expression"

# set working dir again
setwd(project_dir)

# project name
project_name <- "GAP"

# output directory for lumi process and plots
processing_date <- format(Sys.Date(), "%d_%m_%Y")
out_dir <- paste(project_dir, "/", project_name, "_lumi_processing_", processing_date, 
    sep = "")

# make project pre-processing directory
make_dir_command <- paste(" if [ ! -e ", out_dir, " ]; then mkdir ", out_dir, 
    "; fi", sep = "")

system(make_dir_command)

# genomestudio reports
gs_report <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt"
gs_probe <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/Group_Probe_Profile_Final_Report.txt"  # genomestudio report
gs_sample <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/sample_table_Final_Report.txt"  # genomestudio report
gs_control <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/control_probe_profile_Final_Report.txt"  # genomestudio report
anno_table <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/probe_annotation_Final_Report.txt"

# sample information FILE NAME must contain :
# Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID
pheno_file <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/pheno_info.txt"

# batch information
tech_pheno_file <- "/media/D/expression/GAP_Expression/final_reports_genomestudio/batch_info.txt"

# detection call rate threshold
probe_det <- 80
sample_det <- 80

# flag for gender and sampleNetwork
sex_check <- 1  # DO THIS!! I'm not providing an option to skip this
iac_check <- 1  # DO THIS!! I'm not providing an option to skip this
iac_sd_thrs <- 2  # 2 or 3

# Model based background correction method (MLE as default) All data should
# be background correceted. The recomended methods is MBCB (Model-based
# Background Correction for Beadarray) URL
# http://www.bioconductor.org/packages/release/bioc/html/MBCB.html
mbcb_method <- "MLE"

# Transformation method
transform_method <- "log2"  ## 'vst' # log2, vst or both

# Normalisation method
norm_method <- "rsn"  ## 'rsn' # quantile, rsn, or both

# Folks that done stuff
analyst_email <- "stephen.newhouse@kcl.ac.uk"
analyst_name <- "Stephen Newhouse"
lab_contact_email <- "charles.curtis@kcl.ac.uk"
lab_contact_name <- "Charle Curtis"


# write settings to file
project_settings <- data.frame(project_dir = project_dir, project_name = project_name, 
    out_dir = out_dir, gs_report = gs_report, gs_probe = gs_probe, gs_sample = gs_sample, 
    gs_control = gs_control, anno_table = anno_table, pheno_file = pheno_file, 
    tech_pheno_file = tech_pheno_file, probe_det = probe_det, sample_det = sample_det, 
    sex_check = sex_check, iac_check = iac_check, iac_sd_thrs = iac_sd_thrs, 
    mbcb_method = mbcb_method, transform_method = transform_method, norm_method = norm_method, 
    analyst_email = analyst_email, analyst_name = analyst_name, lab_contact_email = lab_contact_email, 
    lab_contact_name = lab_contact_name)

# some data wrangling
project_settings <- as.data.frame(t(project_settings))
colnames(project_settings) <- "Project_Setting"
project_settings$Project_Variable <- rownames(project_settings)
project_settings <- project_settings[, c("Project_Variable", "Project_Setting")]

# write table to out_dir
write.table(project_settings, file = paste(out_dir, "/", project_name, ".project_settings.csv", 
    sep = ""), row.names = FALSE, quote = FALSE, sep = ",")

# check settings
project_settings
```

```
##                    Project_Variable
## project_dir             project_dir
## project_name           project_name
## out_dir                     out_dir
## gs_report                 gs_report
## gs_probe                   gs_probe
## gs_sample                 gs_sample
## gs_control               gs_control
## anno_table               anno_table
## pheno_file               pheno_file
## tech_pheno_file     tech_pheno_file
## probe_det                 probe_det
## sample_det               sample_det
## sex_check                 sex_check
## iac_check                 iac_check
## iac_sd_thrs             iac_sd_thrs
## mbcb_method             mbcb_method
## transform_method   transform_method
## norm_method             norm_method
## analyst_email         analyst_email
## analyst_name           analyst_name
## lab_contact_email lab_contact_email
## lab_contact_name   lab_contact_name
##                                                                                                                  Project_Setting
## project_dir                                                                                   /media/D/expression/GAP_Expression
## project_name                                                                                                                 GAP
## out_dir                                                        /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014
## gs_report         /media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt
## gs_probe                      /media/D/expression/GAP_Expression/final_reports_genomestudio/Group_Probe_Profile_Final_Report.txt
## gs_sample                            /media/D/expression/GAP_Expression/final_reports_genomestudio/sample_table_Final_Report.txt
## gs_control                  /media/D/expression/GAP_Expression/final_reports_genomestudio/control_probe_profile_Final_Report.txt
## anno_table                       /media/D/expression/GAP_Expression/final_reports_genomestudio/probe_annotation_Final_Report.txt
## pheno_file                                          /media/D/expression/GAP_Expression/final_reports_genomestudio/pheno_info.txt
## tech_pheno_file                                     /media/D/expression/GAP_Expression/final_reports_genomestudio/batch_info.txt
## probe_det                                                                                                                     80
## sample_det                                                                                                                    80
## sex_check                                                                                                                      1
## iac_check                                                                                                                      1
## iac_sd_thrs                                                                                                                    2
## mbcb_method                                                                                                                  MLE
## transform_method                                                                                                            log2
## norm_method                                                                                                                  rsn
## analyst_email                                                                                         stephen.newhouse@kcl.ac.uk
## analyst_name                                                                                                    Stephen Newhouse
## lab_contact_email                                                                                       charles.curtis@kcl.ac.uk
## lab_contact_name                                                                                                   Charle Curtis
```

****

BEGIN PRE-PROCESSING
=====================

Raw Expression Set
-------------------

## 1. read raw gene expression data 

```r

# raw input This is the 1) Probe Profile, 2) Control Probe Profile and 3)
# Sample Table, Final Reports exported from GenomeStudio, all concatenated
if (is.na(gs_report)) stop(" WARNING!: YOU HAVENT PROVIDED ANY DATA TO READ")

# read raw gene expression data from genomestudio reports and create
# ExpressionSet

eset_raw <- lumiR(paste(gs_report), lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, 
    detectionTh = 0.01, convertNuID = TRUE, inputAnnotation = TRUE, annotationColumn = c("PROBE_ID", 
        "CHROMOSOME", "SYMBOL", "DEFINITION", "ACCESSION", "ENTREZ_GENE_ID", 
        "PROBE_TYPE", "PROBE_START", "PROBE_SEQUENCE", "PROBE_CHR_ORIENTATION", 
        "PROBE_COORDINATES", "CHROMOSOME", "TRANSCRIPT", "ILMN_GENE", "REFSEQ_ID", 
        "UNIGENE_ID", "SYMBOL", "PROTEIN_PRODUCT"), QC = TRUE)
```

```
## Inputting the data ...
## Perform Quality Control assessment of the LumiBatch object ...
## Directly converting probe sequence to nuIDs ...
```

```r

# check it
eset_raw
```

```
## Summary of data information:
## 	 Data File Information:
## 		GSGX Version	1.9.0
## 		Report Date	29/10/2013 14:56:38
## 		Project	BRC_GAP_Expression_02
## 		Group Set	BRC_GAP_Expression
## 		Analysis	BRC_GAP_Expression_nonorm_nobkgd
## 		Normalization	none
## 
## Major Operation History:
##             submitted            finished
## 1 2014-01-18 12:06:13 2014-01-18 12:11:45
## 2 2014-01-18 12:06:13 2014-01-18 12:11:45
##                                                                                                                    command
## 1 lumiR("/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt", 
## 2                                                       lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, QC = TRUE, 
##   lumiVersion
## 1      2.14.1
## 2      2.14.1
## ...
##             submitted            finished
## 8 2014-01-18 12:11:45 2014-01-18 12:15:55
## 9 2014-01-18 12:15:55 2014-01-18 12:16:02
##                                                                       command
## 8        lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)
## 9 addNuID2lumi(x.lumi = x.lumi, lib.mapping = lib.mapping, verbose = verbose)
##   lumiVersion
## 8      2.14.1
## 9      2.14.1
## 
## Object Information:
## LumiBatch (storageMode: lockedEnvironment)
## assayData: 47231 features, 618 samples 
##   element names: beadNum, detection, exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: 9020374058_A 9020374058_B ... 9249907052_L (618
##     total)
##   varLabels: sampleID
##   varMetadata: labelDescription
## featureData
##   featureNames: Ku8QhfS0n_hIOABXuE fqPEquJRRlSVSfL.8A ...
##     N8t5EuJCr0Tk9.zHno (47231 total)
##   fvarLabels: ProbeID TargetID ... DEFINITION (17 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation: lumiHumanAll.db 
## Control Data: Available
## QC information: Please run summary(x, 'QC') for details!
```

```r

# n_expression_chips
n_expression_chips <- dim(eset_raw)[2]
cat("  WARNING!: The number of expression chips=[", dim(eset_raw)[2], "]", "\r", 
    "\n")
```

```
##   WARNING!: The number of expression chips=[ 618 ] 
```


## 2. read in sample and batch information

```r

# gs_sample
if (is.na(gs_sample)) stop(" WARNING!: YOU HAVENT PROVIDED ANY SAMPLE INFORMATION!!!")
gs_sample
```

```
## [1] "/media/D/expression/GAP_Expression/final_reports_genomestudio/sample_table_Final_Report.txt"
```

```r

gs_sample_data <- read.table(paste(gs_sample), skip = 8, as.is = T, fill = T, 
    head = T, sep = "\t")
rownames(gs_sample_data) <- gs_sample_data$Sample.ID
gs_sample_data <- gs_sample_data[, names(gs_sample_data) != "X"]  # added this as genomestudio likes to add mystery columns to the end of this report 
n_samples <- dim(gs_sample_data)[1]  # number of rows ie samples

save(gs_sample_data, file = paste(out_dir, "/", project_name, ".gs_sample_data.RData", 
    sep = ""))

# pheno_file
if (is.na(pheno_file)) stop(" WARNING!: YOU HAVENT PROVIDED ANY PHENOTYPE INFORMATION!!!")
pheno_file
```

```
## [1] "/media/D/expression/GAP_Expression/final_reports_genomestudio/pheno_info.txt"
```

```r

pheno_dat <- read.table(paste(pheno_file), as.is = T, fill = T, head = T, sep = "\t")
save(pheno_dat, file = paste(out_dir, "/", project_name, ".pheno_dat.RData", 
    sep = ""))
has_pheno_cols <- c("Sample.ID", "SEX", "GROUPS", "TISSUE", "PHENOTYPE", "Study_ID") %in% 
    names(pheno_dat)
missing_pheno_cols <- "FALSE" %in% has_pheno_cols
if (missing_pheno_cols == "TRUE") stop(" WARNING!: YOU ARE MISSING ESSENTIAL SAMPLE INFORMATION! MAKE SURE YOUR PHENO_FILE HAS:- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID !!!")
raw_n_pheno_dat <- dim(pheno_dat)[1]  # number of rows ie samples

cat(" Running toupper() on PHENOTYPE, GROUP AND TISSUE variables to fix potential case issues", 
    "\r")
```

```
##  Running toupper() on PHENOTYPE, GROUP AND TISSUE variables to fix potential case issues 
```

```r
# fix case
pheno_dat$PHENOTYPE <- toupper(pheno_dat$PHENOTYPE)
pheno_dat$SEX <- toupper(pheno_dat$SEX)
pheno_dat$GROUPS <- toupper(pheno_dat$GROUPS)
pheno_dat$TISSUE <- toupper(pheno_dat$TISSUE)
# a quick looksee at counts
table(pheno_dat$PHENOTYPE)
```

```
## 
##    CASE CONTROL UNKNOWN 
##     393     190      25
```

```r
table(pheno_dat$GROUPS)
```

```
## 
##    CASE CONTROL UNKNOWN 
##     393     190      25
```

```r
table(pheno_dat$TISSUE)
```

```
## 
## BLOOD 
##   608
```

```r
table(pheno_dat$SEX)
```

```
## 
##  FEMALE    MALE UNKNOWN 
##     196     319      93
```

```r

# tech_pheno_file
if (is.na(tech_pheno_file)) stop(" WARNING!: YOU HAVENT PROVIDED ANY BATCH INFORMATION!!!")
tech_pheno_file
```

```
## [1] "/media/D/expression/GAP_Expression/final_reports_genomestudio/batch_info.txt"
```

```r

tech_pheno <- read.table(paste(tech_pheno_file), head = T, sep = "\t")
tech_pheno$Sentrix.Barcode <- as.character(tech_pheno$Sentrix.Barcode)
rownames(tech_pheno) <- tech_pheno$Sample.ID
colnames(tech_pheno) <- paste("tech.", names(tech_pheno), sep = "")
colnames(tech_pheno) <- c("Sample.ID", names(tech_pheno[, -1]))
save(tech_pheno, file = paste(out_dir, "/", project_name, ".tech_pheno.RData", 
    sep = ""))

# get pData()
eset_samples <- pData(eset_raw)
# add chip order and flad for 'has expression data' to eset_samples (pData)
eset_samples$has_expression <- 1
eset_samples$chip_order <- 1:dim(eset_samples)[1]
save(eset_samples, file = paste(out_dir, "/", project_name, ".eset_samples.RData", 
    sep = ""))

# col names
names(eset_samples)
```

```
## [1] "sampleID"       "has_expression" "chip_order"
```

```r
names(gs_sample_data)
```

```
##  [1] "Index"                 "Sample.ID"            
##  [3] "Sample.Group"          "Sentrix.Barcode"      
##  [5] "Sample.Section"        "Detected.Genes..0.01."
##  [7] "Detected.Genes..0.05." "Signal.Average"       
##  [9] "Signal.P05"            "Signal.P25"           
## [11] "Signal.P50"            "Signal.P75"           
## [13] "Signal.P95"            "BIOTIN"               
## [15] "CY3_HYB"               "HOUSEKEEPING"         
## [17] "LABELING"              "LOW_STRINGENCY_HYB"   
## [19] "NEGATIVE..background." "Noise"
```

```r
names(pheno_dat)
```

```
## [1] "Sample.ID"       "SEX"             "GROUPS"          "TISSUE"         
## [5] "PHENOTYPE"       "Study_ID"        "Study_ID_Count"  "GROUPS_ORIGINAL"
```

```r
names(tech_pheno)
```

```
##  [1] "Sample.ID"                           
##  [2] "tech.Sentrix.Barcode"                
##  [3] "tech.SampleSection"                  
##  [4] "tech.Batch"                          
##  [5] "tech.Date_out"                       
##  [6] "tech.Date_extraction"                
##  [7] "tech.person"                         
##  [8] "tech.Conc_Nanodrop"                  
##  [9] "tech.Date_Dilutionand_Amplification" 
## [10] "tech.Date_cRNApurification"          
## [11] "tech.Date_Quantitation_by_RiboGreen" 
## [12] "tech.Eluted_Total_labelled_cRNA"     
## [13] "tech.labelled_cRNA_Yield"            
## [14] "tech.concentration_of_labelled_cRNA" 
## [15] "tech.Date_labelled_cRNA"             
## [16] "tech.Date_Hybridization_for_15_hours"
## [17] "tech.Date_Washing_and_scanning"
```

```r

# head
head(eset_samples)
```

```
##                  sampleID has_expression chip_order
## 9020374058_A 9020374058_A              1          1
## 9020374058_B 9020374058_B              1          2
## 9020374058_C 9020374058_C              1          3
## 9020374058_D 9020374058_D              1          4
## 9020374058_E 9020374058_E              1          5
## 9020374058_F 9020374058_F              1          6
```

```r
head(gs_sample_data)
```

```
##              Index    Sample.ID Sample.Group Sentrix.Barcode
## 9020374058_A     1 9020374058_A 9020374058_A        9.02e+09
## 9020374058_B     2 9020374058_B 9020374058_B        9.02e+09
## 9020374058_C     3 9020374058_C 9020374058_C        9.02e+09
## 9020374058_D     4 9020374058_D 9020374058_D        9.02e+09
## 9020374058_E     5 9020374058_E 9020374058_E        9.02e+09
## 9020374058_F     6 9020374058_F 9020374058_F        9.02e+09
##              Sample.Section Detected.Genes..0.01. Detected.Genes..0.05.
## 9020374058_A              A                  5416                 10351
## 9020374058_B              B                  7670                 11874
## 9020374058_C              C                  6186                 10444
## 9020374058_D              D                  7718                 12057
## 9020374058_E              E                  7877                 11742
## 9020374058_F              F                  8623                 12306
##              Signal.Average Signal.P05 Signal.P25 Signal.P50 Signal.P75
## 9020374058_A            148         79         83         88        101
## 9020374058_B            152         79         83         87        100
## 9020374058_C            119         75         78         82         90
## 9020374058_D            157         79         84         88        104
## 9020374058_E            140         77         81         85         96
## 9020374058_F            179         82         87         92        112
##              Signal.P95 BIOTIN CY3_HYB HOUSEKEEPING LABELING
## 9020374058_A        279   7315    4016         2714     84.7
## 9020374058_B        313   7178    3970         2844     83.6
## 9020374058_C        189   7129    3893         1769     80.1
## 9020374058_D        330   7181    3835         3144     82.7
## 9020374058_E        255   7195    3832         2275     81.2
## 9020374058_F        396   6936    3721         4006     91.0
##              LOW_STRINGENCY_HYB NEGATIVE..background. Noise
## 9020374058_A               3948                  85.5  16.8
## 9020374058_B               3880                  83.6   6.9
## 9020374058_C               3800                  79.5   6.9
## 9020374058_D               3737                  85.0   8.9
## 9020374058_E               3706                  81.8   6.1
## 9020374058_F               3704                  88.0   6.8
```

```r
head(pheno_dat)
```

```
##      Sample.ID     SEX  GROUPS TISSUE PHENOTYPE Study_ID Study_ID_Count
## 1 9020374058_A    MALE CONTROL  BLOOD   CONTROL  SGAP393              1
## 2 9020374058_B UNKNOWN CONTROL  BLOOD   CONTROL  SGAP377              1
## 3 9020374058_C  FEMALE    CASE  BLOOD      CASE  SGAP464              1
## 4 9020374058_D    MALE CONTROL  BLOOD   CONTROL  SGAP331              1
## 5 9020374058_E    MALE CONTROL  BLOOD   CONTROL   GAP625              1
## 6 9020374058_F    MALE    CASE  BLOOD      CASE  SGAP116              1
##   GROUPS_ORIGINAL
## 1        baseline
## 2        baseline
## 3        baseline
## 4        baseline
## 5        baseline
## 6        baseline
```

```r
head(tech_pheno)
```

```
##                 Sample.ID tech.Sentrix.Barcode tech.SampleSection
## 9020374058_A 9020374058_A           9020374058                  A
## 9020374058_B 9020374058_B           9020374058                  B
## 9020374058_C 9020374058_C           9020374058                  C
## 9020374058_D 9020374058_D           9020374058                  D
## 9020374058_E 9020374058_E           9020374058                  E
## 9020374058_F 9020374058_F           9020374058                  F
##              tech.Batch tech.Date_out tech.Date_extraction tech.person
## 9020374058_A          1    13/05/2013           14/05/2013           1
## 9020374058_B          1    13/05/2013           14/05/2013           1
## 9020374058_C          1    21/05/2013           22/05/2013           1
## 9020374058_D          1    21/05/2013           22/05/2013           1
## 9020374058_E          1    13/05/2013           14/05/2013           1
## 9020374058_F          1    13/05/2013           14/05/2013           1
##              tech.Conc_Nanodrop tech.Date_Dilutionand_Amplification
## 9020374058_A              62.89                          16/07/2013
## 9020374058_B              52.75                          16/07/2013
## 9020374058_C              98.97                          16/07/2013
## 9020374058_D              78.19                          16/07/2013
## 9020374058_E              56.07                          16/07/2013
## 9020374058_F              37.58                          16/07/2013
##              tech.Date_cRNApurification
## 9020374058_A                 17/07/2013
## 9020374058_B                 17/07/2013
## 9020374058_C                 17/07/2013
## 9020374058_D                 17/07/2013
## 9020374058_E                 17/07/2013
## 9020374058_F                 17/07/2013
##              tech.Date_Quantitation_by_RiboGreen
## 9020374058_A                          23/07/2013
## 9020374058_B                          23/07/2013
## 9020374058_C                          23/07/2013
## 9020374058_D                          23/07/2013
## 9020374058_E                          23/07/2013
## 9020374058_F                          23/07/2013
##              tech.Eluted_Total_labelled_cRNA tech.labelled_cRNA_Yield
## 9020374058_A                              40                     9884
## 9020374058_B                              40                    13085
## 9020374058_C                              40                    22013
## 9020374058_D                              40                    17625
## 9020374058_E                              40                    16468
## 9020374058_F                              40                    10738
##              tech.concentration_of_labelled_cRNA tech.Date_labelled_cRNA
## 9020374058_A                               247.1              25/07/2013
## 9020374058_B                               327.1              25/07/2013
## 9020374058_C                               550.3              25/07/2013
## 9020374058_D                               440.6              25/07/2013
## 9020374058_E                               411.7              25/07/2013
## 9020374058_F                               268.5              25/07/2013
##              tech.Date_Hybridization_for_15_hours
## 9020374058_A                           25/07/2013
## 9020374058_B                           25/07/2013
## 9020374058_C                           25/07/2013
## 9020374058_D                           25/07/2013
## 9020374058_E                           25/07/2013
## 9020374058_F                           25/07/2013
##              tech.Date_Washing_and_scanning
## 9020374058_A                     26/07/2013
## 9020374058_B                     26/07/2013
## 9020374058_C                     26/07/2013
## 9020374058_D                     26/07/2013
## 9020374058_E                     26/07/2013
## 9020374058_F                     26/07/2013
```

```r

# quick check these should all have the same number of rows or samples!
dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(gs_sample_data)
```

```
## [1] 618  20
```

```r
dim(pheno_dat)
```

```
## [1] 608   8
```

```r
dim(tech_pheno)
```

```
## [1] 641  17
```

```r

# Venn of Sample.ID
ex <- eset_samples$sampleID
pp <- pheno_dat$Sample.ID
tt <- tech_pheno$Sample.ID
venninput <- list(ArrayExpression = ex, Batch_Info = tt, Pheno_Info = pp)
venn(venninput)
```

![plot of chunk readSamplePhenoAndBatchInfo](figure/readSamplePhenoAndBatchInfo.pdf) 

```r
# dev.off()
```


## 2.1 check for duplicate Study_ID

```r
# check for duplicate Study_ID
tab_id <- table(pheno_dat$Study_ID)
tab_id_df <- as.data.frame(tab_id)
colnames(tab_id_df) <- c("Study_ID", "Freq")
dupe_samples <- subset(tab_id_df, tab_id_df$Freq >= 2)
cat("  WARNING!: You have duplicate Study_IDs. N=[", dim(dupe_samples)[1], "]", 
    "\r", "\n")
```

```
##   WARNING!: You have duplicate Study_IDs. N=[ 59 ] 
```

```r
dupe_samples
```

```
##     Study_ID Freq
## 7    CGAP130    2
## 32   EUGE237    2
## 63    GAP555    2
## 78   GAP584L    2
## 94    GAP619    2
## 98    GAP626    2
## 110   GAP645    2
## 117   GAP660    2
## 118   GAP661    2
## 119   GAP663    2
## 126  GAP674L    2
## 128   GAP680    2
## 132   GAP689    2
## 141   GAP712    2
## 144   GAP716    2
## 151   GAP729    2
## 163  GAP736L    2
## 165   GAP741    2
## 168  GAP743C    2
## 191  GAP796C    2
## 256   GAP863    2
## 267   GAP892    2
## 269  GAP893C    2
## 280   GAP902    2
## 283   GAP904    2
## 285  GAP906C    2
## 291  GAP909C    2
## 328   GAP959    2
## 330   GAP962    2
## 335   GAP970    2
## 344   GAP990    2
## 349  LGAP103    2
## 351  LGAP134    2
## 361  LGAP171    2
## 369  SGAP136    2
## 381  SGAP161    2
## 382  SGAP163    2
## 384  SGAP169    2
## 387  SGAP175    2
## 396  SGAP197    2
## 404  SGAP208    2
## 408  SGAP232    2
## 417  SGAP250    2
## 421  SGAP260    2
## 423  SGAP265    2
## 439  SGAP309    2
## 442  SGAP316    2
## 444  SGAP320    2
## 445  SGAP322    2
## 466  SGAP348    2
## 469  SGAP353    2
## 472  SGAP358    2
## 475  SGAP363    2
## 496  SGAP391    2
## 502  SGAP399    2
## 511  SGAP408    2
## 520  SGAP421    2
## 522  SGAP424    2
## 546  SGAP473    2
```

```r
write.table(dupe_samples, file = paste(out_dir, "/", project_name, ".dupe_Study_IDs.txt", 
    sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
# n_unique_study_id
n_unique_study_id <- length(unique(tab_id_df$Study_ID))
cat("  WARNING!: The number of unique Study_Ids=[", n_unique_study_id, "]", 
    "\r", "\n")
```

```
##   WARNING!: The number of unique Study_Ids=[ 549 ] 
```


## 3. check eset_samples, sample & batch info Sample.ID's match in names and numbers & merge all

```r
cat(" Megreing pdata, pheno data, batch information adn genomestudio samople data", 
    "\r", "\n")
```

```
##  Megreing pdata, pheno data, batch information adn genomestudio samople data 
```

```r
# 1. merge eset_samples with pheno_dat.  Keep ALL overlaps only
eset_pheno_merge <- merge(eset_samples, pheno_dat, by.x = "sampleID", by.y = "Sample.ID")
eset_pheno_merge <- eset_pheno_merge[order(eset_pheno_merge$chip_order), ]
dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(eset_pheno_merge)  # check size
```

```
## [1] 608  10
```

```r

# 2. merge eset_pheno_merge with tech_pheno
eset_pheno_batch_merge <- merge(eset_pheno_merge, tech_pheno, by.x = "sampleID", 
    by.y = "Sample.ID")
eset_pheno_batch_merge <- eset_pheno_batch_merge[order(eset_pheno_batch_merge$chip_order), 
    ]
dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(eset_pheno_merge)
```

```
## [1] 608  10
```

```r
dim(eset_pheno_batch_merge)  # check size
```

```
## [1] 608  26
```

```r

# 3. merge all with genomestudio final report
eset_pheno_batch_gs_merge <- merge(eset_pheno_batch_merge, gs_sample_data, by.x = "sampleID", 
    by.y = "Sample.ID")
eset_pheno_batch_gs_merge <- eset_pheno_batch_gs_merge[order(eset_pheno_batch_gs_merge$chip_order), 
    ]

# final look at numbers in each merged data set
dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(eset_pheno_merge)
```

```
## [1] 608  10
```

```r
dim(eset_pheno_batch_merge)
```

```
## [1] 608  26
```

```r
dim(eset_pheno_batch_gs_merge)
```

```
## [1] 608  45
```

```r

# names
names(eset_pheno_batch_gs_merge)
```

```
##  [1] "sampleID"                            
##  [2] "has_expression"                      
##  [3] "chip_order"                          
##  [4] "SEX"                                 
##  [5] "GROUPS"                              
##  [6] "TISSUE"                              
##  [7] "PHENOTYPE"                           
##  [8] "Study_ID"                            
##  [9] "Study_ID_Count"                      
## [10] "GROUPS_ORIGINAL"                     
## [11] "tech.Sentrix.Barcode"                
## [12] "tech.SampleSection"                  
## [13] "tech.Batch"                          
## [14] "tech.Date_out"                       
## [15] "tech.Date_extraction"                
## [16] "tech.person"                         
## [17] "tech.Conc_Nanodrop"                  
## [18] "tech.Date_Dilutionand_Amplification" 
## [19] "tech.Date_cRNApurification"          
## [20] "tech.Date_Quantitation_by_RiboGreen" 
## [21] "tech.Eluted_Total_labelled_cRNA"     
## [22] "tech.labelled_cRNA_Yield"            
## [23] "tech.concentration_of_labelled_cRNA" 
## [24] "tech.Date_labelled_cRNA"             
## [25] "tech.Date_Hybridization_for_15_hours"
## [26] "tech.Date_Washing_and_scanning"      
## [27] "Index"                               
## [28] "Sample.Group"                        
## [29] "Sentrix.Barcode"                     
## [30] "Sample.Section"                      
## [31] "Detected.Genes..0.01."               
## [32] "Detected.Genes..0.05."               
## [33] "Signal.Average"                      
## [34] "Signal.P05"                          
## [35] "Signal.P25"                          
## [36] "Signal.P50"                          
## [37] "Signal.P75"                          
## [38] "Signal.P95"                          
## [39] "BIOTIN"                              
## [40] "CY3_HYB"                             
## [41] "HOUSEKEEPING"                        
## [42] "LABELING"                            
## [43] "LOW_STRINGENCY_HYB"                  
## [44] "NEGATIVE..background."               
## [45] "Noise"
```


## 4. Subset raw ExpressionSet to matched/complete Sample.IDs & Update pData() slot.

```r

# samples in gene expression data
samples_eset <- pData(eset_raw)$sampleID
length(samples_eset)
```

```
## [1] 618
```

```r

# samples with complete data
samples_complete_data <- eset_pheno_batch_gs_merge$sampleID
length(samples_complete_data)
```

```
## [1] 608
```

```r

# samples to remove
samples_to_remove <- (samples_eset %in% samples_complete_data) == FALSE
samples_to_remove <- pData(eset_raw)$sampleID[samples_to_remove]
length(samples_to_remove)
```

```
## [1] 10
```

```r

# rename eset_raw & save
eset_raw_preqc <- eset_raw
save(eset_raw_preqc, file = paste(out_dir, "/", project_name, ".eset_raw_preqc.RData", 
    sep = ""))

# subset eset_raw
eset_raw <- removeSamples_eset_lumi(eset = eset_raw_preqc, sampleRemove = samples_to_remove)
```

```
## The sample names in the controlData don't match sampleNames(object).
```

```r
eset_raw
```

```
## Summary of data information:
## 	 Data File Information:
## 		GSGX Version	1.9.0
## 		Report Date	29/10/2013 14:56:38
## 		Project	BRC_GAP_Expression_02
## 		Group Set	BRC_GAP_Expression
## 		Analysis	BRC_GAP_Expression_nonorm_nobkgd
## 		Normalization	none
## 
## Major Operation History:
##             submitted            finished
## 1 2014-01-18 12:06:13 2014-01-18 12:11:45
## 2 2014-01-18 12:06:13 2014-01-18 12:11:45
##                                                                                                                    command
## 1 lumiR("/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt", 
## 2                                                       lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, QC = TRUE, 
##   lumiVersion
## 1      2.14.1
## 2      2.14.1
## ...
##              submitted            finished
## 9  2014-01-18 12:15:55 2014-01-18 12:16:02
## 10 2014-01-18 12:16:47 2014-01-18 12:16:50
##                                                                        command
## 9  addNuID2lumi(x.lumi = x.lumi, lib.mapping = lib.mapping, verbose = verbose)
## 10                                                     Subsetting 618 samples.
##    lumiVersion
## 9       2.14.1
## 10      2.14.1
## 
## Object Information:
## LumiBatch (storageMode: lockedEnvironment)
## assayData: 47231 features, 608 samples 
##   element names: beadNum, detection, exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: 9020374058_A 9020374058_B ... 9249907052_L (608
##     total)
##   varLabels: sampleID
##   varMetadata: labelDescription
## featureData
##   featureNames: Ku8QhfS0n_hIOABXuE fqPEquJRRlSVSfL.8A ...
##     N8t5EuJCr0Tk9.zHno (47231 total)
##   fvarLabels: ProbeID TargetID ... DEFINITION (17 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation: lumiHumanAll.db 
## Control Data: Available
## QC information: Please run summary(x, 'QC') for details!
```

```r

# update pData
old_pdata <- pData(eset_raw)
old_pdata$old_order <- 1:dim(old_pdata)[1]

# merge with eset_pheno_batch_gs_merge
new_pdata <- merge(old_pdata, eset_pheno_batch_gs_merge, by.x = "sampleID", 
    by.y = "sampleID", all = TRUE, sort = FALSE)
new_pdata <- new_pdata[order(new_pdata$old_order), ]

# remove columns old_order has_expression chip_order
new_pdata <- new_pdata[, -c(2, 3, 4)]

# update rownames
rownames(new_pdata) <- new_pdata$sampleID

# update pData slot
pData(eset_raw) <- new_pdata
dim(pData(eset_raw))
```

```
## [1] 608  43
```

```r

# check it
eset_raw
```

```
## Summary of data information:
## 	 Data File Information:
## 		GSGX Version	1.9.0
## 		Report Date	29/10/2013 14:56:38
## 		Project	BRC_GAP_Expression_02
## 		Group Set	BRC_GAP_Expression
## 		Analysis	BRC_GAP_Expression_nonorm_nobkgd
## 		Normalization	none
## 
## Major Operation History:
##             submitted            finished
## 1 2014-01-18 12:06:13 2014-01-18 12:11:45
## 2 2014-01-18 12:06:13 2014-01-18 12:11:45
##                                                                                                                    command
## 1 lumiR("/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt", 
## 2                                                       lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, QC = TRUE, 
##   lumiVersion
## 1      2.14.1
## 2      2.14.1
## ...
##              submitted            finished
## 9  2014-01-18 12:15:55 2014-01-18 12:16:02
## 10 2014-01-18 12:16:47 2014-01-18 12:16:50
##                                                                        command
## 9  addNuID2lumi(x.lumi = x.lumi, lib.mapping = lib.mapping, verbose = verbose)
## 10                                                     Subsetting 618 samples.
##    lumiVersion
## 9       2.14.1
## 10      2.14.1
## 
## Object Information:
## LumiBatch (storageMode: lockedEnvironment)
## assayData: 47231 features, 608 samples 
##   element names: beadNum, detection, exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: 9020374058_A 9020374058_B ... 9249907052_L (608
##     total)
##   varLabels: sampleID SEX ... Noise (43 total)
##   varMetadata: labelDescription
## featureData
##   featureNames: Ku8QhfS0n_hIOABXuE fqPEquJRRlSVSfL.8A ...
##     N8t5EuJCr0Tk9.zHno (47231 total)
##   fvarLabels: ProbeID TargetID ... DEFINITION (17 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation: lumiHumanAll.db 
## Control Data: Available
## QC information: Please run summary(x, 'QC') for details!
```

```r

## n_expression_chips_with_data
n_expression_chips_with_data <- dim(eset_raw)[2]
```


## 5. Add nuID to fData

```r
# Add nuID to fData
fData(eset_raw)$nuID <- rownames(fData(eset_raw))
head(fData(eset_raw))
```

```
##                    ProbeID TargetID  TRANSCRIPT ILMN_GENE   REFSEQ_ID
## Ku8QhfS0n_hIOABXuE 6450255      7A5 ILMN_183371       7A5 NM_182762.2
## fqPEquJRRlSVSfL.8A 2570615     A1BG ILMN_175569      A1BG NM_130786.2
## ckiehnugOno9d7vf1Q 6370619     A1BG  ILMN_18893      A1BG NM_130786.2
## x57Vw5B5Fbt5JUnQkI 2600039     A1CF  ILMN_18532      A1CF NM_138932.1
## ritxUH.kuHlYqjozpE 2650615     A1CF   ILMN_7300      A1CF NM_014576.2
## QpE5UiUgmJOJEkPXpc 5340672     A1CF ILMN_165661      A1CF NM_138933.1
##                    UNIGENE_ID ENTREZ_GENE_ID   ACCESSION SYMBOL
## Ku8QhfS0n_hIOABXuE                    346389 NM_182762.2    7A5
## fqPEquJRRlSVSfL.8A                         1 NM_130786.2   A1BG
## ckiehnugOno9d7vf1Q                         1 NM_130786.2   A1BG
## x57Vw5B5Fbt5JUnQkI                     29974 NM_138932.1   A1CF
## ritxUH.kuHlYqjozpE                     29974 NM_014576.2   A1CF
## QpE5UiUgmJOJEkPXpc                     29974 NM_138933.1   A1CF
##                    PROTEIN_PRODUCT     PROBE_ID PROBE_TYPE PROBE_START
## Ku8QhfS0n_hIOABXuE     NP_877439.2 ILMN_1762337          S        2725
## fqPEquJRRlSVSfL.8A     NP_570602.2 ILMN_2055271          S        3151
## ckiehnugOno9d7vf1Q     NP_570602.2 ILMN_1736007          S        2512
## x57Vw5B5Fbt5JUnQkI     NP_620310.1 ILMN_2383229          A        1826
## ritxUH.kuHlYqjozpE     NP_055391.2 ILMN_1806310          A        1893
## QpE5UiUgmJOJEkPXpc     NP_620311.1 ILMN_1779670          I         278
##                    CHROMOSOME PROBE_CHR_ORIENTATION PROBE_COORDINATES
## Ku8QhfS0n_hIOABXuE          7                     - 20147187-20147236
## fqPEquJRRlSVSfL.8A         19                     - 63548541-63548590
## ckiehnugOno9d7vf1Q         19                     - 63549180-63549229
## x57Vw5B5Fbt5JUnQkI         10                     - 52566586-52566635
## ritxUH.kuHlYqjozpE         10                     - 52566495-52566544
## QpE5UiUgmJOJEkPXpc         10                     - 52610479-52610528
##                                                                                         DEFINITION
## Ku8QhfS0n_hIOABXuE                          Homo sapiens putative binding protein 7a5 (7A5), mRNA.
## fqPEquJRRlSVSfL.8A                               Homo sapiens alpha-1-B glycoprotein (A1BG), mRNA.
## ckiehnugOno9d7vf1Q                               Homo sapiens alpha-1-B glycoprotein (A1BG), mRNA.
## x57Vw5B5Fbt5JUnQkI Homo sapiens APOBEC1 complementation factor (A1CF), transcript variant 2, mRNA.
## ritxUH.kuHlYqjozpE Homo sapiens APOBEC1 complementation factor (A1CF), transcript variant 1, mRNA.
## QpE5UiUgmJOJEkPXpc Homo sapiens APOBEC1 complementation factor (A1CF), transcript variant 3, mRNA.
##                                  nuID
## Ku8QhfS0n_hIOABXuE Ku8QhfS0n_hIOABXuE
## fqPEquJRRlSVSfL.8A fqPEquJRRlSVSfL.8A
## ckiehnugOno9d7vf1Q ckiehnugOno9d7vf1Q
## x57Vw5B5Fbt5JUnQkI x57Vw5B5Fbt5JUnQkI
## ritxUH.kuHlYqjozpE ritxUH.kuHlYqjozpE
## QpE5UiUgmJOJEkPXpc QpE5UiUgmJOJEkPXpc
```


## 6. Save updated raw ExpressionSet eset_raw

```r
# Save updated raw ExpressionSet eset_raw
save(eset_raw, file = paste(out_dir, "/", project_name, ".eset_raw.RData", sep = ""))
```


## 7. Write data files to out_dir for eset_raw

```r
# Write data files to out_dir for eset_raw
write_expression_files(eset = eset_raw, outfile = paste(out_dir, "/", project_name, 
    ".eset_raw", sep = ""))
```

```
##  Writing probe exprs matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.exprs_matrix.txt ]  
##  Writing probe se.exprs matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.se.exprs_matrix.txt ]  
##  Writing probe detection matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.detection_matrix.txt ]  
##  Writing probe beadNum matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.beadNum_matrix.txt ]  
##  Writing probe PCA matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.pca_matrix.txt ]  
##  Writing pData slot of eset and adding PCA data to [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.pData.txt ]  
##  Writing fData slot of eset [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.fData.txt ] 
```


## 8. Basic QC plots on eset_raw

```r
gx_qc_plots_lumi_2(eset_raw)
```

```
##  starting qc plots  
##  setting up data for qc plots  
##  get pheno data  
##  basic colours  
##  batch pheno data  
##  expression matrix and IAC  
##  fundamentalNetworkConcepts  
##  flashClust  
##  beging plotting boxplot 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw1.pdf) 

```
##  beging plotting density 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw2.pdf) 

```
##  beging plotting cv 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw3.pdf) 

```
##  beging plotting outlier 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw4.pdf) 

```
##  beging plotting sampleTree 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw5.pdf) 

```
##  begin PCA plots 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw6.pdf) ![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw7.pdf) ![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw8.pdf) ![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw9.pdf) ![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw10.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Sentrix.Barcode 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw11.pdf) 

```
##  begin looping through batch variable PCA plots  tech.SampleSection 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw12.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_out 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw13.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_extraction 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw14.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_Dilutionand_Amplification 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw15.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_cRNApurification 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw16.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_Quantitation_by_RiboGreen 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw17.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_labelled_cRNA 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw18.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_Hybridization_for_15_hours 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw19.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Date_Washing_and_scanning 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw20.pdf) 

```
##  begin looping through batch variable PCA plots  tech.Conc_Nanodrop 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw21.pdf) 

```
##  begin looping through batch variable PCA plots  tech.labelled_cRNA_Yield 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw22.pdf) 

```
##  begin looping through batch variable PCA plots  tech.concentration_of_labelled_cRNA 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw23.pdf) 

```
##  begin plotDendroAndColors and  heatmap.2 plots 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw24.pdf) ![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw25.pdf) ![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw26.pdf) 

```
##  begin SampleNetwork plots 
```

![plot of chunk basicQCplots_eset_raw](figure/basicQCplots_eset_raw27.pdf) 


## 9. SampleNetwork on eset_raw for all samples as a first pass

```r
datExprs <- exprs(eset_raw)
samle_names <- sampleNames(eset_raw)
IAC = cor(datExprs, method = "p", use = "p")
diag(IAC) = 0
A.IAC = ((1 + IAC)/2)^2  ## ADJACENCY MATRIX
# fundamentalNetworkConcepts
FNC = fundamentalNetworkConcepts(A.IAC)  ## WGCNA
K2 = FNC$ScaledConnectivity
Z.K = round((K2 - mean(K2))/sd(K2), 3)
Z.C = round((FNC$ClusterCoef - mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef), 3)
# cor K,C
rho <- signif(cor.test(Z.K, Z.C, method = "s")$estimate, 2)
rho_pvalue <- signif(cor.test(Z.K, Z.C, method = "s")$p.value, 2)
# Z.K_outliers
Z.K_outliers <- Z.K < -iac_sd_thrs
Z.K_outliers <- names(Z.K_outliers[Z.K_outliers == TRUE])
n_outliers <- length(Z.K_outliers)
mean_IAC <- signif(mean(IAC[upper.tri(IAC)]), 2)
min_Z.K <- min(Z.K)
cat(" Number of Z.K outliers=[", n_outliers, "]", "\r", "\n")
```

```
##  Number of Z.K outliers=[ 22 ] 
```

```r
cat(" mean_IAC=[", mean_IAC, "]", "\r", "\n")
```

```
##  mean_IAC=[ 0.93 ] 
```

```r
cat(" cor(Z.k,Z.C)=[", rho, "] P=[", rho_pvalue, "]", "\r", "\n")
```

```
##  cor(Z.k,Z.C)=[ 0.29 ] P=[ 7.2e-13 ] 
```

```r
# print chip ids
Z.K_outliers
```

```
##  [1] "9020374071_I" "9020374072_F" "9031356054_A" "9031356100_H"
##  [5] "9031356100_K" "9216457012_A" "9216457023_A" "9216457023_I"
##  [9] "9216457029_K" "9216457032_J" "9216457033_L" "9234921070_L"
## [13] "9234921082_J" "9234921083_J" "9234921100_A" "9234921100_E"
## [17] "9235792061_J" "9235792095_J" "9249896091_B" "9249896091_C"
## [21] "9249907011_A" "9249907031_F"
```

```r
# get these bad samples from pData and update pdata to include these metrics
pData(eset_raw)$Z.K_eset_raw <- Z.K
pData(eset_raw)$Z.C_eset_raw <- Z.C
pData(eset_raw)$cor_Z.K.Z.C_eset_raw <- rho
pData(eset_raw)$cor_p_Z.K.Z.C_eset_ra <- rho_pvalue
# take a look at these outliers
samples_out <- pData(eset_raw[, Z.K_outliers])
head(samples_out)
```

```
##                  sampleID     SEX  GROUPS TISSUE PHENOTYPE Study_ID
## 9020374071_I 9020374071_I UNKNOWN CONTROL  BLOOD   CONTROL   GAP932
## 9020374072_F 9020374072_F    MALE    CASE  BLOOD      CASE  SGAP415
## 9031356054_A 9031356054_A  FEMALE CONTROL  BLOOD   CONTROL   GAP578
## 9031356100_H 9031356100_H    MALE CONTROL  BLOOD   CONTROL  SGAP406
## 9031356100_K 9031356100_K  FEMALE UNKNOWN  BLOOD   UNKNOWN   GAP949
## 9216457012_A 9216457012_A    MALE    CASE  BLOOD      CASE  GAP651L
##              Study_ID_Count GROUPS_ORIGINAL tech.Sentrix.Barcode
## 9020374071_I              1        baseline           9020374071
## 9020374072_F              1        baseline           9020374072
## 9031356054_A              1        baseline           9031356054
## 9031356100_H              1        baseline           9031356100
## 9031356100_K              1 MotherofProband           9031356100
## 9216457012_A              1      12FollowUp           9216457012
##              tech.SampleSection tech.Batch tech.Date_out
## 9020374071_I                  I          1    22/05/2013
## 9020374072_F                  F          1    21/05/2013
## 9031356054_A                  A          1    25/06/2013
## 9031356100_H                  H          1    21/05/2013
## 9031356100_K                  K          1    22/05/2013
## 9216457012_A                  A          1    28/05/2013
##              tech.Date_extraction tech.person tech.Conc_Nanodrop
## 9020374071_I           22/05/2013           2              66.52
## 9020374072_F           22/05/2013           1             204.10
## 9031356054_A           26/06/2013           3              23.00
## 9031356100_H           22/05/2013           1             124.83
## 9031356100_K           23/05/2013           2             439.09
## 9216457012_A           29/05/2013           2              43.29
##              tech.Date_Dilutionand_Amplification
## 9020374071_I                          16/07/2013
## 9020374072_F                          16/07/2013
## 9031356054_A                          24/07/2013
## 9031356100_H                          16/07/2013
## 9031356100_K                          16/07/2013
## 9216457012_A                          19/07/2013
##              tech.Date_cRNApurification
## 9020374071_I                 17/07/2013
## 9020374072_F                 17/07/2013
## 9031356054_A                 25/07/2013
## 9031356100_H                 17/07/2013
## 9031356100_K                 17/07/2013
## 9216457012_A                 20/07/2013
##              tech.Date_Quantitation_by_RiboGreen
## 9020374071_I                          23/07/2013
## 9020374072_F                          23/07/2013
## 9031356054_A                          30/07/2013
## 9031356100_H                          23/07/2013
## 9031356100_K                          23/07/2013
## 9216457012_A                          23/07/2013
##              tech.Eluted_Total_labelled_cRNA tech.labelled_cRNA_Yield
## 9020374071_I                              40                    23060
## 9020374072_F                              40                    27554
## 9031356054_A                              40                    16584
## 9031356100_H                              40                    27428
## 9031356100_K                              40                    48132
## 9216457012_A                              40                     6107
##              tech.concentration_of_labelled_cRNA tech.Date_labelled_cRNA
## 9020374071_I                               576.5              25/07/2013
## 9020374072_F                               688.8              25/07/2013
## 9031356054_A                               414.6              05/08/2013
## 9031356100_H                               685.7              25/07/2013
## 9031356100_K                              1203.3              25/07/2013
## 9216457012_A                               152.7              24/07/2013
##              tech.Date_Hybridization_for_15_hours
## 9020374071_I                           25/07/2013
## 9020374072_F                           25/07/2013
## 9031356054_A                           05/08/2013
## 9031356100_H                           25/07/2013
## 9031356100_K                           25/07/2013
## 9216457012_A                           24/07/2013
##              tech.Date_Washing_and_scanning Index Sample.Group
## 9020374071_I                     26/07/2013    32 9020374071_I
## 9020374072_F                     26/07/2013    41 9020374072_F
## 9031356054_A                     06/08/2013    69 9031356054_A
## 9031356100_H                     26/07/2013   142 9031356100_H
## 9031356100_K                     26/07/2013   145 9031356100_K
## 9216457012_A                     25/07/2013   171 9216457012_A
##              Sentrix.Barcode Sample.Section Detected.Genes..0.01.
## 9020374071_I       9.020e+09              I                  5742
## 9020374072_F       9.020e+09              F                  4814
## 9031356054_A       9.031e+09              A                  5567
## 9031356100_H       9.031e+09              H                  6593
## 9031356100_K       9.031e+09              K                  5959
## 9216457012_A       9.216e+09              A                 10208
##              Detected.Genes..0.05. Signal.Average Signal.P05 Signal.P25
## 9020374071_I                 10288            115         73         77
## 9020374072_F                  9224            102         69         72
## 9031356054_A                  9545            118         73         76
## 9031356100_H                  9810            108         71         74
## 9031356100_K                  9285            109         72         75
## 9216457012_A                 14216            306        103        111
##              Signal.P50 Signal.P75 Signal.P95 BIOTIN CY3_HYB HOUSEKEEPING
## 9020374071_I         80         86        185   7029    3702         1730
## 9020374072_F         75         80        143   6928    3438         1429
## 9031356054_A         79         84        165   7167    3774         2506
## 9031356100_H         77         81        148   6512    3463         1429
## 9031356100_K         78         83        141   6797    3441         1596
## 9216457012_A        120        173        842   7138    3664         6314
##              LABELING LOW_STRINGENCY_HYB NEGATIVE..background. Noise
## 9020374071_I     76.6               3734                  77.7   6.6
## 9020374072_F     70.8               4048                  73.4   5.2
## 9031356054_A     78.7               3623                  77.3   7.1
## 9031356100_H     69.3               3433                  74.8   3.2
## 9031356100_K     72.9               3402                  75.9   4.1
## 9216457012_A    107.5               3845                 111.0   9.2
##              Z.K_eset_raw Z.C_eset_raw cor_Z.K.Z.C_eset_raw
## 9020374071_I       -2.533       -0.203                 0.29
## 9020374072_F       -3.314       -0.391                 0.29
## 9031356054_A       -5.737       -2.536                 0.29
## 9031356100_H       -3.045        0.110                 0.29
## 9031356100_K       -4.800       -0.513                 0.29
## 9216457012_A       -2.197       -1.348                 0.29
##              cor_p_Z.K.Z.C_eset_ra
## 9020374071_I               7.2e-13
## 9020374072_F               7.2e-13
## 9031356054_A               7.2e-13
## 9031356100_H               7.2e-13
## 9031356100_K               7.2e-13
## 9216457012_A               7.2e-13
```

```r
table(samples_out$GROUPS)
```

```
## 
##    CASE CONTROL UNKNOWN 
##      14       7       1
```

```r
table(samples_out$PHENOTYPE)
```

```
## 
##    CASE CONTROL UNKNOWN 
##      14       7       1
```

```r
table(samples_out$SEX)
```

```
## 
##  FEMALE    MALE UNKNOWN 
##       7      12       3
```

```r
table(samples_out$tech.Sentrix.Barcode)
```

```
## 
## 9020374071 9020374072 9031356054 9031356100 9216457012 9216457023 
##          1          1          1          2          1          2 
## 9216457029 9216457032 9216457033 9234921070 9234921082 9234921083 
##          1          1          1          1          1          1 
## 9234921100 9235792061 9235792095 9249896091 9249907011 9249907031 
##          2          1          1          2          1          1
```

```r
# Save updated raw ExpressionSet eset_raw
eset_raw_Z.K_outliers <- Z.K_outliers
save(eset_raw_Z.K_outliers, file = paste(out_dir, "/", project_name, ".eset_raw_Z.K_outliers.RData", 
    sep = ""))
save(eset_raw, file = paste(out_dir, "/", project_name, ".eset_raw.RData", sep = ""))
# Write data files to out_dir for eset_raw
write_expression_files(eset = eset_raw, outfile = paste(out_dir, "/", project_name, 
    ".eset_raw", sep = ""))
```

```
##  Writing probe exprs matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.exprs_matrix.txt ]  
##  Writing probe se.exprs matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.se.exprs_matrix.txt ]  
##  Writing probe detection matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.detection_matrix.txt ]  
##  Writing probe beadNum matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.beadNum_matrix.txt ]  
##  Writing probe PCA matrix [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.pca_matrix.txt ]  
##  Writing pData slot of eset and adding PCA data to [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.pData.txt ]  
##  Writing fData slot of eset [ /media/D/expression/GAP_Expression/GAP_lumi_processing_18_01_2014/GAP.eset_raw.fData.txt ] 
```

