Microarry pre-processing for Illumina BeadArray data
=====================================================

## Project Directory and Files

```r

# Data Directory
project_dir=/media/D/expression/GAP_Expression

# Data Files
control_probe_profile_Final_Report.txt
Group_Probe_Profile_Final_Report.txt
probe_annotation_Final_Report.txt
Sample_and_Control_Probe_Profile_FinalReport.txt # input for lumiR
sample_table_Final_Report.txt
pheno_info.txt
batch_info.txt

```


# Begin R
Start R and call the following


```r
rm(list = ls())
options(stringsAsFactors = FALSE)
```


## setwd
This is the full path to the directory where all the raw genomestudio output is stored


```r
setwd("/media/D/sjnewhouse/GENE_EXPRESSION")
```


## load libs

```r
# Load libraries
library(lumi)
```

```
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Warning: replacing previous import by 'graphics::image' when loading
## 'methylumi'
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
```

```
## Warning: replacing previous import by 'nleqslv::nleqslv' when loading
## 'lumi'
```

```r
library(annotate)
```

```
## Loading required package: AnnotationDbi
```

```r
library(lumiHumanAll.db)
```

```
## Loading required package: org.Hs.eg.db
## Loading required package: DBI
## 
## 
## lumiHumanAll.db is using or is likely to need access to special
##   nuID identifiers.  Users can learn about these identifiers from
##   vignette documentation provided with the lumi package.
```

```r
library(affy)
```

```
## 
## Attaching package: 'affy'
## 
## The following objects are masked from 'package:lumi':
## 
##     MAplot, plotDensity
```

```r
library(cluster)
library(impute)
library(WGCNA)
```

```
## Loading required package: dynamicTreeCut
## Loading required package: flashClust
## 
## Attaching package: 'flashClust'
## 
## The following object is masked from 'package:stats':
## 
##     hclust
## 
## Loading required package: Hmisc
## Loading required package: grid
## Loading required package: lattice
## Loading required package: survival
## Loading required package: splines
## Loading required package: Formula
## Hmisc library by Frank E Harrell Jr
## 
## Type library(help='Hmisc'), ?Overview, or ?Hmisc.Overview')
## to see overall documentation.
## 
## 
## Attaching package: 'Hmisc'
## 
## The following object is masked from 'package:survival':
## 
##     untangle.specials
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     contents
## 
## The following object is masked from 'package:lumi':
## 
##     combine
## 
## The following objects are masked from 'package:Biobase':
## 
##     combine, contents
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     combine
## 
## The following objects are masked from 'package:base':
## 
##     format.pval, round.POSIXt, trunc.POSIXt, units
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

```
## 
## Attaching package: 'WGCNA'
## 
## The following object is masked from 'package:stats':
## 
##     cor
```

```r
library(gplots)
```

```
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
library(limma)
```

```
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
library(vsn)
library(MBCB)
```

```
## Loading required package: tcltk
## Loading required package: tcltk2
## 
## Attaching package: 'tcltk2'
## 
## The following objects are masked from 'package:Hmisc':
## 
##     label, label<-
```

```r
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)
```

```
## Loading required package: MASS
## 
## Attaching package: 'MASS'
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     select
## 
## Loading required package: boot
## 
## Attaching package: 'boot'
## 
## The following object is masked from 'package:survival':
## 
##     aml
## 
## The following object is masked from 'package:lattice':
## 
##     melanoma
## 
## Loading required package: survey
## 
## Attaching package: 'survey'
## 
## The following object is masked from 'package:Hmisc':
## 
##     deff
## 
## The following object is masked from 'package:graphics':
## 
##     dotchart
## 
## Loading required package: mitools
## Loading required package: foreign
## This is the global version of package relaimpo.
## 
## If you are a non-US user, a version with the interesting additional metric pmvd is available
## 
## from Ulrike Groempings web site at prof.beuth-hochschule.de/groemping.
```


## source processing functions

```r
# path to gene expression processing scripts
path_to_scripts <- "/media/D/sjnewhouse/GENE_EXPRESSION"
source(paste(path_to_scripts, "/sjnewhouse_misc_R.R", sep = ""))
ls()
```

```
##  [1] "basic_sampleNetwork"        "basic_sampleNetworkIterate"
##  [3] "data_summary_plots"         "has_var_probe"             
##  [5] "has_var_probe2"             "max_probe"                 
##  [7] "mean_probe"                 "min_probe"                 
##  [9] "negBeadOutlierRepMean"      "path_to_scripts"           
## [11] "quantfun"                   "removeSamples_eset_lumi"   
## [13] "sd_probe"                   "shuffle_cols"              
## [15] "shuffle_rows"               "var_probe"                 
## [17] "write_expression_files"     "zero_var_probe"
```


## set project settings and I/O
User is asked to manually provide the options


```r
# project directory
project_dir <- "/media/D/expression/GAP_Expression"

# set working dir again
setwd(project_dir)

# project name
project_name <- "GAPtest01"

# output directory for lumi process and plots
out_dir <- paste(project_name, "_lumi_processing", sep = "")

# make project pre-processing directory
make_dir_command <- paste(" if [ ! -e ./", out_dir, " ]; then mkdir ./", out_dir, 
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
norm_method <- "quantile"  ## 'rsn' # quantile, rsn, or both

# write settings to file
project_settings <- data.frame(project_dir = project_dir, project_name = project_name, 
    out_dir = out_dir, gs_report = gs_report, gs_probe = gs_probe, gs_sample = gs_sample, 
    gs_control = gs_control, anno_table = anno_table, pheno_file = pheno_file, 
    tech_pheno_file = tech_pheno_file, probe_det = probe_det, sample_det = sample_det, 
    sex_check = sex_check, iac_check = iac_check, iac_sd_thrs = iac_sd_thrs, 
    mbcb_method = mbcb_method, transform_method = transform_method, norm_method = norm_method)

# write table to out_dir
write.table(project_settings, file = paste(out_dir, "/", project_name, ".project_settings", 
    sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

# check settings
t(project_settings)
```

```
##                  [,1]                                                                                                            
## project_dir      "/media/D/expression/GAP_Expression"                                                                            
## project_name     "GAPtest01"                                                                                                     
## out_dir          "GAPtest01_lumi_processing"                                                                                     
## gs_report        "/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt"
## gs_probe         "/media/D/expression/GAP_Expression/final_reports_genomestudio/Group_Probe_Profile_Final_Report.txt"            
## gs_sample        "/media/D/expression/GAP_Expression/final_reports_genomestudio/sample_table_Final_Report.txt"                   
## gs_control       "/media/D/expression/GAP_Expression/final_reports_genomestudio/control_probe_profile_Final_Report.txt"          
## anno_table       "/media/D/expression/GAP_Expression/final_reports_genomestudio/probe_annotation_Final_Report.txt"               
## pheno_file       "/media/D/expression/GAP_Expression/final_reports_genomestudio/pheno_info.txt"                                  
## tech_pheno_file  "/media/D/expression/GAP_Expression/final_reports_genomestudio/batch_info.txt"                                  
## probe_det        "80"                                                                                                            
## sample_det       "80"                                                                                                            
## sex_check        "1"                                                                                                             
## iac_check        "1"                                                                                                             
## iac_sd_thrs      "2"                                                                                                             
## mbcb_method      "MLE"                                                                                                           
## transform_method "log2"                                                                                                          
## norm_method      "quantile"
```


BEGIN PRE-PROCESSING
=====================

Raw Expression Set
-------------------

## 1. read raw gene expression data 

```r

# raw input
gs_report
```

```
## [1] "/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt"
```

```r

# read raw gene expression data from genomestudio reports and create
# ExpressionSet

eset_raw <- lumiR(paste(gs_report), lib.mapping = "lumiHumanIDMapping", annotationColumn = c("PROBE_ID", 
    "CHROMOSOME", "SYMBOL", "DEFINITION", "ACCESSION", "ENTREZ_GENE_ID", "PROBE_TYPE", 
    "PROBE_START", "PROBE_SEQUENCE", "PROBE_CHR_ORIENTATION", "PROBE_COORDINATES", 
    "CHROMOSOME", "TRANSCRIPT", "ILMN_GENE", "REFSEQ_ID", "UNIGENE_ID", "SYMBOL", 
    "PROTEIN_PRODUCT"))
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
## 1 2014-01-08 18:18:18 2014-01-08 18:23:44
## 2 2014-01-08 18:18:18 2014-01-08 18:23:44
##                                                                                                                                       command
## 1                    lumiR("/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt", 
## 2     annotationColumn = c("/media/D/expression/GAP_Expression/final_reports_genomestudio/Sample_and_Control_Probe_Profile_FinalReport.txt", 
##   lumiVersion
## 1      2.14.1
## 2      2.14.1
## ...
##             submitted            finished
## 7 2014-01-08 18:23:44 2014-01-08 18:27:51
## 8 2014-01-08 18:27:51 2014-01-08 18:27:57
##                                                                       command
## 7        lumiQ(x.lumi = x.lumi, detectionTh = detectionTh, verbose = verbose)
## 8 addNuID2lumi(x.lumi = x.lumi, lib.mapping = lib.mapping, verbose = verbose)
##   lumiVersion
## 7      2.14.1
## 8      2.14.1
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

raw_sample_data <- read.table(paste(gs_sample), skip = 8, as.is = T, fill = T, 
    head = T, sep = "\t")
rownames(raw_sample_data) <- raw_sample_data$Sample.ID
raw_sample_data <- raw_sample_data[, names(raw_sample_data) != "X"]  # added this as genomestudio likes to add mystery columns to the end of this report 
raw_n_samples <- dim(raw_sample_data)[1]  # number of rows ie samples
save(raw_sample_data, file = paste(out_dir, "/", project_name, ".raw_sample_data.RData", 
    sep = ""))
```

```
## Warning: cannot open compressed file
## 'GAPtest01_lumi_processing/GAPtest01.raw_sample_data.RData', probable
## reason 'No such file or directory'
```

```
## Error: cannot open the connection
```

```r

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
```

```
## Warning: cannot open compressed file
## 'GAPtest01_lumi_processing/GAPtest01.pheno_dat.RData', probable reason 'No
## such file or directory'
```

```
## Error: cannot open the connection
```

```r
has_pheno_cols <- c("Sample.ID", "SEX", "GROUPS", "TISSUE", "PHENOTYPE", "Study_ID") %in% 
    names(pheno_dat)
missing_pheno_cols <- "FALSE" %in% has_pheno_cols
if (missing_pheno_cols == "TRUE") stop(" WARNING!: YOU ARE MISSING ESSENTIAL SAMPLE INFORMATION! MAKE SURE YOUR PHENO_FILE HAS:- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID !!!")
raw_n_pheno_dat <- dim(pheno_dat)[1]  # number of rows ie samples

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
colnames(tech_pheno) <- paste("tech.", names(tech_pheno), sep = "")
save(tech_pheno, file = paste(out_dir, "/", project_name, ".tech_pheno.RData", 
    sep = ""))
```

```
## Warning: cannot open compressed file
## 'GAPtest01_lumi_processing/GAPtest01.tech_pheno.RData', probable reason
## 'No such file or directory'
```

```
## Error: cannot open the connection
```

```r

# get pData()
eset_samples <- pData(eset_raw)
# add chip order and flad for 'has expression data' to eset_samples (pData)
eset_samples$has_expression <- 1
eset_samples$chip_order <- 1:dim(eset_samples)[1]
save(eset_samples, file = paste(out_dir, "/", project_name, ".eset_samples.RData", 
    sep = ""))
```

```
## Warning: cannot open compressed file
## 'GAPtest01_lumi_processing/GAPtest01.eset_samples.RData', probable reason
## 'No such file or directory'
```

```
## Error: cannot open the connection
```

```r

# col names
names(eset_samples)
```

```
## [1] "sampleID"       "has_expression" "chip_order"
```

```r
names(raw_sample_data)
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
## [1] "Sample.ID" "GROUPS"    "SEX"       "TISSUE"    "PHENOTYPE" "Study_ID"
```

```r
names(tech_pheno)
```

```
##  [1] "tech.Sample.ID"                      
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

# quick check these should all have the same number of rows or samples!
dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(raw_sample_data)
```

```
## [1] 618  20
```

```r
dim(pheno_dat)
```

```
## [1] 618   6
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
tt <- tech_pheno$tech.Sample.ID
gs <- raw_sample_data$Sample.ID
venninput <- list(ArrayExpression = ex, Batch_Info = tt, Pheno_Info = pp, GenomeStudio = gs)
venn(venninput)
```

![plot of chunk readSamplePhenoAndBatchInfo](figure/readSamplePhenoAndBatchInfo.png) 


## 3. check eset_samples, sample & batch info Sample.ID's match in names and numbers & merge all

```r

# 1. merge eset_samples with pheno_dat. Keep ALL overlaps only
eset_pheno_merge <- merge(eset_samples, pheno_dat, by.x. = "sampleID", by.y = "Sample.ID")
```

```
## Error: 'by.x' and 'by.y' specify different numbers of columns
```

```r
eset_pheno_merge <- eset_pheno_merge[order(eset_pheno_merge$chip_order), ]
```

```
## Error: object 'eset_pheno_merge' not found
```

```r
dim(eset_samples)
```

```
## [1] 618   3
```

```r
eset_pheno_merge  # check size
```

```
## Error: object 'eset_pheno_merge' not found
```

```r

# 2. merge eset_pheno_merge with tech_pheno
eset_pheno_batch_merge <- merge(eset_pheno_merge, tech_pheno, by.x. = "Sample.ID", 
    by.y = "tech.Sample.ID")
```

```
## Error: object 'eset_pheno_merge' not found
```

```r
eset_pheno_batch_merge <- eset_pheno_batch_merge[order(eset_pheno_batch_merge$chip_order), 
    ]
```

```
## Error: object 'eset_pheno_batch_merge' not found
```

```r
dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(eset_pheno_merge)
```

```
## Error: object 'eset_pheno_merge' not found
```

```r
dim(eset_pheno_batch_merge)  # check size
```

```
## Error: object 'eset_pheno_batch_merge' not found
```

```r

# 3. merge all with genomestudio final report
eset_pheno_batch_gs_merge <- merge(eset_pheno_batch_merge, gs_sample, by.x. = "Sample.ID", 
    by.y = "Sample.ID")
```

```
## Error: object 'eset_pheno_batch_merge' not found
```

```r
eset_pheno_batch_gs_merge <- eset_pheno_batch_gs_merge[order(eset_pheno_batch_gs_merge$chip_order), 
    ]
```

```
## Error: object 'eset_pheno_batch_gs_merge' not found
```

```r

dim(eset_samples)
```

```
## [1] 618   3
```

```r
dim(eset_pheno_merge)
```

```
## Error: object 'eset_pheno_merge' not found
```

```r
dim(eset_pheno_batch_merge)
```

```
## Error: object 'eset_pheno_batch_merge' not found
```

```r
dim(eset_pheno_batch_gs_merge)  # check size
```

```
## Error: object 'eset_pheno_batch_gs_merge' not found
```

```r

```


## 4. Update pData() slot and subset raw ExpressionSet to matched/complete Sample.IDs



## 5. Save raw ExpressionSet eset_raw



## 6. Write data files to out_dir for eset_raw



## 7. Basic QC and plots on eset_raw





##



MBCB (Model-based Background Correction for Beadarray)
-------------------------------------------------------




