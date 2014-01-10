Microarry pre-processing for Illumina BeadArray data.
========================================================

- Authors : Dr Stephen J Newhouse  
- Email :   stephen.j.newhouse@gmail.com  
- Verison:  1.00  
- Source:   Genomestudio  

# REQUIRED INPUT FILES

1. Firnal reports from Genomestudio:  
    * Sample  
    * Group Probe   
    * Control Probe   
    * Probe annotations  
    
2. Phenotype file:  must contain the following columns:- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID  

3. Technical information file with batch processing and lab realted data.  This must contain the following columns:- Sample.ID followed by any and all processing data information, RIN, RNA yeilds and concetrations.

4. Config File

5. sjnewhouse_misc_R.R

# PROJECT DIRECTORY
This should conatin all required input files listed above

```
# move to project direcotory
cd /media/D/expression/GAP_Expression;
# start R
R
```

### Load Libs

```r
# Load Libs
rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("/media/D/expression/GAP_Expression")
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
##     anyDuplicated, as.data.frame, cbind, colnames, duplicated,
##     eval, Filter, Find, get, intersect, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rep.int, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist
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

```r
source("sjnewhouse_misc_R.R")
```

```
## Warning: cannot open file 'sjnewhouse_misc_R.R': No such file or directory
```

```
## Error: cannot open the connection
```


Email Dr Stephen J Newhouse @ stephen.newhouse@kcl.ac.uk for source code **sjnewhouse\_misc\_R.R**.
This contains a number of fucntion required for running the pipeline. It Will be on git soon.

# CONFIG FILE
User is to set this up in excel. Defines working directory, file names and pre-processing options.
Make sure 


```r
# Read Config File
config_file <- "GAP.project.config"
project_settings <- read.table(config_file, head = TRUE, sep = "\t")  ## as.is=T, fill=TRUE)
```

```
## Warning: cannot open file 'GAP.project.config': No such file or directory
```

```
## Error: cannot open the connection
```

```r
# check config file
t(project_settings)
```

```
## Error: object 'project_settings' not found
```


### Getting Project settings and file names

```r
# Getting Project settings and file names
project_dir <- project_settings$project_dir  ## /scratch/project/pipelines/DATA/expression/GAP_Expression
```

```
## Error: object 'project_settings' not found
```

```r
gs_report <- project_settings$gs_report  ## 'Sample_and_Control_Probe_Profile_FinalReport.txt'
```

```
## Error: object 'project_settings' not found
```

```r
gs_probe <- project_settings$gs_probe  ## 'Group_Probe_Profile_Final_Report.txt' # genomestudio report
```

```
## Error: object 'project_settings' not found
```

```r
gs_sample <- project_settings$gs_sample  ## 'sample_table_Final_Report.txt' # genomestudio report
```

```
## Error: object 'project_settings' not found
```

```r
gs_control <- project_settings$gs_control  ## 'control_probe_profile_Final_Report.txt' # genomestudio report
```

```
## Error: object 'project_settings' not found
```

```r
anno_table <- project_settings$anno_table  ## annotation table for project
```

```
## Error: object 'project_settings' not found
```

```r
pheno_file <- project_settings$pheno_file  ## FILE NAME must contain : Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID
```

```
## Error: object 'project_settings' not found
```

```r
project_name <- project_settings$project_name  ## 'GAP' # GAP
```

```
## Error: object 'project_settings' not found
```

```r
out_dir <- paste(project_name, "_lumi_processing", sep = "")  # output dir for lumi process and plots
```

```
## Error: object 'project_name' not found
```

```r
probe_det <- project_settings$probe_det  ## 80 # 50, 80, 90,100
```

```
## Error: object 'project_settings' not found
```

```r
sample_det <- project_settings$sample_det  ## 80 # 50,80,90,100
```

```
## Error: object 'project_settings' not found
```

```r
sex_check <- project_settings$sex_check  ## 1 # 1 or 0
```

```
## Error: object 'project_settings' not found
```

```r
iac_check <- project_settings$iac_check  ## 1 # 1 or 0
```

```
## Error: object 'project_settings' not found
```

```r
iac_sd_thrs <- project_settings$iac_sd_thrs  ## 2 #
```

```
## Error: object 'project_settings' not found
```

```r
norm_method <- project_settings$norm_method  ## 'rsn' # quantile, rsn, or both
```

```
## Error: object 'project_settings' not found
```

```r
transform_method <- project_settings$transform_method  ## 'vst' # log2, vst or both
```

```
## Error: object 'project_settings' not found
```

```r
mbcb_method <- project_settings$mbcb_method  ## NP or MLE
```

```
## Error: object 'project_settings' not found
```

```r
tech_pheno_file <- project_settings$tech_pheno_file  ## 'sample_tech_info.txt'
```

```
## Error: object 'project_settings' not found
```


### Set working directory

```r
# set working directory
setwd(project_dir)
```

```
## Error: object 'project_dir' not found
```

```r
cat(" Making processing directory:- ", paste(" mkdir ./", out_dir, sep = ""), 
    "\r", "\n")
```

```
## Error: object 'out_dir' not found
```

```r
make_dir_command <- paste(" if [ ! -e ./", out_dir, " ]; then mkdir ./", out_dir, 
    "; fi", sep = "")
```

```
## Error: object 'out_dir' not found
```

```r
system(make_dir_command)
```

```
## Error: object 'make_dir_command' not found
```



# THE END



