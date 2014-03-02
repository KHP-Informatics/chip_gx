Illumina Gene Expression Analysis Workflow
===========================================

Illumina Chip based Gene Expression Pipeline (QC, Normalization, Differential Expression)  

This is an example of a typical workflow we apply at the NIHR BRC-MH Bioinformatics unit.  

Please note that this is under constant development and tweaking.  

### Basic routines included in workflow 

- Genomestudio Final Reports read into ***LumiBatch***  using `lumiR()`
- QC Plots: Boxplots, PCA Plots, hierarchical clustering, Coloured Dendograms...
- Probe Detection based on expression level greater than the mean on the NEGATIVE beads
- Gender Checks : XIST Probe expression
- SampleNetwork Outliers : REF TO COME
- Analysis of ***Batch Effects***
- SVA, ComBat or `rlm` PCA Batch regressions to help remove potential batch effects  

