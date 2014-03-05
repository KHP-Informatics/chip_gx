Illumina Gene Expression Analysis Workflow
===========================================
Illumina Chip based Gene Expression Pipeline (QC, Normalization, Differential Expression)  

This is an example of a typical workflow we apply at the NIHR BRC-MH Bioinformatics unit.  

***Please note that this is under constant development and may need tweaking for your needs***

### Basic routines included in workflow 

- Genomestudio Final Reports read into ***LumiBatch***  using `lumi`
- ***QC Plots***: Boxplots, PCA Plots, hierarchical clustering, Coloured Dendograms...
- ***Probe Detection*** : based on expression level greater than the mean on the NEGATIVE beads
- ***Gender Checks*** : XIST Probe expression
- ***SampleNetwork Outliers*** : (based on : BMC Syst Biol. 2012 Jun 12;6:63. doi: 10.1186/1752-0509-6-63.
Network methods for describing sample relationships in genomic datasets: application to Huntington's disease.
Oldham MC1, Langfelder P, Horvath S. URL http://ccforum.com/1752-0509/6/63)
- Analysis of ***Batch Effects***
- Removal of ***Batch Effects*** using `sva`, `ComBat` or `lm` PCA Batch regressions

More soon...

******

Making Lumi input files in Genomestudio
=========================================

## Creating a new gene expression project in GenomeStudio

- Select “File” > “New Project” > “Gene Expression” and follow the GenomeStudio Project Wizard.
- Select “Direct Hyb”.
- Locate the “GenomeStudio_project” folder created in step one under “Project Repository”.
- Create a “project name” as [PROJECT_NAME]_expression_[DATE]_01.
- In the next window, under “Repository” specify the location of the “Data” folder created in step one. Once selected this should show all available chips within the “Sentrix Array Products” sub-window.
- Select and highlight the chip barcodes which require processing. Use the shift key to select multiple barcodes.
- Use the “Add Selected Samples to Project” icon to add the chip barcodes into the “Project Data” window. Click “Next”.
- Under “Groupset” repeat the project name followed by “raw” e.g project_expression_131101_01_raw.
- In the “Sentrix Array Products” sub-window select and highlight the chip barcodes which require processing.
- Use the “Add Selected Samples to Project” icon to add the chip barcodes into the “Project Data” window. Click “Next”.
- Under the “Name” tab select “Default” and click “Finish”.
- A message window will appear stating: “GenomeStudio detected that some samples have missing bead types. Would you like to impute missing data?” select “No”. Poor performing samples will be removed at a later step, after which the samples can be imputed.

***NOTE:*** The data will take a long time to load, check using Windows Task Manager (Performance panel) that at least one core is maxed out (otherwise it could mean GS has hung unexpectedly, restart the whole process if that's the case

*******

## Initial quality control

- The screen will present you next with the various tables of data. Select "Sample Table"
- Click on the Scatter Graph, then look at the following plots
- x=index, y=P05, labels=Sample Id
- x=index, y=Average Signal, labels=Sample Id
- x=index, y=Detected Genes, labels=Sample Id
- In each case note any samples which deviate substantially from the main chord of data points, some discretion is required here. Refer back to the gene expression log spread sheet to investigate deviated samples. This spread sheet is created by the geneticist conducting the gene expression assay, and can indicate samples of low amplification, gene expression failure or any technical problems associated with the assay.

*******

## Excluding samples and imputing

- Select “Analysis” > “Manage Group Sets...”
- Under “Groupset” create a project name in the format of: [PROJECT_NAME]_expression_[DATE]_02 e.g “project_expression_130111_02”
- From the “Sentrix Array Products” sub window select and highlight the chip barcodes which require processing.
- Use the “Create a group for each selected sample” icon to add the chip barcodes into the “Project Groups” window. This will show each individual sample from the chips selected. Using the “Remove selected groups and samples from the project” icon, remove all samples which have failed the initial quality control step. Select “Next”
Under the “Name” tab select “Default”. Click “Finish”.
- A message window will appear stating: “GenomeStudio detected that some samples have missing bead types. Would you like to impute missing data?” select “Yes”.

*******


