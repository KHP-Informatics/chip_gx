#########################################################################
# -- Author: Amos Folarin                                               #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################


if(!require("limma")){cat("ERROR: cannot load limma package!"); quit(save="no", status=1)};
if(!require("lumiHumanAll.db")){cat("ERROR: cannot load lumiHumanAll.db package"); quit(save="no", status=1)};
if(!require("annotate")){cat("ERROR: cannot load annotate package"); quit(save="no", status=1)};

# For expressionsets post processing: bgcorrect, transform, normalise
# Calculate diff expression for each one, based on the exp design provided
# eset_mbcbNPLOG2TQN
# eset_mbcnNPLOG2RSN
# eset_mbcbNPVSTQN
# eset_mbcnNPVSTRSN 


#------------------------------------------------------------------------
# DESC: Carry out Illumina diff expression using limma package, by default
# we adopt a factorial design and pairwise group contrasts are extracted. 
# If a design matrix is provided (.csv) then more complex experimental 
# designs are possible
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# USAGE:
# Designed to be sourced from the master script, to avoid re-reading large 
# variables that already exist
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# ARGS:
# ******** These are required to already exist in the R env. *********
#1) sampleTypes: a vector of group type strings
#2) designMatrix: design matrix csv file
#3) expression sets with of processed data data
#------------------------------------------------------------------------

# check exists in calling env 
if(! exists(sampleTypes) && !exists(designMatrix) )
{
    print("no parameters: sampleTypes and or designMatrix");
    quit(save="no", status=1);
}
# check expression sets also exists #TODO
#..


#------------------------------------------------------------------------
# utility function for generating pairwise contrast vector 
#------------------------------------------------------------------------
make.pairwise.contrasts <- function(groupnames)
{
    m <- groupnames;
    l = NULL;
    contrasts <- NULL;
    for(i in 1:(length(groupnames)-1))
    {
        m <- groupnames[(i+1):length(groupnames)];  #trim 
        l <- paste(groupnames[i], m, sep="-");
        contrasts <- c(contrasts,l);
    }

return(contrasts)
}



#------------------------------------------------------------------------
# Function to carry out differential expression using Limma package
#------------------------------------------------------------------------

diff.exp.calc <- function(eset)
{
    #subset probes based on detection
    dataMatrix <- exprs(eset);
    presentCount <- detectionCall(eset);
    selDataMatrix <- dataMatrix[presentCount > 0,];
    probeList <- rownames(selDataMatrix);
    
    # Validate the design matrix dimensions
    d.d <- dim(designMatrix);
    if( length(levels(factor(sampleTypes))) == d.d[2] && ncol(pData(eset)) == d.d[1] )
    {
        designMat.check <- TRUE;
    }


    ## construct the data matrix from groups for factorial design
    if(!exists(designMatrix))
    {
        f <- factor(sampleTypes, levels=unique(sampleTypes));
        design <- model.matrix(~0 f);
        colnames(design) <- levels(f);
        fit <- lmFit(selDataMatrix, design);
        contrasts.vec <- make.pairwise.contrasts(levels(f));
        contrast.matrix <- makeContrasts(contrast.vec, levels=design);
        fit2 <- contrasts.fit(fit, contrast.matrix);
        fit2 <- eBayes(fit);
        
        ## Add gene symbols to gene properties
        geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db');
        geneName <- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1]);
        fit2$genes <- data.frame(ID= probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE);
    }
    else ## use the design matrix provided which should already list the contrasts required)
    {
        design <- read.csv(designMatrix, stringsAsFactors=FALSE);
        #TODO
        #...
        #...

        
    }






    ## print the top 10 genes
    print(topTable(fit, coef='95:5-100:0', adjust='fdr', number=10));

    ## get significant gene list with FDR adjusted p.values less than 0.01
    p.adj <- p.adjust(fit$p.value[,2]);              
    sigGene.adj <- probeList[ p.adj < 0.01];
    ## without FDR adjustment
    sigGene <- probeList[ fit$p.value[,2] < 0.01];



    #return  the results:


}


#------------------------------------------------------------------------
# Generate the Differential Expression top-tables per eset object  ## TODO: Maybe move this block into the main calling script then  explicitly add the sampleTypes and designMatrix as actual parameters!!! 
#------------------------------------------------------------------------
eset_mbcbNPLOG2TQN.tt <- diff.exp.calc(eset_mbcbNPLOG2TQN);
eset_mbcnNPLOG2RSN.tt <- diff.exp.calc(eset_mbcnNPLOG2RSN);
eset_mbcbNPVSTQN.tt <- diff.exp.calc(eset_mbcnNPLOG2RSN);
eset_mbcnNPVSTRSN.tt  <- diff.exp.calc(eset_mbcnNPLOG2RSN);


#------------------------------------------------------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#------------------------------------------------------------------------





