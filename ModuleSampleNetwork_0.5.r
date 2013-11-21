## Last updated: 10.27.11

## Recent changes:
## 1) Fixed cexlabels1 so it works.

## Recent changes:
## 1) Increased bottom margin size for some plots.
## 2) Fixed bug so infinite p-values dont break plots.

## This document contains R functions for constructing moduel sample networks from genomic (or other) datasets.
## The ModuleSampleNetwork R function is similar to the SampleNetwork R function, but does not provide interactivity and the abilities to identify and remove outliers or perform normalization.
## Instead, this function is designed to enable the comparison of sample network properties for subgroups of samples and subsets of features.
## See Figures 4 and 6 from the corresponding journal article for examples of ModuleSampleNetwork output.
## In our applications, subsets of feature are often defined as modules of coexpressed genes, but other definitions could be used.
## A tutorial describing the usage of these functions is available on our web site: http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/SampleNetwork.
## The ModuleSampleNetwork function was written by Mike Oldham (oldhamm@stemcell.ucsf.edu).
## DISCLAIMER: THIS CODE IS FREELY PROVIDED FOR ACADEMIC USE WITH ABSOLUTELY NO WARRANTIES.  PLEASE CONTACT MIKE OLDHAM WITH BUG REPORTS OR SUGGESTIONS.
## The ModulePrinComps1 function was written by Steve Horvath (shorvath@mednet.ucla.edu).

## To cite this code or the methods contained therein, please use the following references:

# 1) Oldham MC, Langfelder P, Horvath S (2011) "Sample Networks for Enhancing Cluster Analysis of Genomic Data: Application to Huntington's Disease".  Submitted.

## The function "ModuleSampleNetwork" takes as input three files:
## 1) a feature activity file (argument "datExprT"; rows = features (e.g. probe sets), columns = feature information (e.g. probe set IDs, gene symbols, etc.) + samples);
## 2) a feature category file (argument "modules1"; rows = features (e.g. probe sets), columns = feature category information (e.g. probe set IDs [or the unique feature identifier that appears in datExprT] and module assignments [e.g. colors, numbers, etc.]);
## 3) a sample information file (argument "sampleinfo1"; rows = samples, columns = sample traits);
## Additional arguments:
## "method1" = distance measure to be used fo constructing sample networks.  choices are "correlation" and "euclidean".
## "impute1" = logical (T/F): should missing values (must be coded as "NA") be imputed (default = FALSE);
## "skip1" = integer describing the number of feature info columns in the expression matrix (there must be at least one);
## "indices1" = a list of vectors to subset the columns in the expression matrix (each vector in this list defines a separate group of samples for processing, e.g. "CTRL", "EXP", etc.);
## "featurelabels1" = a vector equal in length to the number of features, containing convenient feature labels, which will appear in some plots;
## "subgroup1" = integer that points to the column in sampleinfo1 that specifies the sample subgroups to be compared; these subgroups will be colored separately in plots;
## "samplelabels1" = integer that points to the column in sampleinfo1 with the sample labels that will appear in plots (note: these must be identical to the sample column headers in datExprT);
## "grouplabels1" = integer that points to the column in the sample information file that containes the group labels (note: number of groups must equal number of indices);
## "fitmodels1" = logical (T/F): should model-fitting be performed?  if FALSE (default), standard clustering by 1-ISA with accompanying Z.K plot will be produced for outlier removal;
## "whichmodel1" = should a univariate (default) or multivariate linear regression model be used?
## "whichfit1" = variable that specifies the quantity (y-variable) that will be regressed during model-fitting; choices are "pc1" (default; the 1st principal component of the entire dataset [defined by subset1]), "mean" (sample mean), or "K" (sample connectivity);
## "btrait1" = a vector of integers that identifies the columns (traits) in sampleinfo1 that should be included in model-fitting (only specified if fitmodels1 = TRUE; otherwise, NULL [default]);
## "trait1" =  a vector of integers that identifies the columns (traits) in sampleinfo1 that should be tested for the significance of individual factor levels in model-fitting; trait1 must be a subset of btrait1 and specify only categorical variables; only specified if fitmodels1 = TRUE; otherwise, NULL [default];
## "asfactors1" = a vector of integers that identifies the columns (variables) in sampleinfo1 that have been assigned to btrait1 and should be treated as factors;
## "projectname1" = a character label for the project that will appear in some output files;
## "cexlabels1" = a scaling factor for the sample labels in plots (default = 0.7);
## "verbose" = logical (T/F): if TRUE (default), additional network metrics will be exported as text on screen, in plots, and in "..._metrics.csv" output file.  Metrics are produced using the fundamentalNetworkConcepts function from WGCNA library;
## "minGpSize1" = integer that specifies the minimum number of samples that must comprise a subgroup (default = 10);
## "corType1" = specifies the type of correlation to be used when calculating cor(K,C); choices are "s" (Spearman; default), "p" (Pearson), and "k" (Kendall);
## "logT" = logical (T/F): should feature activity levels be log2-transformed?
## "removeGrey" = logical (T/F): by WGCNA convention, unassigned features (e.g. genes) are denoted by the color "grey"; should "grey" features be excluded from analysis? default = TRUE;

## Note: ModuleSampleNetwork requires the following libraries: geneplotter, gtools, qvalue, and WGCNA;
## Note: ModuleSampleNetwork also requires the following R functions, which are included below:
## "ModulePrinComps1" (by Steve Horvath); "plot.mat" (by Sandrine Dudoit; slightly modified by Mike Oldham); and "MedianNegs" (by Mike Oldham)


library(geneplotter)
library(gtools)
library(qvalue)
library(WGCNA)


## ModuleSampleNetwork function (by Mike Oldham)

if(exists("ModuleSampleNetwork")) rm(ModuleSampleNetwork);
ModuleSampleNetwork=function(datExprT,method1="correlation",impute1=FALSE,skip1,indices1,modules1,featurelabels1,sampleinfo1,subgroup1,samplelabels1,grouplabels1,fitmodels1=FALSE,whichmodel1="univariate",whichfit1="pc1",btrait1=NULL,trait1=NULL,asfactors1=NULL,projectname1,cexlabels1=1,verbose=TRUE,minGpSize1=10,corType1="s",logT,removeGrey=TRUE){  

  ## Check for numeric data:
  checknumeric=c()
  for(e in c((skip1+1):length(datExprT[1,]))){
  	checknumeric=c(checknumeric,is.numeric(datExprT[,e]))
  	}
  if(all(checknumeric)!=TRUE){
  	stop("expression matrix contains non-numeric data")
  	}
  
  ## Enforce concordance between datExprT and sampleinfo1:
  rownames(sampleinfo1)=c(1:dim(sampleinfo1)[1])
  groups1=unique(sampleinfo1[,grouplabels1])
  if(length(indices1)[1]!=length(groups1)){
    stop("number of indices does not match number of groups")
	}
  maxindex=c()
  for(i in c(1:length(indices1)[1])){
    maxindex=max(maxindex,max(indices1[[i]]))
	}
  matchlabels=data.frame(dimnames(datExprT)[[2]][(skip1+1):maxindex],sampleinfo1[,samplelabels1])
  dimnames(matchlabels)[[2]]=c("datExprT","sampleinfo1")
  Match=matchlabels[,1]==matchlabels[,2]
  if(length(Match[Match])!=length(dimnames(datExprT)[[2]][(skip1+1):maxindex])){
	stop("Sample labels in datExprT and sampleinfo1 do not match!")
	}
  
  ## Enforce order of group indices listed in indices1 (these must be listed in the same order that groups first appear in datSample1, datExprT):
  gpordervec=c()
  for(f in c(1:length(indices1))){
    gpordervec=c(gpordervec,indices1[[f]][1])
	}
  rankorder=rank(gpordervec)
  indices1=indices1[rankorder]
  
  ## Enforce corType1 = "s" (spearman; default), "p" (pearson), or "k" (kendall); used when correlating Z.C vs Z.K values for each subgroup in each module.
  
  if(length(intersect(corType1,c("s","p","k")))!=1){
	stop("Error! corType1 must be 's', 'p', or 'k'")
	}
  
  ## Impute missing data:
  if(impute1==TRUE){
    missingtable=table(is.na(datExprT[,c((skip1+1):maxindex)]))
	if(is.na(missingtable[2])){
	  print("No missing data...no imputation")
	  } else {
	    print("Imputing missing data...")
	    datimpute=impute.knn(as.matrix(datExprT[,c((skip1+1):maxindex)]))
	    datExprT=data.frame(datExprT[,c(1:skip1)],datimpute$data)
	    }
	  }	    
  
  ## Exclude probe sets with 0 variance:
  excludevec=rep(1,length(datExprT[,1]))
  for(a in c(1:length(indices1)[1])){
    variance=apply(datExprT[,indices1[[a]]],1,var)
	restvar=variance==0
	excludevec[restvar]=0
	}
  excludevec=as.logical(excludevec)
  datExprT=datExprT[excludevec,]
  if(length(excludevec[excludevec])>0){
    print(paste("Note: ",length(excludevec[!excludevec])," probe sets had 0 variance and were excluded",sep=""))
	}
  collectGarbage()

   ## Enforce asfactors1, whichmodel1, and whichfit1:
  if(fitmodels1==TRUE){
  	if(length(asfactors1)>0){
      for(h in c(1:length(asfactors1))){
        sampleinfo1[,asfactors1[h]]=factor(sampleinfo1[,asfactors1[h]])
        }
      }
  	if(length(grep("univariate",whichmodel1))!=1&length(grep("multivariate",whichmodel1))!=1){
	  stop("whichmodel1 must equal 'univariate' or 'multivariate'")
	  }
	if(length(grep("K",whichfit1))!=1&length(grep("mean",whichfit1))!=1&length(grep("pc1",whichfit1))!=1){
  	  stop("whichfit1 must equal 'K', 'mean', or 'pc1'")
  	  }
  	}
  
   ## Create ModuleSampleNetwork root directory, if necessary:
	if(length(grep(paste(projectname1,"_ModuleSampleNetworks",sep=""),getwd()))==1){
	  breakdir=strsplit(getwd(),split="/")
	  SNroot=c(1:length(breakdir[[1]]))[is.element(breakdir[[1]],paste(projectname1,"_ModuleSampleNetworks",sep=""))]
	  setwd(paste(breakdir[[1]][1:SNroot],collapse="/"))
	  } else {
	  dir.create(paste(projectname1,"_ModuleSampleNetworks",sep=""))
	  setwd(paste(getwd(),"/",projectname1,"_ModuleSampleNetworks",sep=""))
	  }
	SNrootDir=getwd()	
  
  ## MakePlots2 function:
  MakePlots2=function(datexpr2,IAC2,adj2,sampleinfo2,meansample2,indexgp2,group2,subgroup2,btrait2,trait2,asfactors2,whichmodule2,fitmodels2,whichmodel2,whichfit2,verbose2){
	
	datexprgp=datexpr2[,indexgp2]
  	meansample2=apply(datexprgp,2,mean,na.rm=T)
	allsamples=c(1:length(indices1[[r]]))
	
	## For coloring subgroups in plots:
	subgpcolors=c("black","red","turquoise","blue","brown","darkgreen","orange","purple","salmon","tan","magenta","grey60","skyblue","darkred","lightgreen","midnightblue","greenyellow")
	colorvec=rep("black",length(datexprgp[1,]))
	if(!is.null(subgroup2)){
	  whichsubgroups=sort(unique(sampleinfo2[indexgp2-skip1,subgroup2]))
	  for(b in c(1:length(whichsubgroups))){
	    colorvec[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]=subgpcolors[b]
		}
	  }
	
	## To set colored leaf labels (adapted from dendrapply function help):
    local({
    colLab <<- function(n,treeorder) {
      if(is.leaf(n)) {
        a <- attributes(n)
        b <<- b+1
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colorvec[treeorder][b], lab.font= i%%3))
        }
      n
     }
    b <- 0
    })
	  
	FNC=fundamentalNetworkConcepts(adj2)
	K2=FNC$ScaledConnectivity
	pcs1all=ModulePrinComps1(datexpr=as.matrix(t(datexprgp)),couleur=rep("gold",length(datexprgp[,1])))
	pc2=pcs1all$PrinComps$PCgold
	meanIAC=mean(IAC[upper.tri(IAC)],na.rm=T)
    cluster1=hclust(as.dist(1-adj2),method="average")
	cluster1order=cluster1$order
	cluster2=as.dendrogram(cluster1,hang=0.1)
	cluster2=dendrapply(cluster2,colLab,cluster1order)
    Z.K=(K2-mean(K2))/sd(K2)
	Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
	
	GeneCor=cor(t(datexprgp),method="p",use="p")
	diag(GeneCor)=0
	GeneAdj=((1+GeneCor)/2)^2
	diag(GeneAdj)=0
	kin=apply(GeneAdj,1,sum)
		
	  sampleinfobtrait=data.frame(sampleinfo2[indexgp2-skip1,btrait2])
      onelevelcheckbtrait=c()
	  for(z in c(1:dim(sampleinfobtrait)[2])){
      	onelevelcheckbtrait=c(onelevelcheckbtrait,length(unique(sampleinfobtrait[,z]))>1)
      	}
	  if(!all(onelevelcheckbtrait)){
      	btrait2=btrait2[onelevelcheckbtrait]
      	}	
      if(length(trait2)>0){
	    sampleinfotrait=data.frame(sampleinfo2[indexgp2-skip1,trait2])
		onelevelchecktrait=c()
		for(y in c(1:dim(sampleinfotrait)[2])){
		  onelevelchecktrait=c(onelevelchecktrait,length(unique(sampleinfotrait[,y]))>1)
      	  }
		if(!all(onelevelchecktrait)){
		  trait2=trait2[onelevelchecktrait]
      	  }
		}
	  alltraits=union(btrait2,trait2)
      allmodelfactors=is.element(alltraits,asfactors2)
      allmodelterms=paste("sampleinfo2[indexgp2-skip1,",alltraits,"]",sep="")
      allmodelterms[allmodelfactors]=paste("factor(",allmodelterms[allmodelfactors],")",sep="")
      origmodelterms=allmodelterms
      allmodelterms=paste(allmodelterms,collapse="+")
     
	   if(whichfit2=="mean"){
        if(whichmodel2=="univariate"){
		  allsigniflist=c()
		  alltempout=c()
		  for(f in c(1:length(origmodelterms))){
		    allmodelname=paste("meansample2~",origmodelterms[f],sep="")
			alltemp=anova(lm(as.formula(allmodelname)))
			alltempout=c(alltempout,gsub(" ","",rownames(alltemp)[1]))
			allsigniflist=c(allsigniflist,alltemp[1:(dim(alltemp)[1]-1),5])
			}
		  allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
		  } else {
		  allmodelname=paste("meansample2~",allmodelterms,sep="")  
          alltemp=anova(lm(as.formula(allmodelname)))
		  alltempout=gsub(" ","",rownames(alltemp))
          allsigniflist=alltemp[1:(dim(alltemp)[1]-1),5]
          allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
          }
		bplabel="mean expression"
		outcomevec=meansample2
		}
		
      if(whichfit2=="K"){
      	if(whichmodel2=="univariate"){
		  allsigniflist=c()
		  alltempout=c()
		  for(f in c(1:length(origmodelterms))){
		    allmodelname=paste("K2~",origmodelterms[f],sep="")
			alltemp=anova(lm(as.formula(allmodelname)))
			alltempout=c(alltempout,gsub(" ","",rownames(alltemp)[1]))
			allsigniflist=c(allsigniflist,alltemp[1:(dim(alltemp)[1]-1),5])
			}
		  allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
		  } else {
		  allmodelname=paste("K2~",allmodelterms,sep="")  
		  alltemp=anova(lm(as.formula(allmodelname)))
		  alltempout=gsub(" ","",rownames(alltemp))
		  allsigniflist=alltemp[1:(dim(alltemp)[1]-1),5]
		  allsigniflist[is.nan(allsigniflist)]=1
		  lmoutput=lm(as.formula(allmodelname))
		  }
		bplabel="K"
		outcomevec=K2
		}	
		
      if(whichfit2=="pc1"){
      	if(whichmodel2=="univariate"){
		  allsigniflist=c()
		  alltempout=c()
		  for(f in c(1:length(origmodelterms))){
		    allmodelname=paste("pc2~",origmodelterms[f],sep="")
			alltemp=anova(lm(as.formula(allmodelname)))
			alltempout=c(alltempout,gsub(" ","",rownames(alltemp)[1]))
			allsigniflist=c(allsigniflist,alltemp[1:(dim(alltemp)[1]-1),5])
			}
		  allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
		  } else {
		  allmodelname=paste("pc2~",allmodelterms,sep="")  
      	  alltemp=anova(lm(as.formula(allmodelname)))
		  alltempout=gsub(" ","",rownames(alltemp))
          allsigniflist=alltemp[1:(dim(alltemp)[1]-1),5]
          allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
          }
		bplabel="pc1"
		outcomevec=pc2
		}
		     	    
      if(length(trait2)>0){
      varlist=vector(mode="list",length=length(trait2))
      for(q in c(1:length(trait2))){
          traitdf=data.frame(table(sampleinfo2[indexgp2-skip1,trait2[q]]))
		  varlist[[q]]=sort(unique(traitdf$Var1[traitdf$Freq>1]))
          }
      signiflist=vector(mode="list",length=length(trait2))
      for(j in c(1:length(trait2))){
        btrait3=setdiff(btrait2,trait2[j])
        whichmodelfactors=is.element(btrait3,asfactors2)
        modelterms=paste("sampleinfo2[indexgp2-skip1,",btrait3,"]",sep="")
        modelterms[whichmodelfactors]=paste("factor(",modelterms[whichmodelfactors],")",sep="")
        modelterms=paste("+",paste(modelterms,collapse="+"),sep="")
        if(length(btrait3)==0|whichmodel2=="univariate"){
        	modelterms=NULL
        	}
        whichtestfactors=is.element(trait2,asfactors2)
	
		if(whichfit2=="mean"){
          if(whichtestfactors[j]==TRUE){
            for(k in c(1:length(varlist[[j]]))){
              modelname=paste("meansample2~factor(sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"])",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
            } else {
      	      modelname=paste("meansample2~sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"]",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
			}
        
		if(whichfit2=="K"){
          if(whichtestfactors[j]==TRUE){
            for(k in c(1:length(varlist[[j]]))){
              modelname=paste("K2~factor(sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"])",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
            } else {
      	      modelname=paste("K2~sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"]",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
			}
        
		if(whichfit2=="pc1"){
          if(whichtestfactors[j]==TRUE){
            for(k in c(1:length(varlist[[j]]))){
              modelname=paste("pc2~factor(sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"])",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
            } else {
      	      modelname=paste("pc2~sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"]",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
			}
            signiflist[[j]][is.nan(signiflist[[j]])]=1
          } ## end of for(j in c(1:length(trait2))){
        }  ## end of if(length(trait2)>0){
		
	  if(verbose2==TRUE){
	    plotrows=3
		plotcols=ceiling((length(trait2)+5)/3)
		} else {
		  plotcols=ceiling((length(trait2)+3)/3)
		  if(length(trait2)>1){
		    plotrows=3
		    } else {
	          plotrows=2
			  }
		    }
		pdf(file=paste(projectname1,"_",group2,"_x_",colnames(sampleinfo2)[subgroup2],"_SampleNetwork_",whichmodule2,".pdf",sep=""))
		par(mfrow=c(plotrows,plotcols))
        par(mar=c(5,5,4,2))
		plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(adj2[upper.tri(adj2)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
		mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
		par(mar=c(5,5,4,2))
		plot(Z.K,main="Connectivity",ylab="Z.K",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
        text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
        abline(h=-2)
		abline(h=-3)
		if(verbose2==TRUE){
	      par(mar=c(5,5,4,2))
		  plot(Z.C,main="ClusterCoef", ylab="Z.C",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
	      text(Z.C,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
		  abline(h=2)
		  abline(h=3)
		  abline(h=-2)
		  abline(h=-3)
		  par(mar=c(5,5,4,2))
		  plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec,cex.main=1.8,cex.lab=1.4)				
		  if(!is.null(subgroup2)){
			if(length(whichsubgroups)==2){
				for(b in c(1:length(whichsubgroups))){
					abline(lm(Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]~Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]),col=subgpcolors[b],lwd=2)
				}
				mtext(paste("rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$p.value,2),
							"; ","rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$p.value,2),
							sep=""),cex=0.8,line=0.2)
			} else {
				abline(lm(Z.C~Z.K),col="black",lwd=2)
				mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
			}
		  }
		}
			
        if(length(intersect(origmodelterms,alltempout))!=length(origmodelterms)){
        alltraitsbp=alltraits[is.element(origmodelterms,alltempout)]
        } else {
          alltraitsbp=alltraits
          }
		allsigniflist[allsigniflist==0]=1e-300
		par(mar=c(8,5,4,2))
		barplot(-log(allsigniflist,10),names=as.character(dimnames(sampleinfo2)[[2]][alltraitsbp]),las=3,cex.names=0.8,ylab="-Log10 p-value",main=paste("ANOVA of ",bplabel,sep=""),cex.main=1.8,cex.lab=1.4)
        mtext(whichmodel2,cex=0.8,line=0.2)
		abline(h=-log(.05,10),col="blue",lwd=2)
        abline(h=-log(.05/length(alltraitsbp),10),col="red",lwd=2)
		names=dimnames(datexprgp)[[2]]
		ordergenes=order(-kin)
		ordersamples=order(sampleinfo2[,subgroup2])
		datexprgpHM=scale(t(datexprgp[ordergenes,ordersamples]))
		par(mar=c(8,5,4,2))
		plot.mat(t(datexprgpHM),clabels=names[ordersamples],cex.main=1.8,cex.lab=1.4,ccols=colorvec[ordersamples])
		#datexprgpHM=scale(t(datexprgp[ordergenes,]))
		#plot.mat(t(datexprgpHM),clabels=names,cex.main=1.8,cex.lab=1.4,ccols=colorvec)
		title(main="Gene expression",cex.main=1.8)
						
		if(length(trait2)>0){
          for(m in c(1:length(signiflist)[[1]])){
            varlabel=dimnames(sampleinfo2)[[2]][trait2][m]  
            signiflist[[m]][signiflist[[m]]==0]=1e-300
			par(mar=c(8,5,4,2))
			plot(-log(signiflist[[m]],10),type="n",xlab=varlabel,ylab="-Log10 p-value",main=paste(group2," ",whichmodule2," ANOVA of ",bplabel,sep=""),xaxt="n",cex.main=1.8,cex.lab=1.4)
            mtext(whichmodel2,cex=0.8,line=0.2)
			text(-log(signiflist[[m]],10),labels=varlist[[m]])
            abline(h=-log(.05,10),col="blue",lwd=2)
            abline(h=-log(.05/length(varlist[[m]]),10),col="red",lwd=2)
            }
          }
        dev.off()
	  collectGarbage()
	  cat("\n")
          
		    if(length(trait2)>0){
			  if(!all(onelevelchecktrait)){
			    print(paste("Warning:",dimnames(sampleinfotrait)[[2]][!onelevelchecktrait],"contains only one level and was excluded as a trait for model-fitting"))
      	        }
			  }
      	    if(!all(onelevelcheckbtrait)){
      	      print(paste("Warning:",dimnames(sampleinfobtrait)[[2]][!onelevelcheckbtrait],"contains only one level and was excluded as a btrait for model-fitting"))
      	      }	
            if(length(lmoutput$coefficients[is.na(lmoutput$coefficients)])>0){
              print(paste("Warning: the following coefficient was not defined because of a singularity:",names(lmoutput$coefficients)[is.na(lmoutput$coefficients)]))
              }
		modmetrics=data.frame(signif(meanIAC,3),signif(mean(FNC$Connectivity),3),signif(mean(FNC$ScaledConnectivity),3),signif(mean(FNC$ClusterCoef),3),signif(mean(FNC$MAR),3),signif(FNC$Density,3),signif(1-FNC$Centralization,3),signif(1-FNC$Heterogeneity,3),signif(pcs1all$varexplained[1,1],3))
		emptydf=data.frame()
		attr(emptydf,"ModMetrics")=modmetrics
		attr(emptydf,"Outcome")=outcomevec
		attr(emptydf,"Z.Cout")=Z.C
		attr(emptydf,"Z.Kout")=Z.K
		attr(emptydf,"WhichSubGroupsout")=whichsubgroups
		emptydf
	  } ## end of MakePlots2 fx
	
  for(r in c(1:length(indices1)[1])){

	## Create output directory:
	tstamp=format(Sys.time(), "%X")
	tstamp=gsub(":","-",tstamp)
	dir.create(paste(groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",tstamp,sep=""))
	setwd(paste(groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",tstamp,sep=""))
	bymodwd=getwd()

  ## Initialize data frame for export:
  datout=data.frame()
  modmetricsout=data.frame()
  sampleinfomods=sampleinfo1
  
  ## Subset expression data by module:
  if(dim(modules1)[1]!=dim(datExprT)[1]){
	cat("\n")
	cat("\n")
	print("WARNING! Number of rows in modules1 does not equal number of rows in datExprT.")
	cat("\n")
	cat("\n")
	}
  allmods=sort(unique(modules1[,2][!is.na(unique(modules1[,2]))]))
  if(removeGrey==TRUE){
    allmods=allmods[allmods!="grey"]
	}
	  
  dir.create("Sample_network_plots")
  setwd(paste(getwd(),"/Sample_network_plots",sep=""))
    
  MEdf=data.frame(c(1:dim(sampleinfo1)[1]))
  MEdf=MEdf[,-c(1)]
  
  for(m in c(1:length(allmods))){
    
	modpsid=as.character(modules1[!is.na(modules1[,2])&modules1[,2]==allmods[m],1])
	restmod=is.element(datExprT[,1],modpsid)
	datExprMod=datExprT[restmod,]
    print("Building sample network...")
    IAC=cor(datExprMod[,indices1[[r]]],method="p",use="p")
    diag(IAC)=0
	if(method1=="correlation"){
	  A.IAC=((1+IAC)/2)^2
	  }
	if(method1=="euclidean"){
	  D.squared=as.matrix(dist(t(datexprgp),method="euclidean")^2)
	  A.IAC=1-D.squared/max(D.squared)
	  }
    diag(A.IAC)=0
    print(paste(groups1[r]," ",allmods[m]," module",sep=""))  
	datoutnew=MakePlots2(datexpr2=datExprMod,IAC2=IAC,adj2=A.IAC,sampleinfo2=sampleinfo1,indexgp2=indices1[[r]],group2=groups1[r],subgroup2=subgroup1,btrait2=btrait1,trait2=trait1,asfactors2=asfactors1,whichmodule2=allmods[m],fitmodels2=fitmodels1,whichmodel2=whichmodel1,whichfit2=whichfit1,verbose2=verbose)
	Z.C=attr(datoutnew,"Z.Cout")
	Z.K=attr(datoutnew,"Z.Kout")
	if(length(intersect(subgroup1,asfactors1))>0){
	  OutcomeCorP=cor.test(attr(datoutnew,"Outcome"),as.numeric(factor(sampleinfo1[indices1[[r]]-skip1,subgroup1])))$p.value
	  allmodelname=paste("Z.C~Z.K*factor(sampleinfo1[indices1[[r]]-skip1,subgroup1])",sep="")
	  alltemp=coef(summary(lm(as.formula(allmodelname))))
	  } else {
	  OutcomeCorP=cor.test(attr(datoutnew,"Outcome"),sampleinfo1[indices1[[r]]-skip1,subgroup1])$p.value
	  allmodelname=paste("Z.C~Z.K*sampleinfo1[indices1[[r]]-skip1,subgroup1]",sep="")
	  alltemp=coef(summary(lm(as.formula(allmodelname))))
	  }
	
	mod1=data.frame(cbind(OutcomeCorP,t(alltemp[2:dim(alltemp)[1],dim(alltemp)[2]])))
	datout=data.frame(rbind(datout,mod1))
	modmetricsout=data.frame(rbind(modmetricsout,attr(datoutnew,"ModMetrics")))
	sampleinfomods=data.frame(cbind(sampleinfomods,Z.K,Z.C))
	colnames(sampleinfomods)[dim(sampleinfomods)[2]-1]=paste("Z.K_",allmods[m],sep="")
	colnames(sampleinfomods)[dim(sampleinfomods)[2]]=paste("Z.C_",allmods[m],sep="")
	MEdf=data.frame(MEdf,attr(datoutnew,"Outcome"))
	} ## end of for(m in c(1:length(allmods))){

    if(length(intersect(subgroup1,asfactors1))>0){
	  termlabels=rownames(alltemp)[2:dim(alltemp)[1]]
	  termlabels=gsub("factor(sampleinfo1[indices1[[r]] - 1, subgroup1])",colnames(sampleinfo1)[subgroup1],termlabels,fixed=TRUE)
	  } else {
	  termlabels=rownames(alltemp)[2:dim(alltemp)[1]]
	  termlabels=gsub("sampleinfo1[indices1[[r]] - 1, subgroup1]",colnames(sampleinfo1)[subgroup1],termlabels,fixed=TRUE)
	  }

	datout=signif(datout,3)
	datout=data.frame(as.character(allmods),datout)
	colnames(datout)=c("Module","MEcorP",termlabels)
	modmetricsout=data.frame(as.character(allmods),modmetricsout)
	colnames(modmetricsout)=c("Module","Mean_IAC","Mean_Connectivity","Mean_ScaledConnectivity","Mean_ClusterCoef","Mean_MAR","Density","Decentralization","Homogeneity","PC1_VE")
		
	## Export network metrics by module (modmetricsout) and sample (sampleinfomods):
	setwd(bymodwd)
	write.table(modmetricsout,file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_network_metrics_all_samples.csv",sep=""),sep=",",row.names=F)
	write.table(sampleinfomods,file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_network_metrics_by_sample.csv",sep=""),sep=",",row.names=F)
	
	## Export linear model analysis:
	
	print("Building linear models...")
	dir.create("Linear_models")
	setwd(paste(getwd(),"/Linear_models",sep=""))
		
	write.table(datout,file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_LM.csv",sep=""),sep=",",row.names=F)
	
	noplots=dim(datout)[2]-1
	if(noplots>4){
	    newplotrows=3
		newplotcols=ceiling(noplots/3)
		} else {
		newplotrows=2
		newplotcols=2
		}
		
	pdf(file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_LM.pdf",sep=""))
#par(mfrow=c(newplotrows,newplotcols))
	par(mar=c(7,5,4,2))
	  datout[,2][datout[,2]==0]=NA
	  barplot1=barplot(c(-log(datout[,2],10)),col=as.character(datout[,1]),names.arg=as.character(datout[,1]),las=3,cex.axis=1.25,ylab="-Log10 p-value",main=paste("cor(",colnames(sampleinfo1)[subgroup1],",ME)",sep=""),cex.main=1.8,cex.lab=1.4)
	  text(barplot1[is.na(datout[,2])],0.5,"INF",srt=90)
	  abline(h=-log(.05,10),col="blue",lwd=2)
	  abline(h=-log(.05/length(allmods),10),col="red",lwd=2)
	  for(i in c(3:dim(datout)[2])){
	    datout[,i][datout[,i]==0]=NA
		par(mar=c(7,5,4,2))
		barplot2=barplot(c(-log(datout[,i],10)),col=as.character(datout[,1]),names.arg=as.character(datout[,1]),las=3,cex.axis=1.25,ylab="-Log10 p-value",main=paste("lm(Z.C~Z.K*",colnames(sampleinfo1)[subgroup1],")",sep=""),cex.main=1.8,cex.lab=1.4)
	    text(barplot2[is.na(datout[,i])],0.5,"INF",srt=90)
		mtext(colnames(datout)[i],line=0.2,cex=0.8)
		abline(h=-log(.05,10),col="blue",lwd=2)
	    abline(h=-log(.05/length(allmods),10),col="red",lwd=2)
		}
	dev.off()
	
	setwd(bymodwd)
	
	whichsubgroupsout=attr(datoutnew,"WhichSubGroupsout")
	gpcounts=data.frame(table(sampleinfo1[,subgroup1]))
	keptgps=sort(factor(intersect(whichsubgroupsout,gpcounts$Var1[gpcounts$Freq>=minGpSize1])))
	print(data.frame(keptgps))
	print("The groups listed above will be retained. Please place them in the desired order for figures using the row numbers listed above (e.g. 2,1,3).",quote=F)
	answer=readline(prompt="Response: ")
	while(answer==""){
  	  answer=readline(prompt="Response: ")
  	  }
	answer=strsplit(answer,split=",")
	keptgps=keptgps[as.numeric(answer[[1]])]
	
	## Need to make above more robust eventually.
	
	ncomp=(length(keptgps)*(length(keptgps)-1))/2
	seq1=seq(dim(sampleinfo1)[2]+1,dim(sampleinfomods)[2],by=2)
	seq2=seq(dim(sampleinfo1)[2]+2,dim(sampleinfomods)[2],by=2)
	
	## Create a matrix for tracking Z.C vs Z.K correlations, mean Z.C, and p-values:
	
	modsignifout=matrix(nrow=length(allmods),ncol=((length(keptgps)*3)+((length(keptgps)*(length(keptgps)-1)))+3),data=-88888)
		
	## Create a matrix with all possible pairwise combinations of kept groups:
	
	combmatrix=combinations(length(keptgps),2,as.character(keptgps))
	
	## Then cycle through mods:
		
	print("Comparing correlations of Z.C and Z.K...") 
	for(i in c(1:length(allmods))){
		modsignifout[i,1]=cor.test(sampleinfomods[,seq1[i]],sampleinfomods[,seq2[i]],method=corType1)$estimate
		modsignifout[i,2]=cor.test(sampleinfomods[,seq1[i]],sampleinfomods[,seq2[i]],method=corType1)$p.value
		modsignifout[i,(length(keptgps)*2)+dim(combmatrix)[1]+3]=mean(sampleinfomods[,seq2[i]])
		for(j in c(1:length(keptgps))){
			modsignifout[i,1+(2*j)]=cor.test(sampleinfomods[as.character(sampleinfomods[,subgroup1])==keptgps[j],seq1[i]],sampleinfomods[as.character(sampleinfomods[,subgroup1])==keptgps[j],seq2[i]],method=corType1)$estimate
			modsignifout[i,2+(2*j)]=cor.test(sampleinfomods[as.character(sampleinfomods[,subgroup1])==keptgps[j],seq1[i]],sampleinfomods[as.character(sampleinfomods[,subgroup1])==keptgps[j],seq2[i]],method=corType1)$p.value
			modsignifout[i,(length(keptgps)*2)+dim(combmatrix)[1]+3+j]=mean(sampleinfomods[as.character(sampleinfomods[,subgroup1])==keptgps[j],seq2[i]])
			}
		for(k in c(1:ncomp)){
			modsignifout[i,((length(keptgps)*3)+((length(keptgps)*(length(keptgps)-1)/2))+3+k)]=wilcox.test(sampleinfomods[as.character(sampleinfomods[,subgroup1])==combmatrix[k,1],seq2[i]],sampleinfomods[as.character(sampleinfomods[,subgroup1])==combmatrix[k,2],seq2[i]])$p.value
			}
		}

	## Calculate significance of difference between group correlations:
	
	rtoz=function(x){0.5*log((1+x)/(1-x))}
	
	seqcors=seq(3,(length(keptgps)*2)+2,by=2)
	seqsignif=seq(3+length(keptgps)*2,2+length(keptgps)*2+dim(combmatrix)[1],by=1)
	
	zmat=matrix(nrow=length(allmods),ncol=length(keptgps),data=-88888)
	
	for(i in c(1:length(allmods))){
		for(j in c(1:length(keptgps))){
			zmat[i,j]=rtoz(modsignifout[i,seqcors[j]])
			}
		}
	colnames(zmat)=keptgps
	keptgpcounts=gpcounts[is.element(gpcounts[,1],keptgps),]
	for(k in c(1:dim(combmatrix)[1])){
		gp1=combmatrix[k,1]
		gp2=combmatrix[k,2]
		n1=keptgpcounts$Freq[keptgpcounts$Var1==gp1]
		n2=keptgpcounts$Freq[keptgpcounts$Var1==gp2]
		SEcorr=sqrt((1/(n1-3))+(1/(n2-3)))
		zdiff=(zmat[,colnames(zmat)==gp1]-zmat[,colnames(zmat)==gp2])/SEcorr
		PcorDiff=2*pnorm(-abs(zdiff))
		modsignifout[,seqsignif[k]]=PcorDiff
		}
		
	colnames(modsignifout)=c("Cor.CCvsK.All","Pval.CCvsK.All",paste(rep(c("Cor.CCvsK.","Pval.CCvsK."),length(keptgps)),rep(keptgps,each=2),sep=""),paste("Pval.",combmatrix[,1],"cor.vs.",combmatrix[,2],"cor",sep=""),"Mean.CC.All",paste("Mean.CC.",keptgps,sep=""),paste("Pval.",combmatrix[,1],"CC.vs.",combmatrix[,2],"CC",sep=""))
	modsignifout=data.frame(allmods,modsignifout)
	colnames(modsignifout)[1]="Module"
	
	if(corType1=="s"){
		  whichcor="Spearman"
		  }
		if(corType1=="p"){
		  whichcor="Pearson"
		  }
		if(corType1=="k"){
		  whichcor="Kendall"
		  }
	
	combvecmatrix=combinations(length(keptgps),2,c(1:length(keptgps)))
	
	## Export eigensamples and contrast plots:
	
	setwd(bymodwd)
	dir.create("Module_eigensamples")
	setwd(paste(getwd(),"/Module_eigensamples",sep=""))
	
	ValueToColor=function(x, min1=-1, max1=1, breaks1=20){
		MasterKMEcols=data.frame(seq(min1,max1,by=(max1-min1)/breaks1), GetColor(seq(min1,max1,by=(max1-min1)/breaks1),GreenRed=TRUE,DisplayRange=1))
		colnames(MasterKMEcols)=c("RoundeKME","Color")
		rounded=round(x,digits=1)
		roundedcols=as.character(MasterKMEcols[match(as.character(rounded),as.character(MasterKMEcols[,1])),2])
		legend1=barplot(rep(1,dim(MasterKMEcols)[1]),col=rev(as.character(MasterKMEcols[,2])),space=0,yaxt="n",cex.main=1.8,cex.lab=1.4)
		axis(side=1,at=legend1,labels=rev(MasterKMEcols[,1]))
		as.character(roundedcols)
		}
	
	print(keptgps)
	print("Please choose reference group for eigensample contrast plots.",quote=F)
	answer=readline(prompt="Response: ")
	while(answer==""|length(intersect(answer,keptgps))==0){
  	  answer=readline(prompt="Response: ")
  	  }
	reference=answer
	refgp=c(1:length(keptgps))[is.element(keptgps,reference)]
	nonrefgps=c(1:length(keptgps))[!is.element(keptgps,reference)]
	
	for(i in c(1:length(allmods))){
	  whichmod=allmods[i]
	  restmodrow=is.element(modules1[,2],whichmod)
	  datExprRestT=t(datExprT[restmodrow,indices1[[r]]])
	  NewME=ModulePrinComps1(datexpr=datExprRestT,couleur=rep(whichmod,dim(datExprRestT)[2]))
	  NewMEmod=NewME$PrinComps[,1]
	  modkmenew=c()
	  for(j in c(1:dim(datExprRestT)[2])){
			modkmenew=c(modkmenew,cor.test(NewMEmod,datExprRestT[,j])$estimate)
			}
	  modkmenewrank=order(-modkmenew)
	  datExprRest=t(datExprRestT)
	  SElist=vector(mode="list",length=length(keptgps))
	  for(k in c(1:length(keptgps))){
		TempSE=ModulePrinComps1(datexpr=datExprRest[,sampleinfo1[,subgroup1]==as.character(keptgps[k])],couleur=rep(whichmod,length(sampleinfo1[,subgroup1][sampleinfo1[,subgroup1]==as.character(keptgps[k])])))
		SElist[[k]]=TempSE$PrinComps[,1]
		}
	  ContrastList=vector(mode="list",length=length(nonrefgps))
	  for(m in c(1:length(nonrefgps))){
	    ContrastList[[m]]=as.numeric(SElist[[nonrefgps[m]]][modkmenewrank])-as.numeric(SElist[[refgp]][modkmenewrank])
		}
	  modGenes=as.character(featurelabels1[restmodrow])
	  
	  pdf(file=paste(whichmod,"_eigensamples_by_",colnames(sampleinfo1)[subgroup1],".pdf",sep=""),width=21,height=12)
	  nf = layout(matrix(c(1:(length(keptgps)+1)),(length(keptgps)+1),1,byrow=TRUE), widths=c(1), heights=c(0.4,rep(1,length(keptgps))), respect=FALSE)
	  par(pin=c(20,.1))
	  barColors=ValueToColor(modkmenew[modkmenewrank], min1=-1, max1=1, breaks1=20)
	  mtext(text="kME",side=3,line=.2,cex=1.5)
	  for(n in c(1:length(keptgps))){
	    par(pin=c(20,(10/length(keptgps))))
	    barplot(SElist[[n]][modkmenewrank],col=barColors,cex.main=1.8,cex.lab=1.4)
	    mtext(text=paste(whichmod," module eigensample (",colnames(sampleinfo1)[subgroup1]," ",keptgps[n],")",sep=""),side=3,line=-2,cex=1.5)
	    }
	  dev.off()
	 
	  ## Contrast:

	  meanComp=c()
	  sdComp=c()
	  nosd=2
	  restGene=rep(FALSE,length(ContrastList[[1]]))
	  
	  if(length(ContrastList[[1]])>200){
	    for(p in c(1:length(nonrefgps))){
		  meanComp=c(meanComp,mean(ContrastList[[p]]))
		  sdComp=c(sdComp,sd(ContrastList[[p]]))
	      }
		for(q in c(1:length(nonrefgps))){
		  restGene[ContrastList[[q]]>(meanComp[q]+nosd*sdComp[q])|ContrastList[[q]]<(meanComp[q]-nosd*sdComp[q])]=TRUE
		  }
		} else {
		restGene=rep(TRUE,length(ContrastList[[1]]))
		}
	  
	  SEdf=data.frame(c(1:length(modGenes)))
	  SEdf=SEdf[,-c(1)]
	  for(t in c(1:length(ContrastList)[[1]])){
	    SEdf=data.frame(cbind(SEdf,ContrastList[[t]]))
		}
	  
	  pdf(file=paste(whichmod,"_eigensample_contrasts_by_",colnames(sampleinfo1)[subgroup1],".pdf",sep=""),width=21,height=12)
	  par(mar=c(7,5,4,2))
	  matplot(SEdf,type="l",lty=1,lwd=1.5,ylab="Difference between eigensamples",xaxt="n",cex.lab=1.4,main=paste(whichmod,"module"),cex.main=1.8)
	  abline(v=c(1:length(modGenes)),col="darkgrey")
	  axis(side=1,at=c(1:length(modGenes))[restGene],labels=as.character(modGenes[modkmenewrank][restGene]),las=3,cex.axis=1.2)
	  legend("topleft",legend=paste(keptgps[nonrefgps]," - ",keptgps[refgp],sep=""),col=1:length(nonrefgps),lwd=2)
	  dev.off()
	  print(allmods[i])
	  }  ## end of for(i in c(1:length(allmods))){
	
	## Export analysis of differences in sample network concepts between each pair of groups, by module:
	
	print("Comparing sample network concepts...")
	setwd(bymodwd)
	write.table(modsignifout,file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_DSNC_summary.csv",sep=""),sep=",",row.names=F)
	
	seqmeans=c(((length(keptgps)*2)+dim(combmatrix)[1]+5):(dim(modsignifout)[2]-dim(combmatrix)[1]))
	seqmeansignif=c(((dim(modsignifout)[2]-dim(combmatrix)[1])+1):dim(modsignifout)[2])
					
	for(i in c(1:dim(combvecmatrix)[1])){
		pdf(file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",keptgps[as.numeric(combvecmatrix[i,1])],"_vs_",keptgps[as.numeric(combvecmatrix[i,2])],"_DSNC_summary.pdf",sep=""))
		par(mfrow=c(2,2))
		par(mar=c(5,4.5,4,2))
		plot(modsignifout[,seqcors[combvecmatrix[i,1]]+1],modsignifout[,seqcors[combvecmatrix[i,2]]+1],pch=19,col=as.character(modsignifout[,1]),xlim=c(-1,1),ylim=c(-1,1),xlab=paste("cor(K,C)",keptgps[as.numeric(combvecmatrix[i,1])]),ylab=paste("cor(K,C)",keptgps[as.numeric(combvecmatrix[i,2])]),cex.axis=1.25,cex.lab=1.4,main=paste("cor(K,C) by ",colnames(sampleinfomods)[subgroup1],sep=""),cex.main=1.6)
		abline(0,1)
		mtext(whichcor,line=0.2,cex=0.8)
		plotpvalues=c(-log(modsignifout[,intersect(grep(paste(keptgps[as.numeric(combvecmatrix[i,1])],"cor",sep=""),colnames(modsignifout)),grep(paste(keptgps[as.numeric(combvecmatrix[i,2])],"cor",sep=""),colnames(modsignifout)))],10))
		plotpvalues[abs(plotpvalues)=="Inf"]=NA
		par(mar=c(5,4.5,4,2))
		barplot3=barplot(plotpvalues,col=as.character(modsignifout[,1]),names.arg=as.character(modsignifout[,1]),las=3,cex.axis=1.25,ylab="-Log10 p-value",main="p.Diff.cor(K,C)",cex.main=1.6,cex.lab=1.4)
		text(barplot3[is.na(plotpvalues)],0.5,"INF",srt=90)
		mtext(paste(keptgps[as.numeric(combvecmatrix[i,1])]," vs ",keptgps[as.numeric(combvecmatrix[i,2])]),line=0.2,cex=0.8)
		abline(h=-log(.05,10),col="blue",lwd=2)
		abline(h=-log((.05/length(allmods)),10),col="red",lwd=2)
		par(mar=c(7,4.5,4,2))
		barplot4=barplot(c(-log(datout[,2],10)),col=as.character(datout[,1]),names.arg=as.character(datout[,1]),las=3,cex.axis=1.25,ylab="-Log10 p-value",main="p.DE",cex.main=1.6,cex.lab=1.4)
		text(barplot4[is.na(datout[,2])],0.5,"INF",srt=90)
		mtext(paste("cor(",colnames(sampleinfo1)[subgroup1],",ME)",sep=""),line=0.2,cex=0.8)
		abline(h=-log(.05,10),col="blue",lwd=2)
		abline(h=-log(.05/length(allmods),10),col="red",lwd=2)
		par(mar=c(7,4.5,4,2))
		plot(c(-log(datout[,2],10)),plotpvalues,pch=19,col=as.character(modsignifout[,1]),xlab="-Log10 p.DE",ylab="-Log10 p.Diff.cor(K,C)",cex.axis=1.25,cex.lab=1.4,main="Module summary",cex.main=1.6)
		mtext(paste("r = ",signif(cor.test(-log(datout[,2],10),plotpvalues)$estimate,2),"; p = ",signif(cor.test(-log(datout[,2],10),plotpvalues)$p.value,2),sep=""),line=0.2,cex=0.8)
		abline(lm(plotpvalues~c(-log(datout[,2],10))),col="black")
		abline(h=-log(.05/length(allmods),10),col="red",lwd=2)
		abline(v=-log(.05/length(allmods),10),col="red",lwd=2)
		dev.off()
		}
			
	## Export DE and DC analysis:
	
	print("Comparing differential expression and differential connectivity...")
	setwd(bymodwd)
	dir.create("DE_and_DC")
	setwd(paste(getwd(),"/DE_and_DC",sep=""))
	DEandDCwd=getwd()
	
	## Replacing negative values and log transforming for standard screening:
	
	datExprLogT=datExprT[,indices1[[r]]]
	if(logT==TRUE&length(datExprLogT[datExprLogT<0])>0){
		datExprLogT=MedianNegs(datE=datExprLogT)
		} else {
		datExprLogT=datExprLogT
		}
	if(logT==TRUE){
		datExprLogT=t(log(datExprLogT,2))
		}
			
	for(i in c(1:dim(combvecmatrix)[1])){
		dir.create(paste(combmatrix[i,1],"_vs_",combmatrix[i,2],sep=""))
		setwd(paste(getwd(),"/",combmatrix[i,1],"_vs_",combmatrix[i,2],sep=""))
		DEandDCgpsd=getwd()
		OutcomeCorPsubGp=c()
		Allkme1=rep(-88888,length(datExprT[,1]))
		Allkme2=rep(-88888,length(datExprT[,1]))
		restGps=is.element(sampleinfomods[,subgroup1],c(combmatrix[i,1],combmatrix[i,2]))
		print("Performing standard screening...")
		SSBT=standardScreeningBinaryTrait(datE=datExprLogT[restGps,],y=sampleinfomods[restGps,subgroup1],kruskalTest=F)
		SSBT=data.frame(datExprT[,1:skip1],as.character(modules1[,2]),SSBT[,2:dim(SSBT)[2]])
		## Replacing qval=0 with qval=2.2e-16:
		SSBT[,skip1+4][SSBT[,skip1+4]==0]=2.2e-16
		colnames(SSBT)[skip1+1]="Module"
		allmodkmecors=c()
		for(j in c(1:length(allmods))){
			OutcomeCorPsubGp=c(OutcomeCorPsubGp,cor.test(MEdf[restGps,j],as.numeric(factor(sampleinfomods[restGps,subgroup1])))$p.value)
			OutcomeCorPsubGp[OutcomeCorPsubGp==0]=NA
			whichmod=allmods[j]
			restmodrow=is.element(modules1[,2],whichmod)
			datExprRestT=t(datExprT[restmodrow,indices1[[r]]])
			NewME1=ModulePrinComps1(datexpr=datExprRestT[sampleinfo1[,subgroup1]==as.character(combmatrix[i,1]),],couleur=rep(whichmod,dim(datExprRestT)[2]))
			NewME1mod=NewME1$PrinComps[,1]
			NewME2=ModulePrinComps1(datexpr=datExprRestT[is.element(sampleinfo1[,subgroup1],c(as.character(combmatrix[i,1]),as.character(combmatrix[i,2]))),],couleur=rep(whichmod,dim(datExprRestT)[2]))
			NewME2mod=NewME2$PrinComps[,1]
			datExprRestT1=datExprRestT[sampleinfo1[,subgroup1]==as.character(combmatrix[i,1]),]
			datExprRestT2=datExprRestT[is.element(sampleinfo1[,subgroup1],c(as.character(combmatrix[i,1]),as.character(combmatrix[i,2]))),]
			## Need to make sure sign of ME is not arbitrarily flipped beween groups for modules with near-equal numbers of + and - correlated genes:
			corT1=cor(datExprRestT1)
			corT2=cor(datExprRestT2)
			corveccor=cor(corT1[upper.tri(corT1)],corT2[upper.tri(corT2)])
			## Note: will calculate kme1 using all samples from group 1 and kme2 using all samples (combined) from groups 1 and 2.
			## Note: will calculate kme values using original (i.e. non-log-transformed) data.
			for(k in c(1:dim(datExprRestT)[2])){
				Allkme1[restmodrow][k]=cor.test(NewME1mod,datExprRestT1[,k])$estimate
				Allkme2[restmodrow][k]=cor.test(NewME2mod,datExprRestT2[,k])$estimate
				}
			if(cor(Allkme1[restmodrow],Allkme2[restmodrow])<0&corveccor>0){
				Allkme2[restmodrow]=Allkme2[restmodrow]*-1
				}
			allmodkmecors=c(allmodkmecors,cor.test(Allkme1[restmodrow],Allkme2[restmodrow])$estimate)
			print(combmatrix[i,])
			print(allmods[j])
			}  ## end of for(j in c(1:length(allmods))){
		
		Allkme1[Allkme1==-88888]=NA
		Allkme2[Allkme2==-88888]=NA
		SSBT=data.frame(SSBT,Allkme1,Allkme2)
		colnames(SSBT)[11]=paste("kme_",as.character(combmatrix[i,1]),sep="")
		colnames(SSBT)[12]=paste("kme_",as.character(combmatrix[i,1]),"_and_",as.character(combmatrix[i,2]),sep="")
		Z1=rtoz(SSBT[,11])
		Z2=rtoz(SSBT[,12])
		n1=length(NewME1mod)
		n2=length(NewME2mod)
		SEcorr=sqrt((1/(n1-3))+(1/(n2-3)))
		Zdiff=(Z1-Z2)/SEcorr
		DCpval=2*pnorm(-abs(Zdiff))
		DCqval=rep(NA,length(DCpval))
		DCqval[!is.na(DCpval)]=qvalue(DCpval[!is.na(DCpval)])$qvalues
		## Replacing qval=0 with qval=2.2e-16:
		DCqval[DCqval==0]=2.2e-16
		SSBT=data.frame(SSBT,DCpval,DCqval)
		write.table(SSBT,file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_SSBT_",combmatrix[i,1],"_vs_",combmatrix[i,2],".csv",sep=""),sep=",",row.names=F,col.names=T)
		
		## Export DE_vs_DC summary:
		
		pdf(file=paste(groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",combmatrix[i,1],"_vs_",combmatrix[i,2],"_DE_vs_DC_summary.pdf",sep=""))
		par(mfrow=c(2,2))
		par(mar=c(7,5,4,2))
		barplot5=barplot(c(-log(OutcomeCorPsubGp,10)),col=as.character(allmods),names.arg=as.character(allmods),las=3,cex.axis=1.25,ylab="-Log10 p-value",main=paste("cor(",combmatrix[i,1]," & ",combmatrix[i,2],",ME)",sep=""),cex.main=1.8,cex.lab=1.4)
		text(barplot5[is.na(OutcomeCorPsubGp)],0.5,"INF",srt=90)
		abline(h=-log(.05,10),col="blue",lwd=2)
		abline(h=-log(.05/length(allmods),10),col="red",lwd=2)
		par(mar=c(7,5,4,2))
		barplot6=barplot(allmodkmecors,col=as.character(allmods),names=as.character(allmods),las=3,main=paste("kME pres. "," ",combmatrix[i,1]," & ",combmatrix[i,2],sep=""),ylab="Correlation",ylim=c(0,1),cex.lab=1.4,cex.main=1.8)
		abline(h=-log(.05,10),col="blue",lwd=2)
		abline(h=-log(.05/length(allmods),10),col="red",lwd=2)
		if(removeGrey==TRUE){
			nongreymods=modules1[,2]!="grey"
			modules1nogrey=as.character(modules1[nongreymods,2])
			SSBTnogrey=SSBT[nongreymods,]
			par(mar=c(7,5,4,2))
			boxplot(-log(SSBTnogrey[,skip1+4],10)~modules1nogrey,col=as.character(allmods),las=3,notch=T,ylab="-Log10 q-value",main=paste("DE (",combmatrix[i,1]," vs ",combmatrix[i,2],")",sep=""),cex.axis=1,cex.lab=1.4,cex.main=1.8)
			abline(h=-log(.05,10),col="red",lwd=2)
			par(mar=c(7,5,4,2))
			boxplot(-log(SSBTnogrey[,skip1+11],10)~modules1nogrey,col=as.character(allmods),las=3,notch=T,ylab="-Log10 q-value",main=paste("DC (",combmatrix[i,1]," vs ",combmatrix[i,2],")",sep=""),cex.axis=1,cex.lab=1.4,cex.main=1.8)
			abline(h=-log(.05,10),col="red",lwd=2)
			} else {
			par(mar=c(7,5,4,2))
			boxplot(-log(SSBT[,skip1+4],10)~modules1[,2],col=as.character(allmods),las=3,notch=T,ylab="-Log10 q-value",main=paste("DE (",combmatrix[i,1]," vs ",combmatrix[i,2],")",sep=""),cex.axis=1,cex.lab=1.4,cex.main=1.8)
			abline(h=-log(.05,10),col="red",lwd=2)
			par(mar=c(7,5,4,2))
			boxplot(-log(SSBT[,skip1+11],10)~modules1[,2],col=as.character(allmods),las=3,notch=T,ylab="-Log10 q-value",main=paste("DC (",combmatrix[i,1]," vs ",combmatrix[i,2],")",sep=""),cex.axis=1,cex.lab=1.4,cex.main=1.8)
			abline(h=-log(.05,10),col="red",lwd=2)
		}
		dev.off()
		
		## Export DE, DC, and DEvsDC for each module:

		dir.create("DE")
		setwd(DEandDCgpsd)
		dir.create("DC")
		setwd(DEandDCgpsd)
		dir.create("DE_vs_DC")
		setwd(DEandDCgpsd)
		
		for(j in c(1:length(allmods))){
			whichmod=allmods[j]
			restmodrow=is.element(modules1[,2],whichmod)
			setwd(paste(DEandDCgpsd,"/","DE",sep=""))
			signifDE=rep("black",length(featurelabels1[restmodrow]))
			signifDE[SSBT[restmodrow,skip1+4]<.05]="red"
			pdf(file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",allmods[j],"_",combmatrix[i,1],"_vs_",combmatrix[i,2],"_DE.pdf",sep=""))
			par(mar=c(5,5,4,2))
			plot(SSBT[restmodrow,skip1+6],SSBT[restmodrow,skip1+7],col=as.character(signifDE),main=paste(combmatrix[i,1]," vs ",combmatrix[i,2],sep=""),type="n",cex.main=1.8,cex.lab=1.4,cex.axis=1.25,xlab=paste("Log2 expression",combmatrix[i,1]),ylab=paste("Log2 expression",combmatrix[i,2]),xlim=c(floor(min(SSBT[restmodrow,skip1+6])),ceiling(max(SSBT[restmodrow,skip1+6]))))
			text(SSBT[restmodrow,skip1+6],SSBT[restmodrow,skip1+7],labels=as.character(featurelabels1[restmodrow]),cex=cexlabels1,col=as.character(signifDE))
			abline(lm(SSBT[restmodrow,skip1+7]~SSBT[restmodrow,skip1+6]),col="red",lwd=2)
			mtext(paste(as.character(allmods[j]),", r = ",signif(cor.test(SSBT[restmodrow,skip1+6],SSBT[restmodrow,skip1+7])$estimate,3),sep=""),line=0.2,cex=0.8)
			dev.off()
			setwd(paste(DEandDCgpsd,"/","DC",sep=""))
			signifDC=rep("black",length(featurelabels1[restmodrow]))
			signifDC[SSBT[restmodrow,skip1+12]<.05]="red"
			pdf(file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",allmods[j],"_",combmatrix[i,1],"_vs_",combmatrix[i,2],"_DC.pdf",sep=""))
			par(mar=c(5,5,4,2))
			plot(SSBT[restmodrow,skip1+9],SSBT[restmodrow,skip1+10],col=as.character(signifDC),main=paste(combmatrix[i,1]," vs ",combmatrix[i,2],sep=""),type="n",cex.main=1.8,cex.lab=1.4,cex.axis=1.25,xlab=paste("kME",combmatrix[i,1]),ylab=paste("kME",combmatrix[i,1],"+",combmatrix[i,2]),xlim=c(floor(min(SSBT[restmodrow,skip1+9])),ceiling(max(SSBT[restmodrow,skip1+9]))))
			text(SSBT[restmodrow,skip1+9],SSBT[restmodrow,skip1+10],labels=as.character(featurelabels1[restmodrow]),cex=cexlabels1,col=as.character(signifDC))
			abline(lm(SSBT[restmodrow,skip1+10]~SSBT[restmodrow,skip1+9]),col="red",lwd=2)
			mtext(paste(as.character(allmods[j]),", r = ",signif(cor.test(SSBT[restmodrow,skip1+9],SSBT[restmodrow,skip1+10])$estimate,3),sep=""),line=0.2,cex=0.8)
			dev.off()
			setwd(paste(DEandDCgpsd,"/","DE_vs_DC",sep=""))
			pdf(file=paste(projectname1,"_",groups1[r],"_x_",colnames(sampleinfo1)[subgroup1],"_",allmods[j],"_",combmatrix[i,1],"_vs_",combmatrix[i,2],"_DE_vs_DC.pdf",sep=""))
			par(mar=c(5,5,4,2))
			plot(-log(SSBT[restmodrow,skip1+4],10),-log(SSBT[restmodrow,skip1+12],10),main=paste(combmatrix[i,1]," vs ",combmatrix[i,2],sep=""),cex.main=1.8,cex.lab=1.4,cex.axis=1.25,xlab="Differential expression (-Log10 q-value)",ylab="Differential connectivity (-Log10 q-value)",type="n",xlim=c(floor(min(-log(SSBT[restmodrow,skip1+4],10))),ceiling(max(-log(SSBT[restmodrow,skip1+4],10)))))
			text(-log(SSBT[restmodrow,skip1+4],10),-log(SSBT[restmodrow,skip1+12],10),labels=as.character(featurelabels1[restmodrow]),cex=cexlabels1)
			mtext(as.character(allmods[j]),line=0.2,cex=0.8)
			abline(v=-log(.05,10),col="red",lwd=2)
			abline(h=-log(.05,10),col="red",lwd=2)
			dev.off()
			}
			
		setwd(DEandDCwd)
		} ## end of for(i in c(1:dim(combvecmatrix)[1])){
	
		setwd(bymodwd)
	
	if(length(groups1)>1){
      print(paste(groups1[r]," processing complete. Hit enter to proceed to next group.",sep=""))
      answer4=readline(prompt="Press enter: ")
      }	
	  
	  setwd(SNrootDir)
	
	} ## end of for(r in c(1:length(indices1)[1])){
	
  } ## end of function
	
	
	

		
# The function ModulePrinComps1 (by Steve Horvath) finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# This requires the R library impute
# The output is a list with 2 components: 
# 1) a data frame of module eigengenes (PCs), 
# 2) a data frame that lists the percent variance explained by the first 5 PCs of a module
# Options: if removeGrey=T, then no output is generated for the grey module.
# Recall that grey often denotes genes outside proper modules. 
		if (exists("ModulePrinComps1" ) ) rm(ModulePrinComps1);
		ModulePrinComps1=function(datexpr,couleur,removeGrey=F, FiveComponents=F) {
			modlevels=levels(factor(couleur))
			if ( removeGrey ) modlevels=setdiff(modlevels, c("grey") ); 
			if (FiveComponents ) {print("To speed up the calculation, we only compute the five principal components of each module.  Therefore, the estimate of the proportion of variance explained is no longer accurate. If you want an accurate estimate of the proportion of var explained, please choose  the option FiveComponents=F ")  ;} 
			PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
			varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
			names(PrinComps)=paste("PC",modlevels,sep="")
			for(i in c(1:length(modlevels)) ){
				print(i)   
				modulename    = modlevels[i]
				restrict1= as.character(couleur)== modulename
# in the following, rows are genes and columns are samples
				datModule=t(datexpr[, restrict1])
				is.saved = FALSE;
				if (exists(".Random.seed"))
				{
					saved.seed = .Random.seed;
					is.saved = TRUE;
				}
				datModule=impute.knn(as.matrix(datModule))
				if ( (length(datModule)==3) && 
					(!is.null(names(datModule))) &&
					(names(datModule)[1]=="data") )
				datModule = datModule$data;
				if (is.saved) .Random.seed = saved.seed;
				datModule=t(scale(t(datModule)))
				if (FiveComponents ) { svd1 =svd(datModule, nu = 5, nv = 5) } else {svd1=svd(datModule)}
				mtitle=paste("PCs of ", modulename," module", sep="")
				varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
# this is the first principal component
				pc1=svd1$v[,1]
				signh1=sign(sum(cor(pc1,  t(datModule))))
				if (signh1 != 0)  pc1=signh1* pc1
				PrinComps[,i]= pc1
			}
			list(PrinComps=PrinComps, varexplained=varexplained)
		}
		
		
		
## plot.mat function for making heat maps (by Sandrine Dudoit, modified slightly)
		plot.mat=function (x, nrgcols = 50, rlabels = FALSE, clabels = FALSE, 
						   rcols = 1, ccols = 1, title = "", ...)
		{
			n <- nrow(x)
			p <- ncol(x)
			image(1:p, 1:n, t(x[n:1, ]), col = rgcolors.func(nrgcols), 
				  axes = FALSE, xlab = "", ylab = "", ...)
			if (length(ccols) == 1) {
				axis(1, at = 1:p, labels = clabels, las = 2, cex.axis = 0.6, 
					 col.axis = ccols)
			}
			if (length(ccols) == p) {
				cols <- unique(ccols)
				for (i in 1:length(cols)) {
					which <- (1:p)[ccols == cols[i]]
					axis(1, at = which, labels = clabels[which], las = 2, 
						 cex.axis = 0.6, col.axis = cols[i])
				}
			}
			if (length(rcols) == 1) {
				axis(2, at = n:1, labels = rlabels, las = 2, cex.axis = 0.6, 
					 col.axis = rcols)
			}
			if (length(rcols) == n) {
				cols <- unique(rcols)
				for (i in 1:length(cols)) {
					which <- (1:n)[rcols == cols[i]]
					axis(2, at = (n:1)[which], labels = rlabels[which], 
						 las = 2, cex.axis = 0.6, col.axis = cols[i])
				}
			}
			mtext(title, side = 3, line = 3)
			box()
		}
		
		
				
## The function "MedianNegs" (by Mike Oldham) will replace negative values within an expression matrix (e.g. those introduced by ComBat) with the median for the corresponding probe set; required prior to log transformation
## for standardScreeningBinaryTrait.
		
		if(exists("MedianNegs")) rm(MedianNegs);
		MedianNegs=function(datE){
			dat2E=as.matrix(datE)
			negatives=dat2E<0
			for (i in c(1:length(dat2E[,1]))){
				if (length(dat2E[negatives[i,][negatives[i,]]])>0){
					negatives1=dat2E[i,]<0
					dat2E[i,][negatives1]=median(dat2E[i,])
				}
			}
			dat2E
		}
		
