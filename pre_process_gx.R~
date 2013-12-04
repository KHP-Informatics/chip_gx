# Austin Hilliard, White lab UCLA, Sep 2009
#
# group of functions for removing outlier probes and samples in microarray data
#
# removeOutlierProbes removes probes outside of specified stdev range from mean, runs on probesets (rows) or samples (cols) of input data 
#
# removeOutlierProbesIterate runs removeOutlierProbes iteratively until no outliers remain
# 
# removeTooManyNAs looks for probesets (rows) and samples (cols) with more than a specified number of missing values and removes them
# 
# outlierSamples computes inter-sample correlations and performs hierarchical clustering to find sample outliers as taught by Mike Oldham (formerly of the Geschwind lab), 
# written prior to MO giving me beta version of his RemoveOutliers function.
### MO's function RemoveOutliers function does this iteratively, as well as running ComBat, with user interactivity
#
# outlierSamplesIterate runs outlierSamples iteratively until user quits or chooses not to remove any more samples
#
# preProc runs all of the above and saves their output lists, including each iteration of processed data
#
#########################################################################

### libraries
library(lattice)
library(Biobase)
library(affy)
library(limma)
library(vsn)
library(preprocessCore)
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
########################################

############################
## negBeadOutlierRepMean 
#############################
negBeadOutlierRepMean <- function(x) { 
		z_out_samp <- abs(  as.numeric( scale(x) )  ) > 2
		mean_pop <- mean(x[z_out_samp==FALSE])
		sd_pop <- sd(x[z_out_samp==FALSE])
		new_x <- ifelse( abs(  as.numeric( scale(x) )  ) > 2, mean_pop, x ) 
		return(new_x)
}

##############
## quantfun ##
##############
quantfun <- function(x) {
as.integer(cut(x, quantile(x, probs=c(seq(0,1,0.20)  ) ), include.lowest=TRUE)) 
}

##############
## var_gene ##
##############
zero_var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
zero_var <- var_gx==0
return(zero_var)
}

mean_probe <- function(gx_matrix) {
gx <- gx_matrix
mean_gx <- apply(gx,1,mean)
return(mean_gx)
}

sd_probe <- function(gx_matrix) {
gx <- gx_matrix
sd_gx <- apply(gx,1,sd)
return(sd_gx)
}

var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
return(var_gx)
}








########################################

# utility functions

cv = function(x) {
	# compute coefficient of variation
	stdev = sd(as.numeric(x), na.rm=TRUE);
	avg = mean(as.numeric(x), na.rm=TRUE);
	cv = stdev/avg;
	return(cv)
}

# run vsn
runVSN = function(data, plot=T) {
	if (!is.matrix(data)) {data = as.matrix(data)};
	dataSet = new("ExpressionSet", exprs = data);
	dataSetVSN = justvsn(dataSet);
	dataVSN = exprs(dataSetVSN);
	dataVSN = as.data.frame(dataVSN);
	
	if (plot) {par(mfrow=c(1,2)); meanSdPlot(data); meanSdPlot(dataSetVSN)};
	
	return(dataVSN);

	}
##########################################################################

removeOutlierProbes = function(data, deviate=3, rowORcol=1) {
	
	# deviate is number of stdevs away from mean probe must surpass to be outlier.
	# 	could be assessing distribution of single probe across samples,
	#	or distribution of all probes on single sample
	#
	# rowORcol chooses whether computation is over probesets(1) or samples(2)
	# 	if 1, looking for single probe outliers within probeset across samples
	#	if 2, looking for outliers within samples
	
	if(rowORcol==1){label=rownames(data)}; 
	if(rowORcol==2){label=names(data)};
	
	data_dim = dim(data);
	
	# set up data frame to hold mean, sd, cv, and cut-offs for each probeset/sample
	data_stats = data.frame(mean=numeric(data_dim[rowORcol]),
							sd=numeric(data_dim[rowORcol]),
							cv=numeric(data_dim[rowORcol]),
							up_thresh=numeric(data_dim[rowORcol]), 
							low_thresh=numeric(data_dim[rowORcol]), 
							row.names=label 
							);
	
	# compute mean, sd, cv, and cut-offs for each probeset/sample
	data_stats[,1] = apply(data, rowORcol, mean, na.rm=TRUE);
	data_stats[,2] = apply(data, rowORcol, sd, na.rm=TRUE);
	data_stats[,3] = data_stats[,2] / data_stats[,1];
	data_stats[,4] = data_stats[,1] + (deviate * data_stats[,2]);
	data_stats[,5] = data_stats[,1] - (deviate * data_stats[,2]);  
	
	# get indices and values of outliers
	outlier_positions = which(data > data_stats[,4] | data < data_stats[,5], arr.ind=T); 
	outlier_vals = data[outlier_positions];
	outlier_vals = cbind(outlier_positions,outlier_vals);
	
	# get number of outliers and compute percentage of total probes they represent
	total_outliers = dim(outlier_positions)[1];
	total_probes = dim(data)[1] * dim(data)[2];
	percent_outliers = 100 * total_outliers / total_probes;
	
	# make copy of data and insert NAs for outliers
	dataClean = data;
	dataClean[data > data_stats[,4] | data < data_stats[,5]] = NA; 

	out = list(data_stats=data_stats,
				deviate=deviate,
				outlier_positions=outlier_positions,
				outlier_vals=outlier_vals,
				total_probes=total_probes,
				total_outliers=total_outliers,
				percent_outliers=percent_outliers,
				dataClean=dataClean
				);
	return(out);
	 
	}

########################################################
	
removeOutlierProbesIterate = function(data, deviate=3, rowORcol=1) {
	
	if(rowORcol==1){cat('Removing probes >',deviate,'stdevs from probeset mean...\n')}; 
	if(rowORcol==2){cat('Removing probes >',deviate,'stdevs from sample mean...\n')};
	
	single_round_outliers = 1;
	outlier_running_count = 0;
	iteration = 0;
	total_probes = dim(data)[1] * dim(data)[2];
	outlier_positions = list();
	outlier_vals = list(); 
	outliers_by_round = c();
	
	while (single_round_outliers > 0) {
	
		out = removeOutlierProbes(data, deviate, rowORcol);
		single_round_outliers = out$total_outliers;
		outlier_running_count = outlier_running_count + out$total_outliers;
		data = out$dataClean;
		iteration = iteration + 1;
		outlier_positions[[iteration]] = out$outlier_positions;
		outlier_vals[[iteration]] = out$outlier_vals; 
		outliers_by_round[iteration] = single_round_outliers;
		cat('Round',iteration, 'outliers:', single_round_outliers, '\n');
	
		};
	
	names(outlier_positions) = paste('round', c(1:iteration), sep='');
	names(outlier_vals) = paste('round', c(1:iteration), sep='');
	names(outliers_by_round) = paste('round', c(1:iteration), sep='');
	percent_outliers = 100 * outlier_running_count / total_probes;
	dataClean = out$dataClean;
	
	cat('\n');
	cat('Total outliers:', outlier_running_count, '\n');
	cat('Percentage of probes that were outliers:', percent_outliers, '\n');
	
	output = list(dataClean=dataClean, 
					outlier_positions=outlier_positions,
					outlier_vals=outlier_vals,
					outliers_by_round=outliers_by_round,
					total_outliers=outlier_running_count,
					total_probes=total_probes,
					percent_outliers=percent_outliers, 
					sd_cutoff=deviate
					);
	return(output);
	
	}
	
#######################################################

removeTooManyNAs = function (data, probe_thresh=NULL, sample_thresh=NULL) {
	
	if(is.numeric(probe_thresh)){
		probe_thresh=probe_thresh;
		} else {
			probe_thresh = floor(ncol(data)/2);
			}
	cat('\nremoving probes with >', probe_thresh, ' missing measurements\n', sep='');
			
	if(is.numeric(sample_thresh)){
		sample_thresh=sample_thresh;
		} else {
			sample_thresh = floor(nrow(data)/2);
			}
	cat('removing samples with >', sample_thresh, ' missing measurements\n', sep='');
	
	countNAs = apply(is.na(data), 1, sum);
	probes_over_thresh = (1:dim(data)[1])[countNAs > probe_thresh];
	cat('\n');
	cat(length(probes_over_thresh),'probes removed \n');
	
	countNAsSamps = apply(is.na(data), 2, sum);
	samples_over_thresh = (1:dim(data)[2])[countNAsSamps > sample_thresh]; 
	cat(length(samples_over_thresh),'samples removed \n');
	
	dataClean = data;
	if (length(probes_over_thresh) > 0) {dataClean = data[-probes_over_thresh, ] };
	if (length(samples_over_thresh) > 0) {dataClean = data[ , -samples_over_thresh] }; 	
	names(probes_over_thresh) = names(countNAs)[probes_over_thresh];
	names(samples_over_thresh) = names(countNAsSamps)[samples_over_thresh];
	
	output = list(dataClean=dataClean,
					countNAs=countNAs, 
					countNAsSamps=countNAsSamps,
					probes_over_thresh=probes_over_thresh,
					samples_over_thresh=samples_over_thresh,
					probe_thresh=probe_thresh,  
					sample_thresh=sample_thresh 
					);
	
	return(output);
	
	}
	
##########################################################

outlierSamples = function(data, thresh=2, showplots=T) {
	
	if (thresh < 0) {thresh = -thresh};
	
	IAC=cor(data,method="p",use="complete.obs");
	clust=hclust(as.dist(1-IAC),method="average");
	
	meanIAC=mean(IAC[upper.tri(IAC)]);
	meanIACdiag=mean(IAC);
	samplemeanIAC=apply(IAC,2,mean);
	sdCorr=sd(samplemeanIAC);    
	numbersd=(samplemeanIAC-meanIACdiag)/sdCorr;
	
	cat('\n');
	cat('Looking for outlier samples (>',thresh,'stdevs from meanIACdiag)...\n');
	cat('meanIAC =',meanIAC,'\n');
	cat('meanIACdiag =',meanIACdiag,'\n');
	cat('\n') 
	
	over_thresh = numbersd < -thresh | numbersd > thresh; 
	samples_to_remove = numbersd[over_thresh];
	#formatted = data.frame(z.IAC=samples_to_remove);
	
	dataClean = data[, !over_thresh]; 
	
	cat('All samples z.IAC: \n');
	print(numbersd);
	cat('\n\n');
	if (length(samples_to_remove)!=0) {
		cat('Possible outliers: \n')
		#print(formatted);
		print(samples_to_remove);
		cat('\n');
		} else {
			cat('No samples >',thresh,'stdevs from meanIACdiag \n');
			};
		
	output = list(dataClean=dataClean,
					IAC=IAC,
					meanIAC=meanIAC,
					meanIACdiag=meanIACdiag,
					samplemeanIAC=samplemeanIAC,
					numbersd=numbersd,
					clust=clust,
					samples_to_remove=samples_to_remove
					);
					
	if (showplots) {
		par(mfrow=c(1,2)); 
		plot(clust); 
		plot(numbersd, type='n'); 
		text(numbersd, labels=names(numbersd), cex=0.75);
		};
		
	return(output);
	
	}
	
##########################################################

outlierSamplesIterate = function (data, IACthresh=2, showplots=T) {
	
	if (IACthresh < 0) {IACthresh = -IACthresh};
	samples_removed = c();
	temp = data;
	
	while (TRUE) {
		
		out = outlierSamples(temp,as.numeric(IACthresh),showplots);
		to_remove = out$samples_to_remove; 
		if (length(to_remove) < 1) {break}; 
		
		answer_raw = readline(prompt='List samples (0 if none) to remove with single spaces in between (no commas): ');
		if (answer_raw==0) {cat('You didn\'t remove any samples!!! \n'); break};
		answer = strsplit(answer_raw, ' ');
		answer = answer[[1]]; 
		
		temp = temp[, -match(answer, names(temp))];
		samples_removed = c(samples_removed, to_remove[names(to_remove)==answer]); 
		cat('Sample(s)', answer_raw, 'removed \n');
			
		};
	
	cat('\n');
	cat('Any more suspicious samples to remove?');
	choose = menu(c('Yes','No'));
	
	if (choose==1) {
		answer_raw = readline(prompt='List samples to remove with single spaces in between (no commas): ');
		answer = strsplit(answer_raw, ' ');
		answer = answer[[1]]; 
		temp = temp[, -match(answer, names(temp))];
		samples_removed = c(samples_removed, out$numbersd[answer]); 
		cat('Sample(s)', answer_raw, 'removed \n');
		} else {
			cat('Ok... finished\n')
			};

	output = list(dataClean=temp, 
					IAC=out$IAC, 
					meanIACdiag=out$meanIACdiag,
					samplemeanIAC=out$samplemeanIAC,
					numbersd=out$numbersd,
					samples_removed=samples_removed
					);
	
	return(output);
	
	}
	
##########################################################


preProcess = function (datIN,
					removeOutlierProbes=T, deviate=3, rowORcol=1,
					removeTooManyNAs=T, probe_thresh=NULL, sample_thresh=NULL,
					removeOutlierSamples=T, IACthresh=2, showplots=T,
					Qnorm=T,
					vsn=F, vsnPlot=F,
					CVsort=F) {
	
	# check input, if ok, assign input to 'temp'. if not, quit 
	if (is.data.frame(datIN) || is.matrix(datIN)) {
		temp = datIN;
		cat('Input data has',nrow(temp),'rows (genes) and',ncol(temp),'columns (samples) \n\n');
		} else {
			stop('Input data must be in the form of a data frame or matrix!');
			};
	
	# if removeOutlierProbes=T, run on 'temp' and save processed data in 'temp' for further processing
	# if removeOutlierProbes=F, skip. 'temp' continues to hold input data.
	# 	set outlierProbesOUTPUT to NULL
	if (removeOutlierProbes) {
		outlierProbesOUTPUT = removeOutlierProbesIterate(temp, deviate, rowORcol);
		temp = outlierProbesOUTPUT$dataClean;
		cat('Processed data available in output as $data_removedOutlierProbes \n');
		} else {
			cat('Skipping removal of outlier probes, checking for probesets and samples with too much missing data...\n');
			outlierProbesOUTPUT = NULL;
			};
	
	# if removeTooManyNAs=T, run on 'temp' and save processed data in 'temp' for further processing	
	# if removeTooManyNAs=F, skip. 'temp' holds output of removeOutlierProbes or input data.
	# 	set checkMissingDataOUTPUT to NULL
	if (removeTooManyNAs) {
		checkMissingDataOUTPUT = removeTooManyNAs(temp, probe_thresh, sample_thresh);
		temp = checkMissingDataOUTPUT$dataClean;		
		cat('...Processed data ($data_checkedMissingData) has',nrow(temp),'rows and',ncol(temp),'columns \n');
		} else {
			cat('Skipping removal of probesets and samples with too much missing data, checking for outlier samples...\n');
			checkMissingDataOUTPUT = NULL;	
			};
	
	# if removeOutlierSamples=T, run on 'temp' and save processed data in 'temp' for further processing
	# if removeOutlierSamples=F, skip. 'temp' holds output of removeTooManyNAs, removeOutlierProbes, or input 
	# 	set outlierSamplesOUTPUT to NULL
	if (removeOutlierSamples) {
		outlierSamplesOUTPUT = outlierSamplesIterate(temp, IACthresh, showplots);
		temp = outlierSamplesOUTPUT$dataClean;
		cat('...Processed data ($data_removedOutlierSamples) now has',nrow(temp),'rows and',ncol(temp),'columns \n\n');
		} else {
			cat('Skipping removal of outlier samples...\n');
			outlierSamplesOUTPUT = NULL;
			};
	
	# if Qnorm=T, run on 'temp', save output to 'tempQnorm' for function output as 'data_Qnorm'.
	# 	set 'temp' equal to 'tempQnorm' for further processing.
	# 	ask if user wants to re-check for outlier samples but don't run iterative version as it is just a check
	# if Qnorm=F, skip and set data_Qnorm to NULL.
	# 	'temp' holds output of removeOutlierSamples, removeTooManyNAs, removeOutlierProbes, or input  
	if (Qnorm) {
		cat('..........\n');
		cat('Performing quantile normalization...\n')
		data_Qnorm = as.data.frame(normalize.quantiles(as.matrix(temp)));
		names(data_Qnorm) = names(temp); rownames(data_Qnorm) = rownames(temp);
		temp = data_Qnorm;
		cat('Normalized data available in output as $data_Qnorm \n\n');
		cat('Re-check for sample outliers?\n')
		answer = menu(c('Yes','No'));
		if (answer==1) {
			postNormOutlierSamples = outlierSamples(data_Qnorm, IACthresh, showplots);
			cat('\n');
			cat('Don\'t remove suspicious samples after normalizing!!!\n');
			cat('Instead, re-run and remove right before normalizing \n\n');
			} else {
				cat('You may want to re-check for outlier samples \n\n');
				};
		} else {
			cat('Skipping quantile normalization...\n\n');
			data_Qnorm = NULL;
			};
	
	# if vsn=T, run on 'temp', save output to 'tempVSN' for function output as 'data_VSN'.
	# 	set 'temp' equal to 'tempVSN' for further processing
	# if vsn=F, skip and set data_VSN to NULL.
	# 	'temp' holds output of Qnorm, removeOutlierSamples, removeTooManyNAs, removeOutlierProbes, or input
	if (vsn) {
		cat('Performing variance stabilization...\n\n');
		data_VSN = runVSN(temp, vsnPlot);
		temp = data_VSN;
		cat('Variance stabilized data available in output as $data_VSN \n\n');
		} else {
			cat ('Skipping variance stabilization...\n\n');
			data_VSN = NULL;
			};
	
	# if CVsort=T, compute CVs of rows in 'temp' and sort rows of 'temp' by CVs
	# 	save sorted 'temp' as 'tempCVsort' for function output as 'data_Sorted'
	# 	'temp' holds output of vsn, Qnorm, removeOutlierSamples, removeTooManyNAs, removeOutlierProbes, or input
	# if CVsort=F, skip and set 'data_Sorted' to NULL
	if (CVsort) {
		cat('Sorting probes by CV...\n')
		data_CVs = apply(temp, 1, cv);
		data_Sorted = temp[order(data_CVs, decreasing=T), ];
		cat('Sorted data available in output as $data_Sorted, CVs in $data_CVs \n\n');
		} else {
			cat('Skipping CV sort...\n\n');
			data_Sorted = NULL;
			data_CVs = NULL;
			};
	
	cat('Creating list for output...\n');
	output = list(outlierProbesOUTPUT=outlierProbesOUTPUT,
					checkMissingDataOUTPUT=checkMissingDataOUTPUT,
					outlierSamplesOUTPUT=outlierSamplesOUTPUT,
					data_removedOutlierProbes=outlierProbesOUTPUT$dataClean,
					data_checkedMissingData=checkMissingDataOUTPUT$dataClean,
					data_removedOutlierSamples=outlierSamplesOUTPUT$dataClean,
					data_Qnorm=data_Qnorm,
					data_VSN=data_VSN,
					data_CVs=data_CVs,
					data_Sorted=data_Sorted
					);
	
	cat('Write any of the output to .csv file?\n');
	choose = menu(c('Yes','No'));
	if (choose==1) {
		choices = names(output)[4:10];
		cat('Choose one...\n');
		choose = 3 + menu(choices=choices);
		file = readline('Enter a name for the output file, including .csv extension... ');
		write_out = output[[choose]];
		write.csv(write_out, file);
		} else {
			cat('You didn\'t write any of the output to a file\n');
			};
	
	cat('All done!\n');
	return(output); 
	
	}
