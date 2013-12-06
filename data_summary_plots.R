data_summary_plots <- function(data,plotOut) {
## 
## data_summary_plots(data=my.data.frame, plotOut="data_summary_plots.pdf")
##

library(car)

def.par <- par(no.readonly = TRUE) # save default, for resetting...c(5, 4, 4, 2) + 0.1 ## c(bottom, left, top, right)
data_class <- sapply(data ,class) 
class_list <- unique(data_class)
cat(" The following data classes are observed [",unique(data_class),"]","\r","\n")

pdf(file=plotOut,width=11,height=8)

for( class_type in class_list ) { 

cat(" doing ",class_type,"\r","\n")

	      if(class_type =="character") {
					new_data <- data[,data_class==class_type]
					    for(var in names(new_data) ){ 
					    var_table <- sort(table(new_data[,var]),decreasing=FALSE)
					    par(mar=c(5.1, 8, 4, 2.1))
					    barplot(var_table,las=2,main=paste(var), horiz=TRUE)
						    }
	      } 
	      else if(class_type=="factor") {
					      new_data <- data[,data_class==class_type]
					      for(var in names(new_data) ){ 
					      var_table <- sort(table(new_data[,var]),decreasing=FALSE)
					      par(mar=c(5.1, 8, 4, 2.1))
					      barplot(var_table,las=2,main=paste(var), horiz=TRUE)
	      	            	      }
	      } 
	      else if(class_type=="numeric") {
							  new_data <- data[,data_class==class_type]
							  for(var in names(new_data) ){ 
							  nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))    
							  par(mar=c(4.1, 4.1, 1.1, 2.1))		      
							  X <- as.numeric(new_data[,var])
							  boxplot(X, horizontal=TRUE,  outline=TRUE,main=paste(var))      
							  hist(X,xlab=paste(var),breaks=50,main="",prob=TRUE)
							  lines(density(X),col="blue")
							  meanvar <- signif(mean(X,na.rm=TRUE),3)
							  sdvar <- signif(sd(new_data[,var]),3)
							  normtest_sw <- signif(shapiro.test(X)$p.value,3)
							  mtext(paste("Mean:[",meanvar,"] SD:[",sdvar,"] normP:[",normtest_sw,"]",sep=" "), side = 3, adj=1)
							  par(def.par)
							  qqPlot(X,main=paste(var),pch=20,ylab=paste(var),col="blue")
							  							  }
      	      } 
      	      else if(class_type=="integer") {
							  new_data <- data[,data_class==class_type]
							  new_data <- apply(new_data,2,as.numeric)
							  for(var in names(new_data) ){ 
							  nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))    
							  par(mar=c(4.1, 4.1, 1.1, 2.1))		      
							  X <- as.numeric(new_data[,var])
							  boxplot(X, horizontal=TRUE,  outline=TRUE,main=paste(var))		      
							  hist(X,xlab=paste(var),breaks=50,main="",prob=TRUE)
							  lines(density(X),col="blue")
							  meanvar <- signif(mean(X,na.rm=TRUE),3)
							  sdvar <- signif(sd(new_data[,var]),3)
							  normtest_sw <- signif(shapiro.test(X)$p.value,3)
							  mtext(paste("Mean:[",meanvar,"] SD:[",sdvar,"] normP:[",normtest_sw,"]",sep=" "), side = 3, adj=1)
							  par(def.par)
							  qqPlot(X,main=paste(var),pch=20,ylab=paste(var),col="blue")
							  }
	      }
	      else if(class_type=="logical") {
					      new_data <- data[,data_class==class_type]
					      for(var in names(new_data) ){ 
					      var_table <- sort(table(new_data[,var]),decreasing=FALSE)
					      par(mar=c(5.1, 8, 4, 2.1))
					      barplot(var_table,las=2,main=paste(var), horiz=TRUE)
					      par(def.par)
	      	      }
	      }
	  }
	  dev.off()
  
}