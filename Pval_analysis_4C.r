# R script to calculate P values using binomial test


##########################################################
## Functions

#This function calculates the p value for each window.
#The probability is the % of RE sites in the chromosome covered with reads and for each window a p value is assigned using the binom.test function.
pValCal<-function (windows_data,windows_data_full)
{
	pVal<-c()
	chromosomes <- unique(windows_data[,1])
	for (i in 1:length(chromosomes)) 
	{
		tmp.data<-windows_data[windows_data[,1]==i,]
		tmp.data_full<-windows_data_full[windows_data_full[,1]==i,]
		prob<-sum(tmp.data[,4])/sum(tmp.data_full[,4])
		tmp<-mapply(binom.test, x=as.vector(tmp.data[,4]), n=as.vector(tmp.data_full[,4]), p=prob, alternative="greater")
		pVal<-c(pVal,as.numeric(tmp[3,]))
	}
	#We save both the p-values and the p-scores which are -log10 of the p-values.
	pVal2<- -log10(pVal);pVal2[pVal2 == Inf]<-max(pVal2[pVal2!=Inf])
	list(pVal,pVal2)
}

#################################################################

## Constants
WIND<- N/2 # Half of the size of the window (N) used for the analysis.
data_full<-read.delim("RE_sites_genome_N_windows_with_RE", header=FALSE) #RE windows with the number of total RE per window.

## Calcualte per sample
data_bait< - read.delim("RE_sites_genome_N_windows_with_Sample_RE", header=FALSE) #RE windows with the number of valid RE from the sample per window.

#Calcualte p-values
pValues<-pValCal(data_bait,data_full)
#Save the windows with their p scores for genome browser and other analyses.
Sample_N_pScore<-cbind(data_bait[,1:3],pValues[[2]])
#Each window is represented by the RE site it is spanning
Sample_N_pScore[,3]<-Sample_N_pScore[,3]-WIND
Sample_N_pScore[,2]<-Sample_N_pScore[,3]-4
#save the file
write.table(Sample_N_pScore,"Sample_N_pScore.bedgraph",quote=F,row.names=F,col.names=F,append=F, sep="\t")


