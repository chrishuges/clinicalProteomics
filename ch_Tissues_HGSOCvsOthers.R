# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(limma)
#
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/Hughes-Coscia_Expression")
#grab the protein files
cpt = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput/ch_feb2017_OvC_cptac_proteinSet_processed.rds')
hug = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput/ch_feb2017_OvC_tissues_proteinSet-processed.rds')

#merge the data
allData = merge(cpt[,c(1,30:113)],hug[,c(1,4:9)],by='Gene')
#log transform the data
allData[,2:91] = log2(allData[,2:91])


#want to make a new data frame for each case where you compare it vs. all of the CPTAC ones
#make a data holder
hgsCompare = list()
#loop over the HGS cases
for (i in 1:6){
	eSet = allData[,c(1,i+85,2:85)]
	#find the cases that have the highest correlation with the input case
	corMap = sort(apply(eSet[,3:86],2,function(x) cor(x,eSet[,2],use='pairwise.complete.obs',method='spearman')),decreasing=TRUE)[1:5]
	#filter out the 5 cases with the closest identity
	pSet = as.data.frame(cbind(eSet[,1:2],eSet[,colnames(eSet) %in% names(corMap)]))
	cSet = eSet[,!(colnames(eSet) %in% names(corMap))]
	cSet = cSet[,c(1,3:81)]
	#combine the two sets
	expSet = merge(pSet,cSet,by='Gene')
	#get the differential expression
	limmaSet = as.matrix(expSet[,c(2:86)])
	design <- cbind(HGSOC=c(rep(1,6),rep(0,79)),CPTAC=c(rep(0,6),rep(1,79)))
	fit <- lmFit(limmaSet,design)
	cont.matrix <- makeContrasts(HGSOCvsCPTAC=HGSOC - CPTAC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	expSet$logFC = fit2$coef
	expSet$pVal = fit2$p.value
	expSet$AdjPVal = p.adjust(fit2$p.value, method="BH")
	#order the data
	expSet = expSet[order(expSet$AdjPVal),]
	#input into the list
	hgsCompare[[i]] = expSet[,c(1,87:89)]
	#names(hgsCompare)[i] = paste('hgs',i,sep='')	
}

#write out the data
for (i in 1:length(hgsCompare)){
	write.table(hgsCompare[[i]],paste('ch_Tissues_HGStumour_',i,'.txt',sep=''),quote=FALSE,sep='\t')
}


########plot the two cases of interest
proh = hgsCompare[[6]]
#make a plot of the data
proSD<-sd(proh$logFC, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))

###make an initial plot with all points
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_HGSvCCC_Volcano.pdf')
xCol = col2rgb(ifelse(proh$Gene %in% markCCC, cols[1], ifelse(proh$Gene %in% markHGS, cols[6],'gray80')))
xCex = ifelse(proh$Gene %in% markHGSCCC, 2, 1)
plot(proh$logFC,
		-log10(proh$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type comparison',
		xlim = c(-5,5),
		ylim = c(0,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(proh),sep=""),cex=1.25)
text(4,2.3,paste('p<0.05'),cex=1.25)


