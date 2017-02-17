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
	expSet = expSet[order(-expSet$logFC),]
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
#proh = hug[,c(1,9)]
#proh[,2] = log2(proh[,2])
#proh = proh[order(proh$hgs6),]
#make a target gene list
geneSet = c('SMTN','FN1','TAGLN','CTHRC1','MMP2','ITGA5','CRABP2','CRYAB','KRT17','KRT19','CDH1','MSLN','MUC1','CTH','QPCT')
#assign size and colors
cols = brewer.pal(6,'YlOrRd')
proh$colors = 'grey80'
proh$sizes = 2
geneSpots = proh$Gene %in% geneSet
proh[geneSpots,5] = cols[5]
proh[geneSpots,6] = 4
###make an initial plot with all points
pdf('ch_OvC-Tissues_Patient-6_RankedAbundance.pdf')
plot(proh[,2],
		col = proh$colors,
		cex = proh$sizes,
		pch = 20,
		ylim = c(-4,4)
)
text(proh[geneSpots,2], proh[geneSpots,1])
dev.off()



#make ranked heatmap for expression in single patient
proh = hug[,c(1,9)]
proh[,2] = log2(proh[,2])
proh = proh[order(-proh$hgs6),]
proh = proh[!is.na(proh$hgs6),]
geneSpots = proh$Gene %in% geneSet
proh[!geneSpots,1] = ''
rankIN = as.matrix(cbind(proh$hgs6,proh$hgs6))

mybreaks = seq(-2,2,by=0.05)
pdf('ch_OvC-Tissues_Patient-6_RankedAbundanceHeat.pdf')
heatmap.2(rankIN,
		col= rev(colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1)),
		symkey=TRUE,
		Rowv=FALSE,
		Colv=FALSE,
		labRow = proh[,1],
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.25,
		cexCol=0.75
		
)
dev.off()


		
		
		
		