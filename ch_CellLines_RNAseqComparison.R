# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(mygene)
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/nci60_RNA")
#read in the RNAseq data from http://research-pub.gene.com/KlijnEtAl2014/...'A comprehensive transcriptional portrait of human cancer cell lines'
rna<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/nci60_RNA/140331_RPKMExport.txt", header=TRUE, sep='\t')
#select only the cell lines we are using here
rna.s = rna[, c('GeneID', 'NIH.OVCAR.3', 'OVCAR.5', 'OVSAHO', 'OVISE')]
colnames(rna.s) = c('Gene','OVCAR3','OVCAR5','OVSAHO','OVISE')
#convert the Entrez gene numbers to IDs
gn = queryMany(rna.s$Gene,species='human')
colnames(gn)[6] = 'Gene'
#merge the data to annotation the rnaseq
rna.gn = as.data.frame(merge(rna.s,gn,by='Gene'))
#keep only the data we need
rnaSet = unique(rna.gn[,c(7,2:5)])
colnames(rnaSet)[1] = 'Gene'
#save the RNAseq data
saveRDS(rnaSet,'ch_CancerCellLine_rnaSEQ_RPKM.rds')



######get the protein and rna data
rna = readRDS("/Users/cshughes/Documents/projects/clinicalProteomics/nci60_RNA/ch_CancerCellLine_rnaSEQ_RPKM.rds")
#subset the rna to use genes with an RPKM above 1
rna = subset(rna, rowMeans(rna[,2:5],na.rm=TRUE)>1)
#replaces zeroes with NA
rna[,2:5][rna[,2:5]==0]<-NA
rna[,2:5] = log2(rna[,2:5])
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/Hughes-Coscia_Expression")
#make a data holder
ids = list()
#read in the Hughes data
ids[[1]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_celllines_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_1.rds')
ids[[2]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_celllines_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_2.rds')
#read in the Coscia data
ids[[3]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Coscia_2016/Routput/ch_feb2017_Coscia_cell-line_proteinSet.rds')
names(ids) = c('h1','h2','c1')

#need to log transform the Hughes data first
for (i in 1:2){
	pep = as.data.frame(ids[[i]])
	pep[,4:9] = log2(pep[,4:9])
	ids[[i]] = pep
}
#iterate through the data and filter out missing values
for (i in 1:length(ids)){
	pep = as.data.frame(ids[[i]])
	pep.s = pep %>%
			select(Gene,OVCAR3,OVCAR5,OVSAHO,OVISE) %>%
			filter(!is.na(OVCAR3),!is.na(OVCAR5),!is.na(OVSAHO),!is.na(OVISE))
	ids[[i]] = pep.s
}

######capture the correlation with the RNA data
#data storage object
rpCorMatch = data.frame()
rpCorOther = data.frame()
#loop over the files
for (i in 1:length(ids)){
	rna.pro = merge(rna,ids[[i]],by='Gene')
	message(nrow(rna.pro))
	rpCor = cor(rna.pro[,2:9],use='pairwise.complete.obs',method='spearman')[5:8,1:4]
	#get the matched correlations
	rpCorMatch[i,1] = rpCor[1,1]
	rpCorMatch[i,2] = rpCor[2,2]
	rpCorMatch[i,3] = rpCor[3,3]
	rpCorMatch[i,4] = rpCor[4,4]
	#get the correlations of the others
	rpCorOther[i,1] = mean(rpCor[c(2,3,4),1])
	rpCorOther[i,2] = mean(rpCor[c(1,3,4),2])
	rpCorOther[i,3] = mean(rpCor[c(1,2,4),3])
	rpCorOther[i,4] = mean(rpCor[c(1,2,3),4])
}
colnames(rpCorMatch) = c('OVCAR3','OVCAR5','OVSAHO','OVISE')
#make the colors
cols = c('grey40','grey40','grey80')
#plot the data
pdf('ch_OvC-Cell-Lines_RNACorrelations.pdf')
for (i in 1:nrow(rpCorOther)){
barplot(t(rpCorMatch[i,1:4]),
		#col = rgb(cols2[1,],cols2[2,],cols2[3,],80,maxColorValue=255),
		col = cols[i],
		ylim = c(0,1),
		beside=TRUE,
		space=0.5)
abline(h=rpCorOther[i,1],col='red',lty=2)
abline(h=rpCorOther[i,2],col='red',lty=2)
abline(h=rpCorOther[i,3],col='red',lty=2)
abline(h=rpCorOther[i,4],col='red',lty=2)
}
dev.off()










