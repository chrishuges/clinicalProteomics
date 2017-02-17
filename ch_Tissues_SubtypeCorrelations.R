# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(reshape2)
library(beeswarm)
library(TCGA2STAT)
####
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/Hughes-Coscia_Expression")
#grab the protein files
cpt = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput/ch_feb2017_OvC_cptac_proteinSet_processed.rds')
hug = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput/ch_feb2017_OvC_tissues_proteinSet-processed.rds')
#want to annotate the TCGA data...so we can order by subtype
cpSubtypes<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/proteomicSubtypes.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
#determine which TCGA cases we have
cpIDs = data.frame()
cpIDs[1:84,1] = colnames(cpt[,30:113])
colnames(cpIDs)[1] = 'Tumor'
pnnlSubtypes = merge(cpIDs,cpSubtypes,by='Tumor',all=TRUE)
#assign a couple of cases that don't seem to be in the annotation
pnnlSubtypes[is.na(pnnlSubtypes[,2]),2] = 'Unknown'
pnnlSubtypes = pnnlSubtypes[which(as.character(pnnlSubtypes$Tumor) %in% colnames(cpt)),]
#reorder the cptac data by column name
cpt.Ordered = cpt[,order(names(cpt))]
#merge the expression data
allData = merge(cpt.Ordered[,c(1,30:113)],hug[,c(1,4:9)],by='Gene')
#log transform the data
allData[,2:91] = log2(allData[,2:91])
#z-score the data
allData[,2:91] = scale(allData[,2:91],center=TRUE,scale=TRUE)



############################
#now lets make a loop that calculates the correlation by subtype
subtypes = c('Differentiated','Immunoreactive','Mesenchymal','Proliferative','Stromal')

#data holder
subCors = data.frame()
caseCors = data.frame()
#loop
for (i in 1:length(subtypes)){
	tumorSet = which(pnnlSubtypes$Proteomic.subtype == subtypes[i])
	#now loop through the cases and get the correlation values
	for (l in 1:6){
		eData = allData[,c(l+85,tumorSet+1)]
		pCor = apply(eData[,2:ncol(eData)],2,function(x) cor(x,eData[,1],use='pairwise.complete.obs',method='spearman'))
		subCors[1:length(pCor),l] = pCor
	}
	subCors$subtype = subtypes[i]
	colnames(subCors) = c('hgs1','hgs2','hgs3','hgs4','hgs5','hgs6','subtype')	
	caseCors = rbind(caseCors,subCors)
}
write.table(caseCors,'ch_SubtypesCorrelation.txt',quote=FALSE,sep='\t')
#melt the data to long format
boxIN = melt(caseCors)
###make the plot
#set the colors
cols = brewer.pal(6,'Accent')
#plot
pdf('ch_OvC-Tissues_SubtypesCorrelation-Boxplot.pdf')
boxplot(boxIN$value~boxIN$subtype+boxIN$variable,
		las=2,
		outline = FALSE,
		boxlwd = 1.5,
		ylim = c(0.65,0.8),
		boxwex=0.5,
		staplelwd=2,
		whisklwd=2,
		col = cols[1:5],
		at =c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35)
)
beeswarm(boxIN$value~boxIN$subtype+boxIN$variable,
		pch=16,
		col='black',
		corral="omit",
		add=TRUE,
		cex=0.35,
		at =c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35))
dev.off()



##########compare the subtypes with the RNAseq data
rnaseq.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM")
rnaEXP = as.data.frame(rnaseq.ov[[1]])
colnames(rnaEXP) = paste(sapply(strsplit(colnames(rnaEXP),'-'),'[',1),'-',sapply(strsplit(colnames(rnaEXP),'-'),'[',2),'-',sapply(strsplit(colnames(rnaEXP),'-'),'[',3),sep='')
#keep only those that are in our annotation table
rnaSUB = rnaEXP[,as.character(names(rnaEXP)) %in% as.character(pnnlSubtypes$Tumor)]
#add a column for gene names to be used in merging later
rnaSUB$Gene = row.names(rnaSUB)
row.names(rnaSUB) = NULL
#subset rows with a mean RPKM <1
rnaSUB = subset(rnaSUB, rowMeans(rnaSUB[,1:57])>1)
#log transform
rnaSUB[,1:57] = log2(rnaSUB[,1:57])
#sort by column names
rnaSUB.Ordered = rnaSUB[,order(names(rnaSUB))]
#keep only the pnnlSubtypes that we have
tcgaSubtypes = pnnlSubtypes[as.character(pnnlSubtypes$Tumor) %in% names(rnaSUB.Ordered),]

#merge with the protein data
hugRNA = merge(rnaSUB.Ordered,hug[,c(1,4:9)],by='Gene')
hugRNA[,2:64] = scale(hugRNA[,2:64],scale=TRUE,center=TRUE)
############################
#now lets make a loop that calculates the correlation by subtype
subtypes = c('Differentiated','Immunoreactive','Mesenchymal','Proliferative','Stromal')

#data holder
subCors = data.frame()
caseCors = data.frame()
#loop
for (i in 1:length(subtypes)){
	tumorSet = which(tcgaSubtypes$Proteomic.subtype == subtypes[i])
	#now loop through the cases and get the correlation values
	for (l in 1:6){
		eData = allData[,c(l+58,tumorSet+1)]
		pCor = apply(eData[,2:ncol(eData)],2,function(x) cor(x,eData[,1],use='pairwise.complete.obs',method='spearman'))
		subCors[1:length(pCor),l] = pCor
	}
	subCors$subtype = subtypes[i]
	colnames(subCors) = c('hgs1','hgs2','hgs3','hgs4','hgs5','hgs6','subtype')	
	caseCors = rbind(caseCors,subCors)
}
#write.table(caseCors,'ch_SubtypesCorrelation_RNA.txt',quote=FALSE,sep='\t')
#melt the data to long format
boxIN = melt(caseCors)
###make the plot
#set the colors
cols = brewer.pal(6,'Accent')
#plot
pdf('ch_OvC-Tissues_SubtypesCorrelation_RNA_Boxplot.pdf')
boxplot(boxIN$value~boxIN$subtype+boxIN$variable,
		las=2,
		outline = FALSE,
		boxlwd = 1.5,
		ylim = c(0.65,0.8),
		boxwex=0.5,
		staplelwd=2,
		whisklwd=2,
		col = cols[1:5],
		at =c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35)
)
beeswarm(boxIN$value~boxIN$subtype+boxIN$variable,
		pch=16,
		col='black',
		corral="omit",
		add=TRUE,
		cex=0.35,
		at =c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35))
dev.off()


