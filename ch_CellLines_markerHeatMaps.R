# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/markers/Routput")
#make a data holder
ids = list()
#read in the Hughes data
ids[[1]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Hughes_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_1.rds')
ids[[2]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Hughes_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_2.rds')
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
#merge the data sets
hData = merge(ids[[1]],ids[[2]],by='Gene')
allData = merge(hData,ids[[3]],by='Gene')
#z-score the data
allData[,2:13] = scale(allData[,2:13],center=TRUE,scale=TRUE)


#######work with the Hughes data first
#bring in some marker data
genesCC = read.table("/Users/cshughes/Documents/projects/clinicalProteomics/markers/Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("/Users/cshughes/Documents/projects/clinicalProteomics/markers/Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
#add an associated identifier
genesCC$type = 'CCC'
genesS$type = 'HGS'
#bind together
hMark = rbind(genesCC,genesS)
colnames(hMark)[1] = 'Gene'
#merge the data with the Hughes signature
hSet = merge(hMark,allData,by='Gene')
#make the heatmap
pdf('ch_OvC-Cell-Lines_Hughes-Marker_HeatMap.pdf')
#make the plot labels and boundaries
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[2],4),rep(brewer.pal(6,'Accent')[3],4),rep(brewer.pal(6,'Accent')[1],4))
#make the correlation heatmap
heatmap.2(
		as.matrix(hSet[3:14]),
		col= rev(colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1)),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram="both",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = x$Gene,
		las=2,
		ColSideColors=ColSideColors,
		## labels
		main='OvC CellLine Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.75,
		cexCol=1.5
)
dev.off()

#######work with the Coscia data
#bring in some marker data
cMark = read.table("/Users/cshughes/Documents/projects/clinicalProteomics/markers/Coscia_67proteinSignature.txt", header=TRUE, sep='\t')
#add an associated identifier
cMark$type = 'COS'
colnames(cMark)[1] = 'Gene'
#merge the data with the Hughes signature
cSet = merge(cMark,allData,by='Gene')
#make the heatmap
pdf('ch_OvC-Cell-Lines_Coscia-Marker_HeatMap.pdf')
#make the plot labels and boundaries
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[2],4),rep(brewer.pal(6,'Accent')[3],4),rep(brewer.pal(6,'Accent')[1],4))
#make the correlation heatmap
heatmap.2(
		as.matrix(cSet[3:14]),
		col= rev(colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1)),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram="both",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = cSet$Gene,
		las=2,
		ColSideColors=ColSideColors,
		## labels
		main='OvC CellLine Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.75,
		cexCol=1.5
)
dev.off()







