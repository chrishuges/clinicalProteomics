# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/Hughes-Coscia_Expression")
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
#build the heatmap with clustering
ovCor = cor(allData[,c(2:13)], use='pairwise.complete.obs', method='pearson')
ovCor[upper.tri(ovCor)] <- NA
#ovCor[lower.tri(ovCor)] <- ovCorB[lower.tri(ovCorB)]
#make the plot
pdf('ch_OvC-Cell-Lines_CorrelationHeatMap.pdf')
#pdf('ch_test.pdf')
#make the plot labels and boundaries
xLabels<- names(ovCor)
mybreaks = seq(0,0.9,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[2],4),rep(brewer.pal(6,'Accent')[3],4),rep(brewer.pal(6,'Accent')[1],4))
#make the correlation heatmap
heatmap.2(
		ovCor,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(ovCor,2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 1.2,
		colsep = 1:12,
		rowsep = 1:12,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		ColSideColors=ColSideColors,
		## labels
		main='OvC Cell Lines',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.2,
		cexCol=1.2
)
dev.off()


