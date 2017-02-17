# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(gplots)
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/identificationMetrics")
#make a list holder
ids = list()
#read in the Hughes data
ids[[1]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Hughes_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_1.rds')
ids[[2]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Hughes_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_2.rds')
#read in the Coscia data
ids[[3]] = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Coscia_2016/Routput/ch_feb2017_Coscia_cell-line_proteinSet.rds')
names(ids) = c('h1','h2','c1')


###########################################
######look at RLE values first
#need to log transform the Hughes data first
for (i in 1:2){
	pep = as.data.frame(ids[[i]])
	pep[,4:9] = log2(pep[,4:9])
	ids[[i]] = pep
}
#iterate through the data and calculate RLE values
for (i in 1:length(ids)){
	pep = as.data.frame(ids[[i]])
	pep.s = pep %>%
			select(Gene,OVCAR3,OVCAR5,OVSAHO,OVISE) %>%
			filter(!is.na(OVCAR3),!is.na(OVCAR5),!is.na(OVSAHO),!is.na(OVISE))
	pep.s$med = apply(pep.s[,2:5],1,function(x) median(x,na.rm=TRUE))
	pep.s[,2:5] = apply(pep.s[,2:5],2,function(x) x - pep.s$med)
	ids[[i]] = pep.s
}
#iterate over to make a plot
#first make some colors
hCols = brewer.pal(6,'Accent')
#now make the plot
pdf('ch_ch_OvC-Cell-Lines_RLEplots.pdf')
for (i in 1:length(ids)){
	cols = rep(hCols[i],4)
	boxplot(ids[[i]][,2:5],
			col = cols,
			pch = 20,
			cex = 0.75,
			ylim = c(-15,15),
			main = names(ids)[i])
	box(lwd=3)
}
dev.off()


###########################################
######do a heat map cluster
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
pdf('ch_OvC-Cell-Lines_HeatMap.pdf')
#make the plot labels and boundaries
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[2],4),rep(brewer.pal(6,'Accent')[3],4),rep(brewer.pal(6,'Accent')[1],4))
#make the correlation heatmap
heatmap.2(
		as.matrix(allData[2:13]),
		col= rev(colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1)),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram="both",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
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
		cexRow=1.5,
		cexCol=1.5
)
dev.off()








