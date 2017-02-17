# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/Hughes-Coscia_Expression")
#grab the protein files
cpt = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput/ch_feb2017_OvC_cptac_proteinSet_processed.rds')
hug = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput/ch_feb2017_OvC_tissues_proteinSet-processed.rds')


#want to annotate the TCGA data...bring in the subtype information
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
#now we need to assign the colors based on subtype - Differentiated, Immunoreactive, Mesenchymal, Proliferative, Stromal, Unknown
cols = brewer.pal(6,'Accent')
pnnlSubtypes$color = 'grey60'
pnnlSubtypes[grepl('Differentiated',pnnlSubtypes[,2]),3] = cols[1]
pnnlSubtypes[grepl('Immunoreactive',pnnlSubtypes[,2]),3] = cols[2]
pnnlSubtypes[grepl('Mesenchymal',pnnlSubtypes[,2]),3] = cols[3]
pnnlSubtypes[grepl('Proliferative',pnnlSubtypes[,2]),3] = cols[4]
pnnlSubtypes[grepl('Stromal',pnnlSubtypes[,2]),3] = cols[5]

#now we are ready to make the heatmap
#####first do some basic clustering of the samples
#merge the data
allData = merge(cpt.Ordered[,c(1,30:113)],hug[,c(1,4:15)],by='Gene')
#log transform the data
allData[,2:97] = log2(allData[,2:97])
#z-score the data
allData[,2:97] = scale(allData[,2:97],center=TRUE,scale=TRUE)

#build the heatmap with clustering
pdf('ch_OvC-Tissues_HeatMap.pdf')
#make the plot labels and boundaries
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(pnnlSubtypes$color, rep(cols[6],6))
#make the correlation heatmap
heatmap.2(
		as.matrix(allData[,2:91]),
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
		main='OvC Tissue Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.5,
		cexCol=0.75
)
dev.off()


#####what if we used the subtyped clusters
go1<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/dnaReplication_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
go2<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/cellcellComm_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
go3<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/cytoSig_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
go4<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/platelet_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
go5<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/ECM_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
go6<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/complement_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
go7<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/metabolism_geneSet.txt", header=TRUE, sep='\t',stringsAsFactors=FALSE)
#####first build a dataset that is ordered based on subtype rather than identifier. The data are already sorted based on ID
allData.s = allData[,c(1,86:97,2:85)]
#need to add a better sorting device
pnnlSubtypes$sorter = 80
pnnlSubtypes[grepl('Differentiated',pnnlSubtypes[,2]),4] = 30
pnnlSubtypes[grepl('Immunoreactive',pnnlSubtypes[,2]),4] = 40
pnnlSubtypes[grepl('Mesenchymal',pnnlSubtypes[,2]),4] = 50
pnnlSubtypes[grepl('Proliferative',pnnlSubtypes[,2]),4] = 60
pnnlSubtypes[grepl('Stromal',pnnlSubtypes[,2]),4] = 70
#now add these numbers as a new row
allData.s[nrow(allData.s)+1,] = c(10,rep(20,6), rep(25,6), pnnlSubtypes$sorter)
#reorder the data by subtype
allData.Ordered = allData.s[,order(allData.s[nrow(allData.s),])]
#remove the last row with the annotation
allData.Ordered = allData.Ordered[1:nrow(allData.Ordered)-1,]
#re-sort the pnnlSubtypes
pnnlReOrder = pnnlSubtypes[order(pnnlSubtypes$sorter),]


#select the gene set you want
goSub = go1[go1$absolute.correlation>=0.5,]
heatIN = allData.Ordered[as.character(allData.Ordered$Gene) %in% as.character(goSub$gene.symbol),]
#build the heatmap with clustering
pdf('ch_OvC-Tissues_HeatMapSubset_dnaReplication.pdf')
#make the plot labels and boundaries
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(cols[6],6), pnnlReOrder$color)
#make the correlation heatmap
heatmap.2(
		as.matrix(heatIN[,2:91]),
		col= rev(colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1)),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=FALSE,
		dendrogram="row",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
		las=2,
		ColSideColors=ColSideColors,
		## labels
		main='OvC Tissue Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.5,
		cexCol=0.75
)
dev.off()


#use these for the below plot
#first choose your gene set
#the CPTAC set
geneSet = as.character(go1$gene.symbol)

#the Hughes set
genesCC = read.table("/Users/cshughes/Documents/projects/clinicalProteomics/markers/Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("/Users/cshughes/Documents/projects/clinicalProteomics/markers/Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
#add an associated identifier
genesCC$type = 'CCC'
genesS$type = 'HGS'
#bind together
hMark = rbind(genesCC,genesS)
geneSet = as.character(hMark[,1])

#the Coscia gene set
cMark = read.table("/Users/cshughes/Documents/projects/clinicalProteomics/markers/Coscia_67proteinSignature.txt", header=TRUE, sep='\t')
geneSet = as.character(cMark[,1])

#then apply to the data
heatIN = allData.Ordered[as.character(allData.Ordered$Gene) %in% geneSet,]
x = cor(heatIN[,2:91],use='pairwise.complete.obs',method='spearman')
y = cor(allData.Ordered[,2:97],use='pairwise.complete.obs',method='spearman')
x[upper.tri(x)] <- y[upper.tri(y)]
#make the plot
pdf('ch_OvC-Tissues_HeatMap_Hughesset_Corr.pdf')
#make the plot labels and boundaries
mybreaks = seq(0.4,1,by=0.06)
ColSideColors = c(rep(cols[6],6), pnnlReOrder$color)
#make the correlation heatmap
heatmap.2(
		x,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		labRow = '',
		labCol = '',
		sepwidth = c(0.05, 0.05),
		sepcolor = 'white',
		colsep = c(6,32,43,67,78,88),
		rowsep = c(6,32,43,67,78,88),
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



#then apply to the data
heatIN = allData.Ordered[as.character(allData.Ordered$Gene) %in% geneSet,]
#make the plot
pdf('ch_OvC-Tissues_HeatMap_CosciaSet-Genes_Corr.pdf')
#make the plot labels and boundaries
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(cols[6],6), pnnlReOrder$color)
#make the correlation heatmap
heatmap.2(
		as.matrix(heatIN[,2:91]),
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
		main='OvC Tissue Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.5,
		cexCol=0.75
)
dev.off()



