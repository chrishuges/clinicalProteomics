# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(RColorBrewer)
library(beeswarm)
####
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/hgsoc-ccoc_markers")
#grab the protein files
cpt = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput/ch_feb2017_OvC_cptac_proteinSet_processed.rds')
hug = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput/ch_feb2017_OvC_tissues_proteinSet-processed.rds')


#set the markers
markers = c('NID2','CRYAB','CRABP2','WT1','MSLN')
#get the marker data out of the hughes data
hMark = hug[as.character(hug$Gene) %in% markers,c(1,4:15)]
#get the marker data out of the cptac data
cMark = cpt[as.character(cpt$Gene) %in% markers,c(1,30:113)]
#merge the data sets
hcMark = merge(hMark,cMark,by='Gene')
#reshape for a boxplot...just rearranging to get the plot nicer
hcMarkT = as.data.frame(t(log2(hcMark[2:97])))
colnames(hcMarkT) = hcMark[1:5,1]
hcMarkT = hcMarkT[,c(1,3,5,2,4)]
hcMarkT$set = 'cptac'
hcMarkT[1:6,6] = 'hughesHGS'
hcMarkT[7:12,6] = 'hughesCCC'
hcMarkT = hcMarkT[c(1:6,13:96,7:12),]
#get ready for boxplot			
boxIN = melt(hcMarkT)
###make the plot
#set the colors
cols = brewer.pal(6,'Accent')
colVect = rep(c(cols[1],cols[2],cols[3]),5) 
#plot
pdf('ch_OvC-Tissues_BonafideMarkers-Boxplot.pdf')
boxplot(boxIN$value~boxIN$set+boxIN$variable,
		las=2,
		outline = FALSE,
		boxlwd = 1.5,
		#ylim = c(0.65,0.8),
		boxwex=0.5,
		staplelwd=2,
		whisklwd=2,
		col = colVect,
		at =c(1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19)
)
beeswarm(boxIN$value~boxIN$set+boxIN$variable,
		pch=16,
		col='black',
		corral="omit",
		add=TRUE,
		cex=0.55,
		at =c(1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19)
		)
dev.off()


























