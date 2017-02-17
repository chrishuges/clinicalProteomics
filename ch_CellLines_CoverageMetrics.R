# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/identificationMetrics")
#read in the Hughes data
h1 = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_celllines_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_1.rds')
h2 = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_celllines_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_2.rds')
#read in the Coscia data
c1 = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Coscia_2016/Routput/ch_feb2017_Coscia_cell-line_proteinSet.rds')


####lets focus on cell lines that are in both data sets for this (e.g. OVCAR3,OVCAR5,OVSAHO,OVISE)
h1s = select(h1, Gene,pepNum,OVCAR3,OVCAR5,OVSAHO,OVISE)
h2s = select(h2, Gene,pepNum,OVCAR3,OVCAR5,OVSAHO,OVISE)
c1s = select(c1, Gene,pepNum,OVCAR3,OVCAR5,OVSAHO,OVISE)
#put the data sets into a list
ov = list()
ov[[1]] = h1s
ov[[2]] = h2s
ov[[3]] = c1s
names(ov) = c('h1','h2','c1')

#process the list samples and make a storage object for metrics
idMissing = data.frame()
idNoMissing = data.frame()
#loop over
for (i in 1:length(ov)){
	pro = as.data.frame(ov[[i]])
	#lets keep only proteins that have a value in all 4 cell lines
	pro = subset(pro, rowSums(is.na(pro[,c(3:6)]))<1)
	idNoMissing[i,1] = names(ov)[i]
	idNoMissing[i,2] = length(which(!is.na(pro[,3]))) 
	idNoMissing[i,3] = length(which(!is.na(pro[,4]))) 
	idNoMissing[i,4] = length(which(!is.na(pro[,5]))) 
	idNoMissing[i,5] = length(which(!is.na(pro[,6])))
	idNoMissing[i,6] = sum(pro$pepNum)
	#input the data back into the list, and the metrics
	ov[[i]] = pro
}



###plot the number of quantified proteins
#make some colors
col = brewer.pal(9,'RdBu')
cols1 = c(rep(col[1],4),rep(col[1],4),rep(col[9],4))
cols2 = col2rgb(c(rep(col[1],4),rep(col[1],4),rep(col[9],4)))
#make a single vector with the values
bIN = as.vector(as.matrix(t(idNoMissing[,2:5])))
#make the plot
pdf('ch_OvC-Cell-Lines_QuantifiedProteins.pdf')
barplot(bIN,
		#col = rgb(cols2[1,],cols2[2,],cols2[3,],80,maxColorValue=255),
		col = c(rep('grey40',8),rep('grey80',4)),
		horiz = TRUE,
		xlim = c(0,10000))
abline(v=mean(bIN),col='red',lty=2)
abline(v=7889,col='blue',lty=2)
abline(v=7911,col='blue',lty=2)
abline(v=6203,col='blue',lty=2)
dev.off()





#####determine scaling ID cutoffs
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/identificationMetrics")
#read in the Hughes data
h1 = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_celllines_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_1.rds')
h2 = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_celllines_2016/Routput/ch_feb2017_OvC_cell-line_proteinSet_2.rds')
#read in the Coscia data
c1 = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/Coscia_2016/Routput/ch_feb2017_Coscia_cell-line_proteinSet.rds')
####lets focus on cell lines that are in both data sets for this (e.g. OVCAR3,OVCAR5,OVSAHO,OVISE)
h1s = select(h1, Gene,pepNum,OVCAR3,OVCAR5,OVSAHO,OVISE)
h2s = select(h2, Gene,pepNum,OVCAR3,OVCAR5,OVSAHO,OVISE)
c1s = select(c1, Gene,pepNum,OVCAR3,OVCAR5,OVSAHO,OVISE)
#put the data sets into a list
ov = list()
ov[[1]] = h1s
ov[[2]] = h2s
ov[[3]] = c1s
names(ov) = c('h1','h2','c1')
#iterate over the lists to do the counting
#make a data holder
idOutput = data.frame()
idCounter = data.frame()
#do the loop
for (i in 1:length(ov)){
	x = ov[[i]]
	idCounter[1,1] = nrow(subset(x, rowSums(is.na(x[,3:6]))<4))
	idCounter[2,1] = nrow(subset(x, rowSums(is.na(x[,3:6]))<3))
	idCounter[3,1] = nrow(subset(x, rowSums(is.na(x[,3:6]))<2))
	idCounter[4,1] = nrow(subset(x, rowSums(is.na(x[,3:6]))<1))
	#output
	idOutput = rbind(idOutput,idCounter)
}


#plot the data
pdf('ch_OvC-CellLines_QuantifiedProteinsThresholds.pdf')
plot(idOutput[1:4,1],
		type='o',
		pch = 21,
		lwd=3,
		ylim = c(5000,10000),
		xlim = c(0.5,4.5))
lines(idOutput[5:8,1],col='red')
points(idOutput[5:8,1],col='red',pch = 21)
lines(idOutput[9:12,1],col='blue')
points(idOutput[9:12,1],col='blue',pch = 21)
dev.off()















