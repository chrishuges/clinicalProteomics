# TODO: Add comment
# 
# Author: cshughes
###############################################################################
######read in the data sets
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/identificationMetrics")
#grab the protein files
cpt = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput/ch_feb2017_OvC_cptac_proteinSet_processed.rds')
hug = readRDS('/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput/ch_feb2017_OvC_tissues_proteinSet-processed.rds')

#####first just get the raw number of IDs
#data storage object
hIds = data.frame()
#loop over the Hughes data
for (i in 4:21){
	hIds[i-3,1] = colnames(hug)[i]
	hIds[i-3,2] = length(which(!is.na(hug[,i])))	
}
#data storage object
cIds = data.frame()
#loop over the CPTAC data
for (l in 30:113){
	cIds[l-29,1] = colnames(cpt)[l]
	cIds[l-29,2] = length(which(!is.na(cpt[,l])))	
}
#bind together
ids = rbind(hIds,cIds)
colnames(ids) = c('Sample','Proteins')


###plot the number of quantified proteins
#make some colors
col = brewer.pal(9,'Accent')
cols1 = c(rep(col[1],6),rep(col[2],6),rep(col[3],6),rep(col[4],84))
cols2 = col2rgb(c(rep(col[1],6),rep(col[2],6),rep(col[3],6),rep(col[4],84)))
#make the plot
pdf('ch_OvC-Tissues_QuantifiedProteins.pdf')
barplot(ids$Proteins,
		#col = rgb(cols2[1,],cols2[2,],cols2[3,],80,maxColorValue=255),
		col = c(rep('grey30',6),rep('grey40',6),rep('grey50',6),rep('grey70',84)),
		horiz = TRUE,
		xlim = c(0,8500),
		space=c(rep(0.5,6),1,rep(0.5,5),1,rep(0.5,5),1,rep(0.5,83))
)
dev.off()



#####determine scaling ID cutoffs
#grab the protein files for the cptac data
infiles = dir("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput", pattern=".rds", full.names=TRUE)
#read in the files into a list
cPro <- lapply(infiles, readRDS)
names(cPro) = c(sub(".*?proteinSet_(.*?)(\\.rds*|$)", "\\1", infiles))
####need to bind everything together
cPro.a = Reduce(function(x, y) merge(x, y, by='Gene',all=TRUE), cPro)
#rearrange the data frame to put peptides at the front
cpPepNums = cPro.a[,c(1:2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94,98,102,106,110)]
cpExpNums = cPro.a[,c(1,3:5,7:9,11:13,15:17,19:21,23:25,27:29,31:33,35:37,39:41,43:45,47:49,51:53,55:57,59:61,63:65,67:69,71:73,75:77,79:81,83:85,87:89,91:93,95:97,99:101,103:105,107:109,111:113)]
#make a total set
cpSet = cbind(cpPepNums,cpExpNums[,2:85])

#grab the protein files for the hughes data
infiles = dir("/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput", pattern=".rds", full.names=TRUE)
#read in the files into a list
hPro <- lapply(infiles, readRDS)
names(hPro) = c(sub(".*?tissues_(.*?)(\\.rds*|$)", "\\1", infiles))
#want to collapse the replicates first
#make a data storage object
hRep1 = data.frame()
for (i in 1:3){
	hRep1 = rbind(hRep1,hPro[[i]])	
	
}
#now for the second replicate
hRep2 = data.frame()
for (i in 4:6){
	hRep2 = rbind(hRep2,hPro[[i]])	
	
}
#aggregate each of them
hRep1.a = aggregate(cbind(pepNum,hgs1,hgs2,hgs3,ccc1,ccc2,ccc3,enoc1,enoc2,enoc3)~Gene,data=hRep1,na.action=na.pass,mean,na.rm=TRUE)
hRep2.a = aggregate(cbind(pepNum,hgs1,hgs2,hgs3,ccc1,ccc2,ccc3,enoc1,enoc2,enoc3)~Gene,data=hRep2,na.action=na.pass,mean,na.rm=TRUE)
#merge the two sets
hReps = merge(hRep1.a,hRep2.a,by='Gene',all=TRUE)
hReps = hReps[,c(1,2,12,3:5,13:15,6:8,16:18,9:11,19:21)]
colnames(hReps) = c('Gene','pepNumA','pepNumB','hgs1','hgs2','hgs3','hgs4','hgs5','hgs6','ccc1','ccc2','ccc3','ccc4','ccc5','ccc6','enoc1','enoc2','enoc3','enoc4','enoc5','enoc6')



#now apply different filters
#do the hughes data first
#make a data holder
hIds = data.frame()
hIds[1,1] = nrow(subset(hReps, rowSums(is.na(hReps[,4:9]))<6))
hIds[2,1] = nrow(subset(hReps, rowSums(is.na(hReps[,4:9]))<5))
hIds[3,1] = nrow(subset(hReps, rowSums(is.na(hReps[,4:9]))<4))
hIds[4,1] = nrow(subset(hReps, rowSums(is.na(hReps[,4:9]))<3))
hIds[5,1] = nrow(subset(hReps, rowSums(is.na(hReps[,4:9]))<2))
hIds[6,1] = nrow(subset(hReps, rowSums(is.na(hReps[,4:9]))<1))

#do the cptac data
#make a data holder
cIds = data.frame()
#keep only the cptac expression data
cpIDset = cpSet[,30:113]
#get the initial number of proteins
cIds[1,1] = nrow(cpIDset)
#want to randomly sample a set of 6 patients
#data holder
idSet = data.frame()
for (i in 1:1000){
	sampSet = cpIDset[ ,sample(ncol(cpIDset), 6)]
	idSet[i,1] = nrow(subset(sampSet, rowSums(is.na(sampSet[,1:6]))<1))	
}
cIds[6,1] = mean(idSet[,1])

#assign annotation to the data sets
hIds$Set = 'hughes'
cIds$Set = 'cptac'
colnames(hIds)[1] = 'Proteins'
colnames(cIds)[1] = 'Proteins'
#bind together
idSet = rbind(hIds,cIds)


#plot the data
pdf('ch_OvC-Tissues_QuantifiedProteinsThresholds.pdf')
plot(idSet[1:6,1],
		type='o',
		lty=2,
		lwd=3,
		ylim = c(4000,12000),
		xlim = c(0.5,6.5))
lines(idSet[7:12,1],col='red')
points(idSet[7:12,1],col='red')
dev.off()


