# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/hughes_tissues_2016/Routput")
#grab the protein files
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
#subset any proteins not present in at least 75% of the samples
hReps.sub = subset(hReps, rowSums(is.na(hReps[,4:21]))<6)
#output the data
saveRDS(hReps.sub,'ch_feb2017_OvC_tissues_proteinSet-processed.rds')















