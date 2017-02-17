# TODO: Add comment
# 
# Author: cshughes
###############################################################################
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
#want to filter data that isn't present in at least 25% of samples
cpSet.f = subset(cpSet, rowSums(is.na(cpSet[,30:113]))<63)
#output the data
saveRDS(cpSet.f,'ch_feb2017_OvC_cptac_proteinSet_processed.rds')







