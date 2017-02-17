# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/Coscia_2016/Routput")
#grab the published data
cos<-read.table("/Users/cshughes/Documents/projects/clinicalProteomics/Coscia_2016/proteinExpression.txt", header=TRUE, sep='\t')


#select the columns with expression data
pro = cos[,c(1,3,5:34,39,41)]
#take the top protein ID and gene ID
pro$Accession = sapply(strsplit(as.character(pro$Majority.protein.IDs), ';'),'[', 1)
pro$Gene = sapply(strsplit(as.character(pro$Gene.names), ';'),'[', 1)
#subset the new columns out
pro = pro[,c(35,36,3:32,33:34)]
#rename some columns
colnames(pro)[33:34] = c('pepNum','MW')
#subset out contaminants
pro = subset(pro, !grepl('CON_',pro$Accession))
#subset proteins that have no expression values
pro = subset(pro, rowSums(is.na(pro[,3:32]))<29)
#remove proteins that are not associated with a gene ID
pro = subset(pro, Gene!='')
#output the data
saveRDS(pro,'ch_feb2017_Coscia_cell-line_proteinSet.rds')


