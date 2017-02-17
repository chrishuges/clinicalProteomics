# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#set the directory
setwd(dir="/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian/Routput")
#grab the protein files
infiles = dir("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian", pattern="Proteins\\.txt", full.names=TRUE)
#read in the files into a list
proSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(proSet) = c(sub(".*?PNNL_(.*?)(\\_.*|$)", "\\1", infiles))
#do the same for the PSM files
infiles = dir("/Users/cshughes/Documents/projects/clinicalProteomics/cptac_ovarian", pattern="PSMs\\.txt", full.names=TRUE)
#read in the files into a list
psmSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(psmSet) = c(sub(".*?PNNL_(.*?)(\\_.*|$)", "\\1", infiles))

#117 is the pooled standard
############process the data into a usable matrix
#make a data holder
pepExp = list()
#assign a generic sample name
sampleNames = c('p1','p2','p3','pool')
#loop over the files
for (i in 1:length(proSet)){
	pogTrack = data.frame()
	#get the peptide data set out
	pep1 = psmSet[[i]] %>%
			select(Annotated.Sequence,Number.of.Protein.Groups,Master.Protein.Accessions,X114,X115,X116,X117) %>%
			filter(Number.of.Protein.Groups == 1) %>%
			filter(Master.Protein.Accessions != 'sp')
	#do some basic text processing on the accession and sequence columns
	pep1$Accession = sapply(strsplit(as.character(pep1$Master.Protein.Accessions), ';'),'[', 1)
	pep1$Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep1$Annotated.Sequence))
	#filter out the old data
	pep2 = pep1 %>%
			select(Accession,Sequence,X114,X115,X116,X117)
	#get the protein data set out
	pro1 = proSet[[i]] %>%
			select(Accession,Description,MW.in.kDa,Number.of.AAs)
	#add the Gene column
	pro1$Gene = sub(".*?GN=(.*?)( .*|$)", "\\1", pro1$Description)
	proIndex = select(pro1, Accession, Gene, MW.in.kDa, Number.of.AAs)
	#merge the data
	pep2$numNA = rowSums(is.na(pep2[,3:6]))
	pep2$meanSN = rowMeans(pep2[,3:6],na.rm=TRUE)
	#filter out redundant data
	pep3 = pep2 %>%
			filter(numNA < 6) %>%
			filter(meanSN > 5)
	#aggregate non-unique peptides
	pep4 = aggregate(cbind(X114,X115,X116,X117)~Accession+Sequence,data=pep3,na.action=na.pass,FUN=mean,na.rm=TRUE)
	pep4$pepNum = 1
	#roll up data into protein values
	pro2 = aggregate(cbind(X114,X115,X116,X117,pepNum)~Accession,data=pep4,na.action=na.pass,FUN=sum,na.rm=TRUE)
	message(nrow(pro2),' proteins in file.')
	#lookup the relevant data from the protein table
	pro3 = merge(proIndex,pro2,by='Accession')
	#merge non-unique genes
	pro3 = aggregate(cbind(X114,X115,X116,X117,pepNum,Number.of.AAs)~Gene,data=pro3,na.action=na.pass,FUN=mean,na.rm=TRUE)
	#finish the sinQ calculation
	for (j in 2:5){
		sinTOT = sum(pro3[,j],na.rm=TRUE)
		pro3[,j] = pro3[,j]/sinTOT
		pro3[,j] = (pro3[,j]/pro3$Number.of.AAs)*1e7
	}
	colnames(pro3)[2:5] = sampleNames
	#replaces zeroes with NA
	pro3[,2:5][pro3[,2:5]==0]<-NA
	#replace the column names with unique ones
	pro4 = select(pro3, Gene,pepNum,p1,p2,p3)
	colnames(pro4) = c('Gene',paste('pepNum_',names(psmSet)[i],sep=''),paste('TCGA-',sapply(strsplit(as.character(infiles[i]), '_'),'[',3),sep=''),paste('TCGA-',sapply(strsplit(as.character(infiles[i]), '_'),'[',4),sep=''),paste('TCGA-',sapply(strsplit(as.character(infiles[i]), '_'),'[',5),sep=''))
	#input back into the list	
	pepExp[[i]] = pro4
	#counter
	message('Finished ',i,' files.')
}
#save the processed data objects
for (i in 1:length(pepExp)){
	saveRDS(as.data.frame(pepExp[[i]]),paste('ch_feb2017_OvC_cptac_proteinSet_',names(psmSet)[i],'.rds',sep=''))	
}




