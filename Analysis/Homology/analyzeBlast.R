# created by Gabriel on 08/02/17
# last modified by Gabriel on 22/02/17

# Here, I read the blast result and remove columns which are not required for OrthoMCl
# (taking the alignment score as distance matrix)

setwd("/path/to/Analysis/Blast/Blastp/")

blast<-read.table("DbMergedSecretomeQuMergedSecretome.csv", stringsAsFactors=FALSE, sep="\t")
blast<-blast[, c(1,2,12)]
rownames(blast)<-NULL
names(blast)<-NULL
write.table(blast, file="DbMergedSecretomeQuMergedSecretomeShort.csv", sep=" ", row.names=FALSE)
