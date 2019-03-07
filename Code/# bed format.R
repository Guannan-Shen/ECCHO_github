# original .bed 
lnpfoa_female_mval <- read.table(file= "~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/lnpfoa_female_mval.bed",
                                 sep = "\t")
lnpfoa_female_mval$V1 <- sort(as.character(lnpfoa_female_mval$V1))
write.table(lnpfoa_female_mval, file= "~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/lnpfoa_female_mval.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

