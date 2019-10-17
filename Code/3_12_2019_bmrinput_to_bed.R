# read in dmrcate input
library(data.table)
f_pfhxs_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_f_pfhxs__CpGs_withChem.csv", header = T)
f_pfoa_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_f_pfoa__CpGs_withChem.csv", header = T)
f_pfos_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_f_pfos__CpGs_withChem.csv", header = T)
head(f_pfhxs_DMP)
# get hg19 annotation
paste("all ready have chr in dmrcate input")
# .bed file
paste("chrom, chromStart, chromEnd, unadjusted, p_adjust, CpG")
# function for .bed file
bedfile <- function(DMRcate_input){
  # library
  library(tidyverse)
  # library(gtools) the mixedsort() is slow  
  gender = unlist(strsplit( deparse(substitute(DMRcate_input )), "_"))[1]
  chemname = unlist(strsplit( deparse(substitute(DMRcate_input )), "_"))[2]
  # sort the chromosome 
  # if there is 'chrM', using custom defined factor levels
  # prevent scientific notation so bedtools can be used
  options(scipen=999)
  chrOrder = c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  bedforsave = DMRcate_input %>% dplyr::mutate(end = pos + 51,
                                                           chrom = factor(CHR, 
                                                                    levels=chrOrder)) %>% 
                                      # mutate(chromStart = format(pos, scientific = F),
                                      #        chromEnd = format(end, scientific = F)) %>% 
                                      dplyr::select(chrom, pos, end, raw) %>% arrange(chrom)
  # # as numeric, then the format is useless 
  # bedforsave[,2:4] = mutate_all(bedforsave[,2:4], function(x) as.numeric(as.character(x)) )
  write.table(bedforsave, file= paste("~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/", 
                                      Sys.Date(), "_", gender, "_",chemname, ".bed", sep = ""), 
              # col.names F for bedtools sortBed
              # but we will need the colnames later on for combp
                                     quote=F, sep="\t", row.names=F, col.names=F)
}
####################### sort by bedtools as command line tools ###############
# sort by chromosome may not be enough
# sudo apt-get install bedtools
# sortBed -i 2019-03-14_f_pfhxs.bed  > sorted_f_pfhxs.bed   # the default sort by chromsome and start position

bedfile(f_pfhxs_DMP)
bedfile(f_pfoa_DMP)
bedfile(f_pfos_DMP)
# 110671 for f_pfhxs_DMP, scientific format

# for male 
m_pfhxs_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_m_pfhxs__CpGs_withChem.csv", header = T)
m_pfoa_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_m_pfoa__CpGs_withChem.csv", header = T)
m_pfos_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_m_pfos__CpGs_withChem.csv", header = T)
bedfile(m_pfhxs_DMP)
bedfile(m_pfoa_DMP)
bedfile(m_pfos_DMP)

########### add column names to sorted .bed #####################3

add_columnname <- function(sorted_bed){
  dir = "/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/"
  data = fread(paste(dir, sorted_bed, sep = ""), header = F) %>% as.data.frame()
  colnames(data) = c("chrom",	"pos",	"end",	"rawp")
  write.table(data, file= paste(dir, Sys.Date(), "_", sorted_bed, sep = ""), 
              quote=F, sep="\t", row.names=F, col.names=T)
}
add_columnname("sorted_f_pfhxs.bed")

###### for new chems PFNA, PFDA ##########
# read in dmrcate input
library(data.table)
f_pfoa_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_f_pfoa__CpGs_withChem.csv",
                     header = T)
f_pfos_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_f_pfos__CpGs_withChem.csv",
                    header = T)
f_pfhxs_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_f_pfhxs__CpGs_withChem.csv",
                    header = T)
f_pfdea_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_f_pfdea__CpGs_withChem.csv",
                     header = T)
f_pfna_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_f_pfna__CpGs_withChem.csv",
                     header = T)

m_pfoa_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_m_pfoa__CpGs_withChem.csv",
                    header = T)
m_pfos_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_m_pfos__CpGs_withChem.csv",
                    header = T)
m_pfhxs_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_m_pfhxs__CpGs_withChem.csv",
                     header = T)
m_pfdea_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_m_pfdea__CpGs_withChem.csv",
                     header = T)
m_pfna_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/more_2019-10-17_m_pfna__CpGs_withChem.csv",
                    header = T)
dim(m_pfna_DMP)

# function for .bed file
bedfile <- function(DMRcate_input){
  # library
  library(tidyverse)
  # library(gtools) the mixedsort() is slow  
  gender = unlist(strsplit( deparse(substitute(DMRcate_input )), "_"))[1]
  chemname = unlist(strsplit( deparse(substitute(DMRcate_input )), "_"))[2]
  # sort the chromosome 
  # if there is 'chrM', using custom defined factor levels
  # prevent scientific notation so bedtools can be used
  options(scipen=999)
  chrOrder = c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  bedforsave = DMRcate_input %>% dplyr::mutate(end = pos + 51,
                                               chrom = factor(CHR, 
                                                              levels=chrOrder)) %>% 
    # mutate(chromStart = format(pos, scientific = F),
    #        chromEnd = format(end, scientific = F)) %>% 
    dplyr::select(chrom, pos, end, raw) %>% arrange(chrom)
  # # as numeric, then the format is useless 
  # bedforsave[,2:4] = mutate_all(bedforsave[,2:4], function(x) as.numeric(as.character(x)) )
  write.table(bedforsave, file= paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/bed_for_combp/", 
                                      Sys.Date(), "_", gender, "_",chemname, ".bed", sep = ""), 
              # col.names F for bedtools sortBed
              # but we will need the colnames later on for combp
              quote=F, sep="\t", row.names=F, col.names=F)
}
###### run ####333
bedfile(f_pfhxs_DMP)
bedfile(f_pfoa_DMP)
bedfile(f_pfos_DMP)
bedfile(f_pfdea_DMP)
bedfile(f_pfna_DMP)

bedfile(m_pfdea_DMP)
bedfile(m_pfna_DMP)
bedfile(m_pfhxs_DMP)
bedfile(m_pfoa_DMP)
bedfile(m_pfos_DMP)

####################### sort by bedtools as command line tools ###############
# sort by chromosome may not be enough
# sudo apt-get install bedtools
# sortBed -i 2019-03-14_f_pfhxs.bed  > sorted_f_pfhxs.bed   # the default sort by chromsome and start position

add_columnname <- function(sorted_bed){
  dir = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/bed_for_combp/"
  data = fread(paste(dir, sorted_bed, sep = ""), header = F) %>% as.data.frame()
  colnames(data) = c("chrom",	"start",	"end",	"rawp")
  write.table(data, file= paste(dir, Sys.Date(), "_", sorted_bed, sep = ""), 
              quote=F, sep="\t", row.names=F, col.names=T)
}
add_columnname("sorted_f_pfoa.bed")
add_columnname("sorted_f_pfos.bed")
add_columnname("sorted_f_pfhxs.bed")
add_columnname("sorted_f_pfdea.bed")
add_columnname("sorted_f_pfna.bed")

add_columnname("sorted_m_pfoa.bed")
add_columnname("sorted_m_pfos.bed")
add_columnname("sorted_m_pfhxs.bed")
add_columnname("sorted_m_pfdea.bed")
add_columnname("sorted_m_pfna.bed")

##################### run combp #############################
# source activate python27
# for combp https://github.com/brentp/combined-pvalues
##### combp commands
# finally run settings 
# -c 4 --seed 1e-1 --dist 750 -p f_pfhxs --anno hg19 sorted_f_pfhxs.bed
# 
# date: 2019-03-20 15:15:50.316587
# version: 0.46

max(f_pfhxs_DMP$pos)
