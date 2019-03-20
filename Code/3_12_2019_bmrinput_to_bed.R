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
                                     quote=F, sep="\t", row.names=F, col.names=T)
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

##################### run combp #############################
# source activate python27
# for combp https://github.com/brentp/combined-pvalues
##### combp commands
# -c 4, means pvalue at 4th column
# To calclulate autocorrelation from 1 to 500 bases
# cd ~/Documents/gitlab/DMR_combp/combined-pvalues-master/
# python cpv/acf.py -d 1:100000:50 -c 4 data/PFOA_OS_HXS/sorted_f_pfhxs.bed > data/PFOA_OS_HXS/2019-03-14_f_pfhxs_acf.txt
max(f_pfhxs_DMP$pos)
