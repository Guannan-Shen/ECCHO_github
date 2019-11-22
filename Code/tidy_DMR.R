## set up workspace
rm(list = ls())
library(data.table)
library(tidyverse)
library(magrittr)
library(stats)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tibble)
library(grid) # low-level grid functions are required
library(openxlsx)
options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
## not in function
'%nin%' <- Negate('%in%')


######### read and sort ###########
dir2 <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/"
files <- c("f_pfdea.anno.hg19.bed", "f_pfna.anno.hg19.bed","f_pfoa.anno.hg19.bed", "f_pfos.anno.hg19.bed","f_pfhxs.anno.hg19.bed",
           "m_pfdea.anno.hg19.bed", "m_pfna.anno.hg19.bed","m_pfoa.anno.hg19.bed", "m_pfos.anno.hg19.bed","m_pfhxs.anno.hg19.bed")

dmr_read_sort <- function(file, dir, z_pcut){
  data = read.delim(paste0(dir, file), header=TRUE, sep="\t")
  dmr_sort = data.frame(data) %>% arrange(z_sidak_p) %>% 
    dplyr::mutate(No. = 1:nrow(data) ) %>% 
    dplyr::filter(z_sidak_p <= z_pcut) %>% 
    dplyr::rename(chrom = X.chrom) %>%
    select(-c("min_p", "z_p") ) %>%
    select(No. , everything())
  write.csv(dmr_sort, file =  paste0(dir, "tidy_dmr/", z_pcut, file, ".csv"),
            row.names = F)
  return(dmr_sort)
}

for ( i in 1:length(files)){
  dmr_read_sort(files[i], dir2, 0.05)
}


############### read DMR to get DMPs (CpGs per DMR) ###################
#########  get DMP
dir <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/"
f_pfdea_DMP <- fread(paste0(dir,"more_2019-10-17_f_pfdea__CpGs_withChem.csv"), header = T)
f_pfna_DMP <- fread(paste0(dir,"more_2019-10-17_f_pfna__CpGs_withChem.csv"), header = T)
f_pfoa_DMP <- fread(paste0(dir,"more_2019-10-17_f_pfoa__CpGs_withChem.csv"), header = T)
f_pfos_DMP <- fread(paste0(dir,"more_2019-10-17_f_pfos__CpGs_withChem.csv"), header = T)
f_pfhxs_DMP <- fread(paste0(dir,"more_2019-10-17_f_pfhxs__CpGs_withChem.csv"), header = T)

m_pfdea_DMP <- fread(paste0(dir,"2019-09-30_m_pfdea__CpGs_withChem.csv"), header = T)
m_pfna_DMP <- fread(paste0(dir,"2019-09-30_m_pfna__CpGs_withChem.csv"), header = T)
m_pfoa_DMP <- fread(paste0(dir,"more_2019-10-17_m_pfoa__CpGs_withChem.csv"), header = T)
m_pfos_DMP <- fread(paste0(dir,"more_2019-10-17_m_pfos__CpGs_withChem.csv"), header = T)
m_pfhxs_DMP <- fread(paste0(dir,"more_2019-10-17_m_pfhxs__CpGs_withChem.csv"), header = T)

########### get DMR
dir2 <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/"

f_pfdea_dmr <- read.delim(paste0(dir2, "f_pfdea.anno.hg19.bed"), 
                         header=TRUE, sep="\t")
f_pfna_dmr <- read.delim(paste0(dir2, "f_pfna.anno.hg19.bed"), 
                         header=TRUE, sep="\t")
f_pfoa_dmr <- read.delim(paste0(dir2, "f_pfoa.anno.hg19.bed"), 
                          header=TRUE, sep="\t")
f_pfos_dmr <- read.delim(paste0(dir2, "f_pfos.anno.hg19.bed"), 
                         header=TRUE, sep="\t")
f_pfhxs_dmr <- read.delim(paste0(dir2, "f_pfhxs.anno.hg19.bed"), 
                         header=TRUE, sep="\t")


m_pfdea_dmr <- read.delim(paste0(dir2, "m_pfdea.anno.hg19.bed"), 
                          header=TRUE, sep="\t")
m_pfna_dmr <- read.delim(paste0(dir2, "m_pfna.anno.hg19.bed"), 
                         header=TRUE, sep="\t")
m_pfoa_dmr <- read.delim(paste0(dir2, "m_pfoa.anno.hg19.bed"), 
                         header=TRUE, sep="\t")
m_pfos_dmr <- read.delim(paste0(dir2, "m_pfos.anno.hg19.bed"), 
                         header=TRUE, sep="\t")
m_pfhxs_dmr <- read.delim(paste0(dir2, "m_pfhxs.anno.hg19.bed"), 
                          header=TRUE, sep="\t")

anno = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

###################################################################
####### get individual probes and mean and max beta fold change ###
###################################################################

getprobes_perDMR = function(dmr, dmp, anno, pcutoff){
  # name contains pcutoff
  name = paste(substitute(dmr), pcutoff, sep = "_")
  dmr_sig = data.frame(dmr) %>% arrange(z_sidak_p) %>% 
    dplyr::mutate(No. = 1:nrow(dmr) )  %>% 
    dplyr::rename(chrom = X.chrom) %>%
    select(-c("min_p", "z_p") ) %>%
    select(No. , everything()) %>% filter(z_sidak_p < pcutoff)
  # chrom, start, end
  chrom = dmr_sig[, 2]; start =  dmr_sig[, 3]; stop =  dmr_sig[, 4]; n_dmr = nrow(dmr_sig)
  max_probes = max(dmr_sig$n_probes)
  
  # 2 means probename and rawp, + 3 means cpg1 and positive rate No_probes
  # 3 means probename and rawp and beta 
  maxlength = (3*max_probes + 3)
  dmrprobes_all = matrix(NA, 0, maxlength)
  n_probes = NULL
  
  for (i in 1:n_dmr) {
    # this is to get cpgs names by their position 
    anno_probes = data.frame(anno) %>% dplyr::filter(chr == chrom[i]) %>%
      dplyr::filter(pos >= start[i]) %>% dplyr::filter (pos <= stop[i]) %>% dplyr::select(Name)
    
    probes_dmp = data.frame(dmp) %>% 
      filter(ID %in% anno_probes$Name)  %>% 
      select(ID, betafc, raw ) %>% dplyr::rename(beta = betafc, raw_p = raw)
    # calculate the maximum and mean of betafold change of all cpgs across all cpgs within a DMR
    # which.max and return the first for ties
    max_beta = probes_dmp$beta[which.max(abs(probes_dmp$beta) )][1]
    mean_beta = base::mean(probes_dmp$beta)
    
    probes = merge(anno_probes, probes_dmp, by.x = "Name", by.y = "ID")
    probes_order = data.frame(probes) %>% dplyr::arrange(raw_p )
    n_probes = c( n_probes, nrow(probes_order) )
    
    if (n_probes[i] == dmr_sig$n_probes[i] ) {
      dmrprobes = matrix(NA, 1, 0)
      for (j in 1: n_probes[i])
      {
        ## beta at column 2 
        # keep beta not probes_order[j, -2]
        dmrprobes = cbind(dmrprobes, probes_order[j, ])
        colnames(dmrprobes)[ (3*j - 2): (3*j )] = c(paste("Name", j, sep = "") ,
                                                    ##  beta means hyper or hypo direction
                                                    paste("direction", j, sep = ""),
                                                    paste("raw_p", j, sep = "") )
      }
      # build a single row df
      dmrprobes = dmrprobes %>% dplyr::mutate(cpg1 = probes_order$Name[1],
                                              positive =  sum(probes_order$beta > 0)/n_probes[i],
                                              No_probes = n_probes[i] ,
                                              maxbeta = max_beta,
                                              meanbeta = mean_beta) %>%
        dplyr::select(cpg1, direction1, 
                      maxbeta, meanbeta, positive, No_probes, everything())
      # Convert a row of a data frame to vector
      dmrprobes = unlist(dmrprobes)
      # cbind or rbind different lengths vectors without repeating the elements of the shorter vectors
      length(dmrprobes) = maxlength
    }
    else {
      # dmrprobes = c(probes_order$Name[1], 0  ,"Wrong Number of probes")
      # names(dmrprobes) <- c("cpg1", "positive", "No_probes")
      # # # cbind or rbind different lengths vectors without repeating the elements of the shorter vectors
      # length(dmrprobes) = maxlength
      
      dmrprobes = matrix(NA, 1, 0)
      for (j in 1: n_probes[i])
      {
        ## beta at column 2 
        # keep beta not probes_order[j, -2]
        dmrprobes = cbind(dmrprobes, probes_order[j, ])
        colnames(dmrprobes)[ (3*j - 2): (3*j )] = c(paste("Name", j, sep = "") ,
                                                    ##  beta means hyper or hypo direction
                                                    paste("direction", j, sep = ""),
                                                    paste("raw_p", j, sep = "") )
      }
      # build a single row df
      dmrprobes = dmrprobes %>% dplyr::mutate(cpg1 = probes_order$Name[1],
                                              positive =  sum(probes_order$beta > 0)/n_probes[i],
                                              No_probes = n_probes[i] ,
                                              maxbeta = max_beta,
                                              meanbeta = mean_beta) %>%
        dplyr::select(cpg1, direction1, 
                      maxbeta, meanbeta, positive, No_probes, everything())
      # Convert a row of a data frame to vector
      dmrprobes = unlist(dmrprobes)
      # cbind or rbind different lengths vectors without repeating the elements of the shorter vectors
      length(dmrprobes) = maxlength
  
    
      
    }
    # combine all dmrs together
    dmrprobes_all = data.frame(rbind(dmrprobes_all, dmrprobes))
    row.names(dmrprobes_all) = NULL
  }
  fulllist = cbind(dmr_sig, dmrprobes_all)
  write.xlsx(fulllist, 
             file =  paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/all_probes/", 
                                     sep = "_", name, pcutoff, "sigDMR_probes.xlsx" ) )
  ####### list focus on top cpg
  cpg1 = data.frame(fulllist) %>% dplyr::mutate(methylation1 = ifelse(direction1 > 0, 1, -1)) %>%
    dplyr::select(c("No.", "chrom", "start", "end",
                    "n_probes", "z_sidak_p", "maxbeta", "meanbeta", "refGene_name", 
                    "positive","cpg1"))
  
  write.xlsx(cpg1, 
             file =  paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/output_combp/", 
                                 sep = "_", pcutoff, name,"top1CpG.xlsx" ) )
  
}

p_cut <- 0.05
getprobes_perDMR(f_pfdea_dmr, f_pfdea_DMP , anno, p_cut)
getprobes_perDMR(f_pfna_dmr, f_pfna_DMP , anno, p_cut)
getprobes_perDMR(f_pfoa_dmr, f_pfoa_DMP , anno, p_cut)
getprobes_perDMR(f_pfos_dmr, f_pfos_DMP , anno, p_cut)
getprobes_perDMR(f_pfhxs_dmr, f_pfhxs_DMP , anno, p_cut)

getprobes_perDMR(m_pfdea_dmr, m_pfdea_DMP , anno, p_cut)
getprobes_perDMR(m_pfna_dmr, m_pfna_DMP , anno, p_cut)
getprobes_perDMR(m_pfoa_dmr, m_pfoa_DMP , anno, p_cut)
getprobes_perDMR(m_pfos_dmr, m_pfos_DMP , anno, p_cut)
getprobes_perDMR(m_pfhxs_dmr, m_pfhxs_DMP , anno, p_cut)
