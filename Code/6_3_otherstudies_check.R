library(data.table)
library(tidyverse)
library(magrittr)
library(openxlsx)
options(stringsAsFactors = F)
########### import dmp results ##################
f_pfoa_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_f_pfoa__CpGs_withChem.csv", header = T)
f_pfos_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_f_pfos__CpGs_withChem.csv", header = T)

m_pfoa_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_m_pfoa__CpGs_withChem.csv", header = T)
m_pfos_DMP <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/2019-03-07_m_pfos__CpGs_withChem.csv", header = T)

######## import top20 cpgs list ###############
dir <- "~/Documents/gitlab/ECCHO_github/DataRaw/3chem_otherstudies/"

# rename the first column to ID 
rename_ID <- function(data) {
  colnames(data)[1] = "ID"
  return(data.frame(data) )
}

pfoa_miura <- read.xlsx( paste(dir, "pofa_miura_cpg.xlsx", sep = "") ) %>% as.data.frame() %>% rename_ID
pfos_miura <- read.xlsx( paste(dir, "pfos_miura_cpg.xlsx", sep = "") ) %>% as.data.frame() %>% rename_ID
pfoa_kingsley <- read.xlsx( paste(dir, "pfoa_kingsley.xlsx", sep = "") ) %>% as.data.frame() %>% rename_ID

## merge by ID and return the ID and raw p and FDR (BH) by gender and chem
rawp_fdr <- function(dmp, list, gender) {
  # get the name 
  chem = (strsplit(as.character(substitute(list)), "_", fixed = T) %>% unlist )[1]
  study = (strsplit(as.character(substitute(list)), "_", fixed = T) %>% unlist )[2]
  
  data = plyr::join( list, dmp, by = "ID") %>% column_to_rownames("ID") %>%
    # beta fc is the direction of the study 
    dplyr::select ( raw, fdr, betafc ) %>% as.matrix() %>% round( ., 7) %>% as.data.frame() %>%
    rownames_to_column("ID")
  
  colnames(data)[2:4] = c( paste("Raw p-value in our study:", gender),
                           paste("FDR (BH) in our study:", gender),
                           paste("Beta value change in our study", gender))
  dir <- "~/Documents/gitlab/ECCHO_github/DataRaw/3chem_otherstudies/"
  write.csv(data.frame(data), paste(dir, chem, study, gender,".csv", sep = ""), row.names = F)
  return(data.frame(data) )
}

rawp_fdr( f_pfoa_DMP, pfoa_miura, "Female")
rawp_fdr( m_pfoa_DMP, pfoa_miura, "Male")

rawp_fdr( f_pfos_DMP, pfos_miura, "Female")
rawp_fdr( m_pfos_DMP, pfos_miura, "Male")

rawp_fdr( f_pfoa_DMP, pfoa_kingsley, "Female")
rawp_fdr( m_pfoa_DMP, pfoa_kingsley, "Male")



