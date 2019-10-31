######### get mval per gender from whole ##########
library(data.table)
library(readxl)
library(tidyverse)

####### import top cpg ################
top_f_pfhxs <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_f_pfhxs_dmr_0.05_topCpG.xlsx")
top_f_pfos <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_f_pfos_dmr_0.05_topCpG.xlsx")
top_f_pfoa <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_f_pfoa_dmr_0.05_topCpG.xlsx")

top_m_pfhxs <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_m_pfhxs_dmr_0.05_topCpG.xlsx")
top_m_pfos <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_m_pfos_dmr_0.05_topCpG.xlsx")
top_m_pfoa <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_m_pfoa_dmr_0.05_topCpG.xlsx")
dim(top_m_pfos)

################## more_pfas ######################
######### get top1 cpg per DMR ############
dir2 <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/output_combp/"
top_f_pfdea <- read_excel(paste0(dir2, "_0.05_f_pfdea_dmr_0.05_top1CpG.xlsx"))
top_f_pfna <- read_excel(paste0(dir2, "_0.05_f_pfna_dmr_0.05_top1CpG.xlsx"))
top_f_pfoa <- read_excel(paste0(dir2, "_0.05_f_pfoa_dmr_0.05_top1CpG.xlsx"))
top_f_pfos <- read_excel(paste0(dir2, "_0.05_f_pfos_dmr_0.05_top1CpG.xlsx"))
top_f_pfhxs <- read_excel(paste0(dir2, "_0.05_f_pfhxs_dmr_0.05_top1CpG.xlsx"))


top_m_pfdea <- read_excel(paste0(dir2, "_0.05_m_pfdea_dmr_0.05_top1CpG.xlsx"))
top_m_pfna <- read_excel(paste0(dir2, "_0.05_m_pfna_dmr_0.05_top1CpG.xlsx"))
top_m_pfoa <- read_excel(paste0(dir2, "_0.05_m_pfoa_dmr_0.05_top1CpG.xlsx"))
top_m_pfos <- read_excel(paste0(dir2, "_0.05_m_pfos_dmr_0.05_top1CpG.xlsx"))
top_m_pfhxs <- read_excel(paste0(dir2, "_0.05_m_pfhxs_dmr_0.05_top1CpG.xlsx"))
########## m val #############
## get mval 
# mval
dtall_f <- fread( "~/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/dt_all_f.csv", header = T)
# mval of male 
clinchem <- read.csv( "/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/clin_chem.csv",
                      header = T)
t_mval <- fread( "~/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/t_mval.csv",
                 header = T)
dt_all <- merge(clinchem, t_mval, by = "pid")

dim(dt_all)
dt_all <- data.table(dt_all)
dtall_m <-  dt_all[infant_sex == "Male"]
dim(dtall_m)

##### get m for topcpg
getM <- function(top_chem, datatable){
  cpgs = top_chem$cpg1
  name = substitute(top_chem)
  topM = datatable[, cpgs, with = FALSE]
  topM = data.frame(topM) %>% dplyr::mutate(pid = datatable$pid) %>% select(pid, everything())
  write.csv(topM, paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/Mval/",
                        sep = "_", name, ".csv"), row.names = F )
}

getM(top_f_pfhxs, dtall_f)
getM(top_f_pfos, dtall_f)
getM(top_f_pfoa, dtall_f)

getM(top_m_pfhxs, dtall_m)
getM(top_m_pfos, dtall_m)
getM(top_m_pfoa, dtall_m)

######## two more chems , mval of top1 CpG per DMR ##########
getM(top_f_pfdea, dtall_f)
getM(top_f_pfna, dtall_f)

getM(top_m_pfdea, dtall_m)
getM(top_m_pfna, dtall_m)

#################  get mvale for fdr 005 dmps ##############
library(data.table)
library(tidyverse)
# import cpgs which are fdr 005 from DMP analysis
dir <- "~/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/fdr005_dmps/"
top_f_pfhxs <- read.csv( paste0(dir, "Female_PFHXS_0.05_fdr_dmps.csv") )
top_f_pfos <- read.csv( paste0(dir, "Female_PFOS_0.05_fdr_dmps.csv") )
top_f_pfoa <- read.csv( paste0(dir, "Female_PFOA_0.05_fdr_dmps.csv") )

# male female seperately
top_m_pfhxs <- read.csv( paste0(dir, "Male_PFHXS_0.05_fdr_dmps.csv") )
top_m_pfos <- read.csv( paste0(dir, "Male_PFOS_0.05_fdr_dmps.csv") )
top_m_pfoa <- read.csv( paste0(dir, "Male_PFOA_0.05_fdr_dmps.csv") )

############### get DMP sig cpgs ##############
dir <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/rank_anno_sig/"
f_pfdea_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_f_pfdea__CpGs_withChem.csv"))
f_pfna_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_f_pfna__CpGs_withChem.csv"))
f_pfhxs_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_f_pfhxs__CpGs_withChem.csv"))
f_pfos_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_f_pfos__CpGs_withChem.csv"))
f_pfoa_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_f_pfoa__CpGs_withChem.csv"))


m_pfdea_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_m_pfdea__CpGs_withChem.csv"))
m_pfna_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_m_pfna__CpGs_withChem.csv"))
m_pfhxs_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_m_pfhxs__CpGs_withChem.csv"))
m_pfos_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_m_pfos__CpGs_withChem.csv"))
m_pfoa_sigDMP <- read.csv(paste0(dir, "sigFDR_more_2019-10-17_m_pfoa__CpGs_withChem.csv"))


## get mval 
# mval
dtall_f <- fread( "~/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/dt_all_f.csv", header = T)

getM <- function(top_chem, datatable){
  cpgs = top_chem$Name
  name = substitute(top_chem)
  # – Select columns named in a variable using with = FALSE
  topM = datatable[, cpgs, with = FALSE]
  topM = data.frame(topM) %>% dplyr::mutate(pid = datatable$pid) %>% select(pid, everything())
  write.csv(topM, 
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/rank_anno_sig/Mval",
                        sep = "_", name, ".csv"), row.names = F )
}

getM(f_pfdea_sigDMP, dtall_f)
getM(f_pfna_sigDMP, dtall_f)
getM(f_pfhxs_sigDMP, dtall_f)
getM(f_pfoa_sigDMP, dtall_f)
getM(f_pfos_sigDMP, dtall_f)

# double check 
Mtop_f_pfhxs <- read.csv(paste0(dir, "Mval_top_f_pfhxs_.csv"))
ncol(Mtop_f_pfhxs) == nrow(top_f_pfhxs) + 1

# 
# mval of male 
clinchem <- read.csv( "/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/clin_chem.csv",
                      header = T)
t_mval <- fread( "~/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/t_mval.csv",
                 header = T)
dt_all <- merge(clinchem, t_mval, by = "pid")

dim(dt_all)
dt_all <- data.table(dt_all)
dtall_m <-  dt_all[infant_sex == "Male"]
dim(dtall_m)

getM(top_m_pfhxs, dtall_m)
getM(top_m_pfos, dtall_m)
getM(top_m_pfoa, dtall_m)

getM(m_pfdea_sigDMP, dtall_m)
getM(m_pfna_sigDMP, dtall_m)
getM(m_pfhxs_sigDMP, dtall_m)
getM(m_pfoa_sigDMP, dtall_m)
getM(m_pfos_sigDMP, dtall_m)


# double check 
Mtop_m_pfhxs <- read.csv(paste0(dir, "Mval_top_m_pfhxs_.csv"))
ncol(Mtop_m_pfhxs) == nrow(top_m_pfhxs) + 1

