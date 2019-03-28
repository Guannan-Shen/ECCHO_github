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
  write.csv(topM, paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/Mval",
                        sep = "_", name, ".csv"), row.names = F )
}
getM(top_f_pfhxs, dtall_f)
getM(top_f_pfos, dtall_f)
getM(top_f_pfoa, dtall_f)

getM(top_m_pfhxs, dtall_m)
getM(top_m_pfos, dtall_m)
getM(top_m_pfoa, dtall_m)


