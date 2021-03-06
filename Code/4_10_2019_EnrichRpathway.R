library(readxl)
library(knitr)
# load top1 cpg associated gene from FDR 0.05 DMR
####### import top cpg ################
top_f_pfhxs <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_f_pfhxs_dmr_0.05_topCpG.xlsx")
top_f_pfos <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_f_pfos_dmr_0.05_topCpG.xlsx")
top_f_pfoa <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_f_pfoa_dmr_0.05_topCpG.xlsx")

top_m_pfhxs <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_m_pfhxs_dmr_0.05_topCpG.xlsx")
top_m_pfos <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_m_pfos_dmr_0.05_topCpG.xlsx")
top_m_pfoa <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/topProbeDMR/_m_pfoa_dmr_0.05_topCpG.xlsx")
########### import top 300 cpgs from DMP ##########
f_pfhxs_DMP_300 <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/ECCHO_3Chems_fe_male_pathway/top300dmp/_f_pfhxs_DMP_300_.xlsx")
f_pfos_DMP_300 <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/ECCHO_3Chems_fe_male_pathway/top300dmp/_f_pfos_DMP_300_.xlsx")
f_pfoa_DMP_300 <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/ECCHO_3Chems_fe_male_pathway/top300dmp/_f_pfoa_DMP_300_.xlsx")

m_pfhxs_DMP_300 <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/ECCHO_3Chems_fe_male_pathway/top300dmp/_m_pfhxs_DMP_300_.xlsx")
m_pfos_DMP_300 <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/ECCHO_3Chems_fe_male_pathway/top300dmp/_m_pfos_DMP_300_.xlsx")
m_pfoa_DMP_300 <- read_excel("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/ECCHO_3Chems_fe_male_pathway/top300dmp/_m_pfoa_DMP_300_.xlsx")

########### significant dmps ################
dir <- 

## Enrich R pathway analysis  
library(enrichR)
library(openxlsx)

dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018",
         "KEGG_2019_Human")

getpathway_top <- function(CpGtable, database, gene_col){
  #########3 CpGtable a string and database a vector ############## 
  ######## gene_col a string
  name = CpGtable
  # run enrichr
  # the usage of call "$" and eval 
  pathways = enrichr( eval( call("$" ,as.name(CpGtable), gene_col) ), database)
  # data
  go_mol = data.frame( pathways[["GO_Molecular_Function_2018"]] ) 
  go_bio = data.frame( pathways[["GO_Biological_Process_2018"]] )
  go_cel = data.frame( pathways[["GO_Cellular_Component_2018"]] )
  kegg = data.frame( pathways[["KEGG_2019_Human"]] )
  
  ## save 
  write.xlsx(kegg, 
             paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/", sep = "_", name, "kegg_2019.xlsx" ) )
  write.xlsx(go_bio, 
             paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/", sep = "_", name, "go_biopro_2018.xlsx" ) )
  write.xlsx(go_mol, 
             paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/", sep = "_", name, "go_mole_2018.xlsx" ) )
  write.xlsx(go_cel, 
             paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/", sep = "_", name, "go_cell_2018.xlsx" ) )
  
  #return
  return( pathways )
}

top1_data <- c("top_f_pfhxs", "top_f_pfos","top_f_pfoa", 
               "top_m_pfhxs", "top_m_pfos","top_m_pfoa")
top300_data <- c("f_pfhxs_DMP_300", "f_pfos_DMP_300", "f_pfoa_DMP_300",
                 "m_pfhxs_DMP_300", "m_pfos_DMP_300", "m_pfoa_DMP_300")

top1_gene <- c("refGene_name")
top300_gene <- c("UCSC_RefGene_Name")
## top1 cpg dmr
f_pfhxs_1 <- getpathway_top( "top_f_pfhxs",  dbs, top1_gene)
kable(f_pfhxs_1$KEGG_2019_Human[1:10, ], caption = "f_pfhxs Top1 per significant DMR" )
f_pfhxs_300 <- getpathway_top( "f_pfhxs_DMP_300",  dbs, top300_gene)

########### significant dmps ################
library(knitr)
## Enrich R pathway analysis  
library(enrichR)
library(openxlsx)
library(tidyverse)

dbs <- c( "GO_Cellular_Component_2018",
         "KEGG_2019_Human")

dir <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/newdmps/"

f_pfhxs_sigDMP <- read.csv(paste0(dir, "_2019-07-01 18:36:07_Female_PFHXS_0.05_fdr_dmps.csv"))
f_pfos_sigDMP <- read.csv(paste0(dir, "_2019-07-01 18:36:15_Female_PFOS_0.05_fdr_dmps.csv"))
f_pfoa_sigDMP <- read.csv(paste0(dir, "_2019-07-01 18:36:11_Female_PFOA_0.05_fdr_dmps.csv"))

m_pfhxs_sigDMP <- read.csv(paste0(dir, "_2019-07-01 18:36:19_Male_PFHXS_0.05_fdr_dmps.csv"))
m_pfos_sigDMP <- read.csv(paste0(dir, "_2019-07-01 18:36:26_Male_PFOS_0.05_fdr_dmps.csv"))
m_pfoa_sigDMP <- read.csv(paste0(dir, "_2019-07-01 18:36:22_Male_PFOA_0.05_fdr_dmps.csv"))

sig_gene <- c("UCSC_RefGene_Name")
sig_data <- c("f_pfhxs_sigDMP", "f_pfos_sigDMP", "f_pfoa_sigDMP",
                 "m_pfhxs_sigDMP", "m_pfos_sigDMP", "m_pfoa_sigDMP")

getpathway_2 <- function(CpGtable, database, gene_col){
  #########3 CpGtable a string and database a vector ############## 
  ######## gene_col a string
  name = CpGtable
  # run enrichr
  # the usage of call "$" and eval 
  pathways = enrichr( eval( call("$" ,as.name(CpGtable), gene_col) ), database)
  # data
  go_cel = data.frame( pathways[["GO_Cellular_Component_2018"]] ) 
             

  
  kegg = data.frame( pathways[["KEGG_2019_Human"]] ) 

  
  ## save 
  write.xlsx(kegg, 
             paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/SigDMP/", sep = "_", name, "kegg_2019.xlsx" ) )
  write.xlsx(go_cel, 
             paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/SigDMP/", sep = "_", name, "go_cell_2018.xlsx" ) )
  
  #return
  return( pathways )
}

#run and save results
getpathway_2 ("m_pfoa_sigDMP", dbs, sig_gene)

%>%
  dplyr::filter(Adjusted.P.value <= 0.2)
%>% 
  dplyr::filter(Adjusted.P.value <= 0.2)

overlap_go_cel = matrix(unlist(strsplit(go_cel$Overlap, "/", fixed = TRUE)), ncol = 2, byrow = T) [ , 1]
go_cel_3 = go_bio[ overlap_go_cel >= 2, ]


overlap_kegg = matrix(unlist(strsplit(kegg$Overlap, "/", fixed = TRUE)), ncol = 2, byrow = T) [ , 1]
kegg_3 = go_bio[ overlap_kegg >= 2, ]


########## filter existing results by overlapping more than 2 genes and fdr < 0.2 ##########
library(tidyverse)

sig_data <- c("f_pfhxs_sigDMP", "f_pfos_sigDMP", "f_pfoa_sigDMP",
              "m_pfhxs_sigDMP", "m_pfos_sigDMP", "m_pfoa_sigDMP")

dir <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/newdmps/"
sig_m_pfoa_kegg <- read.table(paste0(dir, "sig_m_pfoa_KEGG_2019_Human_table.txt"),
                              header = T, sep = "\t") %>% as.data.frame()
kegg <-  sig_m_pfoa_kegg %>%
  dplyr::filter(Adjusted.P.value <= 0.2) %>% dplyr::arrange(P.value)

write.xlsx(kegg, 
           paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/SigDMP/", sep = "_", sig_data[6], "kegg_2019.xlsx" ) )

sig_m_pfoa_gocell <- read.table(paste0(dir, "sig_m_pfoa_GO_Cellular_Component_2018_table.txt"),
                              header = T, sep = "\t") %>% as.data.frame()
go_cell <- sig_m_pfoa_gocell %>%
  dplyr::filter(Adjusted.P.value <= 0.2) %>% dplyr::arrange(P.value)

write.xlsx(go_cell, 
           paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/pathways_eccho_3chem/SigDMP/", sep = "_", sig_data[6], "go_cell_2018.xlsx" ) )



########################## More pfas ###########################################
######################### enrichment analysis and format #################
library(readxl)
library(knitr)

########### significant dmps ################
dir <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/rank_anno_sig/"
  
## Enrich R pathway analysis  
library(enrichR)
library(openxlsx)
library(tidyverse)
library(magrittr)

dbs <- c( "GO_Cellular_Component_2018", "KEGG_2019_Human")

sig_gene <- c("UCSC_RefGene_Name")
sig_data <- c("f_pfdea_sigDMP", "f_pfna_sigDMP", "f_pfhxs_sigDMP", "f_pfos_sigDMP", "f_pfoa_sigDMP",
              "m_pfdea_sigDMP", "m_pfna_sigDMP","m_pfhxs_sigDMP", "m_pfos_sigDMP", "m_pfoa_sigDMP")

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


getpathway_sig <- function(CpGtable, database, gene_col, fdrcut, overlap_cut){
  #########3 CpGtable a string and database a vector ############## 
  ######## gene_col a string
  name = CpGtable
  # run enrichr
  # the usage of call "$" and eval 
  pathways = enrichr( eval( call("$" ,as.name(CpGtable), gene_col) ), database)
  # data
  go_cel = data.frame( pathways[["GO_Cellular_Component_2018"]] ) %>% 
    dplyr::select(-Old.P.value, - Old.Adjusted.P.value) %>% plyr::arrange(P.value) %>%
    dplyr::mutate(FDR_BH = p.adjust(P.value, "BH"))  %>% dplyr::select(-Adjusted.P.value)
  write.csv(go_cel, 
            paste0("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/enrichment/GO_CELL_",
                   CpGtable, ".csv" ), row.names = F)
  
  overlap_go_cel = matrix(unlist(strsplit(go_cel$Overlap, "/", fixed = TRUE)), ncol = 2, byrow = T) [ , 1]
  go_cel2 = go_cel[ overlap_go_cel >= overlap_cut, ]
  
  go_cel_3 = go_cel2 %>%  dplyr::filter(FDR_BH <= fdrcut) 
  
  if (nrow(go_cel_3) != 0){
  write.csv(go_cel_3, 
            paste0("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/enrichment/Filtered_GO_CELL_",
                   CpGtable, ".csv" ), row.names = F)
  }
  
  kegg = data.frame( pathways[["KEGG_2019_Human"]] ) %>% 
         dplyr::select(-Old.P.value, - Old.Adjusted.P.value) %>% plyr::arrange(P.value) %>%
                  dplyr::mutate(FDR_BH = p.adjust(P.value, "BH")) %>% dplyr::select(-Adjusted.P.value)
  write.csv(kegg, 
            paste0("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/enrichment/kegg_",
            CpGtable, ".csv" ), row.names = F)
  
  overlap_kegg = matrix(unlist(strsplit(kegg$Overlap, "/", fixed = TRUE)), ncol = 2, byrow = T) [ , 1]
  kegg2 = kegg[ overlap_kegg >= overlap_cut, ]
  kegg_3 = kegg2 %>%  dplyr::filter(FDR_BH <= fdrcut) 
  
  if (nrow(kegg_3) != 0){
  write.csv(kegg_3, 
            paste0("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/enrichment/Filtered_KEGG_",
                   CpGtable, ".csv" ), row.names = F)
  }
  # return(kegg)
}

# test
eval( call("$" ,as.name("f_pfdea_sigDMP"), "UCSC_RefGene_Name") )
enrichr(f_pfhxs_sigDMP$UCSC_RefGene_Name , dbs)

getpathway_sig ("f_pfdea_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("f_pfna_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("f_pfhxs_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("f_pfoa_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("f_pfos_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)


getpathway_sig ("m_pfdea_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("m_pfna_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("m_pfhxs_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("m_pfoa_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
getpathway_sig ("m_pfos_sigDMP", dbs, sig_gene, fdrcut = 0.2, overlap_cut = 2)
