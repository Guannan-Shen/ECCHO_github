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


