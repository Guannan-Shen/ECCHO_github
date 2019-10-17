## set up workspace
rm(list = ls())
library(data.table)
library(tidyverse)
library(magrittr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

dir1 <-  "~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/" 
file_name1 <- c("2019-03-07_f_pfhxs__CpGs_withChem.csv",
                 "2019-03-07_f_pfoa__CpGs_withChem.csv",
                 "2019-03-07_f_pfos__CpGs_withChem.csv",
                 "2019-03-07_m_pfhxs__CpGs_withChem.csv", 
                 "2019-03-07_m_pfoa__CpGs_withChem.csv",
                 "2019-03-07_m_pfos__CpGs_withChem.csv")

file_name2 <- c("2019-09-30_f_pfdea__CpGs_withChem.csv",
                "2019-09-30_f_pfna__CpGs_withChem.csv",
                "2019-09-30_m_pfdea__CpGs_withChem.csv",
                "2019-09-30_m_pfna__CpGs_withChem.csv")

dir2 <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/"

file_name3 <- c("more_2019-10-17_f_pfdea__CpGs_withChem.csv",
                "more_2019-10-17_f_pfna__CpGs_withChem.csv",
                "more_2019-10-17_f_pfoa__CpGs_withChem.csv",
                "more_2019-10-17_f_pfos__CpGs_withChem.csv",
                "more_2019-10-17_f_pfhxs__CpGs_withChem.csv",
               
                "more_2019-10-17_m_pfdea__CpGs_withChem.csv",
                "more_2019-10-17_m_pfna__CpGs_withChem.csv",
                "more_2019-10-17_m_pfoa__CpGs_withChem.csv",
                "more_2019-10-17_m_pfos__CpGs_withChem.csv",
                "more_2019-10-17_m_pfhxs__CpGs_withChem.csv")

########### read rank anno fdr cutoff of DMPs ##################
dmp_read_rank_anno_save <- function(file, dir, fdrcut){
  # read
  data = fread(paste0(dir, file), header = T) 
  # 
  df = data[, list(ID, pos, CHR, betafc, raw)] %>% as.data.frame() %>% 
           dplyr::mutate(fdr_BH = p.adjust(raw, "BH") ) %>% 
    plyr::arrange(raw) %>%  dplyr::rename(beta_coef = betafc,
                                          start = pos)
  # get the annotation
  #   library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  # the full df of annotation
  anno = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)) %>% as.data.table() 
  sub_anno = anno[, list(Name, pos, UCSC_RefGene_Name, Relation_to_Island)] %>% as.data.frame()
  # merge
  df_anno = merge(sub_anno, df, by.x = "Name", by.y = "ID") %>% na.omit() %>% 
    plyr::arrange(raw) %>% dplyr::select(Name, start, CHR,  beta_coef, raw, 
                                         fdr_BH, UCSC_RefGene_Name, Relation_to_Island )
  # regular expression to select things before the ;
  df_anno$UCSC_RefGene_Name = sub(";.*", "", df_anno$UCSC_RefGene_Name)
  
  ############## significant dmps ###############
  df_sig = df_anno %>% dplyr::filter(fdr_BH < fdrcut)
  
  ############# save files #############
  write.csv(df_anno, paste0(dir, "rank_anno_sig/anno_", file), row.names = F)
  write.csv(df_sig, paste0(dir, "rank_anno_sig/sigFDR_", file), row.names = F)
  # return(list (DMP = df, DMP_ANNO = df_anno, DMP_sig = df_sig) )
}
# for (i in 1:length(file_name1)){
#   dmp_read_rank_anno_save(file_name1[i], dir1, fdrcut = 0.05) 
# }
# for (i in 1:length(file_name2)){
#   dmp_read_rank_anno_save(file_name2[i], dir2, fdrcut = 0.05) 
# }

for (i in 1:length(file_name3)){
  dmp_read_rank_anno_save(file_name3[i], dir2, fdrcut = 0.05)
}

##########################  fetch corresponding cpg beta coef and p of another gender ########3#
# "rank_anno_sig/sigFDR_"
# "rank_anno_sig/anno_"

join_sig_dmp <- function(dir, sig_filename, full_filename){
  chem1 = unlist(strsplit(full_filename, split = "_" ,fixed = TRUE) )[5]
  chem2 = unlist(strsplit(sig_filename, split = "_" ,fixed = TRUE) )[5]
  gen1 = unlist(strsplit(full_filename, split = "_" ,fixed = TRUE) )[4]
  gen2 = unlist(strsplit(sig_filename, split = "_" ,fixed = TRUE) )[4]
  
  if(chem1 != chem2)stop("Compare the same chemical!")
  if(gen1 == gen2)stop("Compare different Genders!")
  
  sig_list = read.csv(paste0(dir, sig_filename ))
  full_list = read.csv(paste0(dir, full_filename ))
  
  full_gender = ifelse(gen1 == "m", "Male", "Female")
  sub_full_list = full_list %>% dplyr::select(Name, beta_coef, raw, fdr_BH) %>%
    set_colnames(c("Name", 
                   paste0("beta_coef_", full_gender), 
                   paste0("raw_", full_gender), 
                   paste0("fdr_BH_", full_gender) )) 
  
  # plyr join Unlike merge, preserves the order of x no matter what join type is used.
  df = plyr::join(sig_list, sub_full_list, type = "inner", by = "Name")
  write.csv(df, paste0(dir, "join_", sig_filename), row.names = F)
  # return(df)
  
}

dir <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/rank_anno_sig/"

############ sig female #############333
join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_f_pfoa__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_m_pfoa__CpGs_withChem.csv"    )
join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_f_pfos__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_m_pfos__CpGs_withChem.csv"    )

join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_f_pfhxs__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_m_pfhxs__CpGs_withChem.csv"    )

join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_f_pfdea__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_m_pfdea__CpGs_withChem.csv"    )
join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_f_pfna__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_m_pfna__CpGs_withChem.csv"    )
############# sig male #################

join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_m_pfoa__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_f_pfoa__CpGs_withChem.csv"    )
join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_m_pfos__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_f_pfos__CpGs_withChem.csv"    )

join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_m_pfhxs__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_f_pfhxs__CpGs_withChem.csv"    )

join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_m_pfdea__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_f_pfdea__CpGs_withChem.csv"    )
join_sig_dmp(dir, 
             sig_filename = "sigFDR_more_2019-10-17_m_pfna__CpGs_withChem.csv",
             full_filename = "anno_more_2019-10-17_f_pfna__CpGs_withChem.csv"    )




