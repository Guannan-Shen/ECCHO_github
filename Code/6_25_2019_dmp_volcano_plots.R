# Manhattan/volcano plots for 
# the association between PFASOA and DNA methylation 
# among males and females separately.
# from dmp results to volcano plot
library(data.table)
library(tidyverse)
library(magrittr)
# library for manhattan plot
library(qqman)
# "geom_text_repel"
require(ggrepel)
options(stringsAsFactors = F)
getwd()
# notes
print("Basically, using the fdr column from the dmrcate DMP analysis results!")
input_names <- c("2019-03-07_f_pfhxs__CpGs_withChem.csv",
                 "2019-03-07_f_pfoa__CpGs_withChem.csv",
                 "2019-03-07_f_pfos__CpGs_withChem.csv",
                 "2019-03-07_m_pfhxs__CpGs_withChem.csv", 
                 "2019-03-07_m_pfoa__CpGs_withChem.csv",
                 "2019-03-07_m_pfos__CpGs_withChem.csv")
# dir = "C:/Users/hithr/Documents/Stats/gitlab/ECCHO_github/DataProcessed/genomewide_chem/"
# Ubuntu
dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/" 
# dir file name 
# function readin dmp results .xlsx file and reformat to betafc rawp value and cpg name 
dmp_read <- function(dmpresult_name){
  # asus
  # dir = "C:/Users/hithr/Documents/Stats/gitlab/ECCHO_github/DataProcessed/genomewide_chem/"
  # Ubuntu
  dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/" 
  data = fread(paste0(dir, dmpresult_name), header = T) 
  temp = data[, BP := seq_len(.N), by = CHR]
  # get chrom in order
  chrOrder = c(paste("chr",1:22,sep=""),"chrX","chrY")
  df = data[, list(ID, CHR, betafc, raw, BP)] %>% as.data.frame() %>% 
                           dplyr::mutate(beta10 = 10*betafc,
                                         CHR = factor(CHR, levels=chrOrder, labels = 1:24)) %>% 
                                 arrange(CHR) %>% 
    # the qqman package requires numeric CHR not x, y, 
                           dplyr::mutate(CHR = as.numeric(CHR))
  
  return(df)
}
#
df = dmp_read(dmpresult_name = input_names[1])
# table(df$CHR)
# df$CHR
# test if the fdr column can be calculated from raw
# test = dmp_read("2019-03-07_m_pfhxs__CpGs_withChem.csv")
# test[p.adjust(test$raw, "BH") != test$fdr, ]
## get annotation 
#get annotation and merge to all dmps
dmp_anno_read <- function(dmpresult_name){
  # asus
  # dir = "C:/Users/hithr/Documents/Stats/gitlab/ECCHO_github/DataProcessed/genomewide_chem/"
  # Ubuntu
  dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/genomewide_chem/" 
  data = fread(paste0(dir, dmpresult_name), header = T) 
  temp = data[, BP := seq_len(.N), by = CHR]
  # get chrom in order
  chrOrder = c(paste("chr",1:22,sep=""),"chrX","chrY")
  df = data[, list(ID, CHR, betafc, raw, BP)] %>% as.data.frame() %>% 
    dplyr::mutate(beta10 = 10*betafc,
                  CHR = factor(CHR, levels=chrOrder, labels = 1:24)) %>% 
    arrange(CHR) %>% 
    # the qqman package requires numeric CHR not x, y, 
    dplyr::mutate(CHR = as.numeric(CHR))
  # merge df with annotation
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)) %>% as.data.table() 
  sub_anno = anno[, list(Name, UCSC_RefGene_Name)] %>% as.data.frame()
  df_anno = merge(sub_anno, df, by.x = "Name", by.y = "ID") %>% na.omit()
  # regular expression to select things before the ;
  df_anno$UCSC_RefGene_Name = sub(";.*", "", df_anno$UCSC_RefGene_Name)
  return(df_anno)
}

dmp_vol <- function(dmpresult_name, p_cutoff_label, fc_cutoff_label){
  ############ get names #############
  g_name = unlist(strsplit( dmpresult_name, "_"))[2]
  # from f, m to Female and Male 
  gender = ifelse(g_name == "f", "Female", "Male")
  # upper case
  chem_name = toupper(unlist(strsplit( dmpresult_name, "_"))[3])
  name = paste(gender, chem_name)
  ### load in dmp results ###
  df = dmp_read(dmpresult_name)
  ######### volcano plot ######
  # ID, betafc, raw
  p = ggplot(df, aes(beta10, -log10(raw) )) +
    geom_point(colour = "gray70") + 
    geom_point(data = subset(df, raw < 5e-2 & beta10 > 0 ),  aes(beta10, -log10(raw) ), colour="black") +
    geom_point(data = subset(df, raw < 5e-2 & beta10 < 0 ),  aes(beta10, -log10(raw) ), colour="gray40") +
    labs(x = paste0("Fold Change in Beta-values per 10 units increase of ", chem_name ),
         y =  expression(paste("-",log[10]," Raw p-value" )),
         caption =  paste0("Genomewide test of associations between ", 
                           chem_name, " and DNA methylation for ", gender) ) +
    # text annotation with CpGs ID
    geom_text_repel(data=subset(df, (raw <= p_cutoff_label) & (abs(beta10) >= fc_cutoff_label) ),
                    aes(label=ID), size = 3,
                    box.padding = unit(0.6, "lines"),
                    point.padding = unit(0.6, "lines")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    annotate(geom="text", x = 0.75, 
             y=-log10(0.02), 
             label="Raw p-value < 0.05",
             color="black") +
    # range of x
    xlim(-1.0, 1.0) +
    theme_bw()
  print(p)
  # save plots
  dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/plots/genome_wide/" 
  ## ggsave with 
  ggsave(paste0(dir, name, "_cpgs.tiff"), dpi=300, compression = "lzw")
}


# table(gwasResults$BP)

gene_vol <- function(dmpresult_name, p_cutoff_label, fc_cutoff_label){
  ############ get names #############
  g_name = unlist(strsplit( dmpresult_name, "_"))[2]
  # from f, m to Female and Male 
  gender = ifelse(g_name == "f", "Female", "Male")
  # upper case
  chem_name = toupper(unlist(strsplit( dmpresult_name, "_"))[3])
  name = paste(gender, chem_name)
  ### load in dmp results ###
  df = dmp_anno_read(dmpresult_name)
  ######### volcano plot ######
  # ID, betafc, raw
  p = ggplot(df, aes(beta10, -log10(raw) )) +
    geom_point(colour = "gray70") + 
    geom_point(data = subset(df, raw < 5e-2 & beta10 > 0 ),  aes(beta10, -log10(raw) ), colour="black") +
    geom_point(data = subset(df, raw < 5e-2 & beta10 < 0 ),  aes(beta10, -log10(raw) ), colour="gray40") +
    labs(x = paste0("Fold Change in Beta-values per 10 units increase of ", chem_name ),
         y =  expression(paste("-",log[10]," Raw p-value" )),
         caption =  paste0("Genomewide test of associations between ", 
                           chem_name, " and DNA methylation for ", gender) ) +
    # text annotation with CpGs ID
    geom_text_repel(data=subset(df, (raw <= p_cutoff_label) & (abs(beta10) >= fc_cutoff_label) ),
                    aes(label=UCSC_RefGene_Name), size = 3,
                    box.padding = unit(0.6, "lines"),
                    point.padding = unit(0.6, "lines")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    annotate(geom="text", x = 0.75, 
             y=-log10(0.02), 
             label="Raw p-value < 0.05",
             color="black") +
    # range of x
    xlim(-1.0, 1.0) +
    theme_bw()
  print(p)
  # save plots
  dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/plots/genome_wide/" 
  ## ggsave with 
  ggsave(paste0(dir, name, "_genes.tiff"), dpi=300, compression = "lzw")
}

# no text annotation
dmp_vol_clean <- function(dmpresult_name){
  ############ get names #############
  g_name = unlist(strsplit( dmpresult_name, "_"))[2]
  # from f, m to Female and Male 
  gender = ifelse(g_name == "f", "Female", "Male")
  # upper case
  chem_name = toupper(unlist(strsplit( dmpresult_name, "_"))[3])
  name = paste(gender, chem_name)
  ### load in dmp results ###
  df = dmp_read(dmpresult_name)
  ######### volcano plot ######
  # ID, betafc, raw
  p = ggplot(df, aes(beta10, -log10(raw) )) +
    geom_point(colour = "gray70") + 
    geom_point(data = subset(df, raw < 5e-2 & beta10 > 0 ),  aes(beta10, -log10(raw) ), colour="black") +
    geom_point(data = subset(df, raw < 5e-2 & beta10 < 0 ),  aes(beta10, -log10(raw) ), colour="gray40") +
    labs(x = paste0("Fold Change in Beta-values per 10 units increase of ", chem_name ),
         y =  expression(paste("-",log[10]," Raw p-value" )),
         caption =  paste0("Genomewide test of associations between ", 
                           chem_name, " and DNA methylation for ", gender) ) +
    # text annotation with CpGs ID
    # geom_text_repel(data=subset(df, (raw <= p_cutoff_label) & (abs(beta10) >= fc_cutoff_label) ),
    #                 aes(label=ID), size = 3,
    #                 box.padding = unit(0.1, "lines"),
    #                 point.padding = unit(0.1, "lines")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    annotate(geom="text", x = 0.75, 
             y=-log10(0.02), 
             label="Raw p-value < 0.05",
             color="black") +
    # range of x
    xlim(-1.0, 1.0) +
    theme_bw()
  print(p)
  # save plots
  dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/plots/genome_wide/" 
  ## ggsave with 
  ggsave(paste0(dir, name, "_cpgs_clean.tiff"), dpi=300, compression = "lzw")
}

gene_vol_clean <- function(dmpresult_name){
  ############ get names #############
  g_name = unlist(strsplit( dmpresult_name, "_"))[2]
  # from f, m to Female and Male 
  gender = ifelse(g_name == "f", "Female", "Male")
  # upper case
  chem_name = toupper(unlist(strsplit( dmpresult_name, "_"))[3])
  name = paste(gender, chem_name)
  ### load in dmp results ###
  df = dmp_anno_read(dmpresult_name)
  ######### volcano plot ######
  # ID, betafc, raw
  p = ggplot(df, aes(beta10, -log10(raw) )) +
    geom_point(colour = "gray70") + 
    geom_point(data = subset(df, raw < 5e-2 & beta10 > 0 ),  aes(beta10, -log10(raw) ), colour="black") +
    geom_point(data = subset(df, raw < 5e-2 & beta10 < 0 ),  aes(beta10, -log10(raw) ), colour="gray40") +
    labs(x = paste0("Fold Change in Beta-values per 10 units increase of ", chem_name ),
         y =  expression(paste("-",log[10]," Raw p-value" )),
         caption =  paste0("Genomewide test of associations between ", 
                           chem_name, " and DNA methylation for ", gender) ) +
    # text annotation with CpGs ID
    # geom_text_repel(data=subset(df, (raw <= p_cutoff_label) & (abs(beta10) >= fc_cutoff_label) ),
    #                 aes(label=UCSC_RefGene_Name), size = 3,
    #                 box.padding = unit(0.1, "lines"),
    #                 point.padding = unit(0.1, "lines")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    annotate(geom="text", x = 0.75, 
             y=-log10(0.02), 
             label="Raw p-value < 0.05",
             color="black") +
    # range of x
    xlim(-1.0, 1.0) +
    theme_bw()
  print(p)
  # save plots
  dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/plots/genome_wide/" 
  ## ggsave with 
  ggsave(paste0(dir, name, "_genes_clean.tiff"), dpi=300, compression = "lzw")
}

## 

# make plots 
########## be careful #################
# for(i in input_names){
#   # dmp volcano plot
#   dmp_vol(dmpresult_name = i, p_cutoff_label = 1e-5, fc_cutoff_label = 4e-1)
#   # gene volcano plot
#   gene_vol(dmpresult_name = i, p_cutoff_label = 1e-5, fc_cutoff_label = 4e-1)
#   #
#   dmp_vol_clean(dmpresult_name = i)
#   gene_vol_clean(dmpresult_name = i)
# }
# 
# ###### changing from ggplot to basic plot #############
# dev.off()
# 
# ## manhattan plot 
# dmp_man <- function(dmpresult_name){
#   ############ get names #############
#   g_name = unlist(strsplit( dmpresult_name, "_"))[2]
#   # from f, m to Female and Male 
#   gender = ifelse(g_name == "f", "Female", "Male")
#   # upper case
#   chem_name = toupper(unlist(strsplit( dmpresult_name, "_"))[3])
#   name = paste(gender, chem_name)
#   ### load in dmp results ###
#   df = dmp_read(dmpresult_name)
#   ## Manhattan plot ##
#   # save plots
#   # dir = "~/Documents/gitlab/ECCHO_github/DataProcessed/3chems_results_data/plots/genome_wide/" 
#   # par( mfrow = c(1,1) )
#   # tiff(paste0(dir, name, "_Manhattan.tiff"), res = 300, compression = "lzw")
#    manhattan(df, chr="CHR",bp="BP", snp="ID", p="raw", cex = 0.6, ymax = 14, 
#                suggestiveline = F, genomewideline = F, 
#               chrlabs = c(1:22, "X", "Y"))
#    title(sub = paste0("Genomewide test of associations between ", 
#                       chem_name, " and DNA methylation for ", gender))
#   #  dev.off()
# }
# 
# 
# # dmp manhattan plot
# # for(i in input_names){
# # dmp_man(dmpresult_name = i)
# # }
# 
# dmp_man(dmpresult_name = input_names[2])
