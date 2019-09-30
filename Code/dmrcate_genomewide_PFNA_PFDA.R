rm(list = ls())
options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
## not in function
'%nin%' <- Negate('%in%')

library(data.table)
library(tidyverse)
library(magrittr)
library(stats)
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tibble)
library(grid) # low-level grid functions are required

################ function wrapper for DMR analysis ###################
DMRcate_wrapper <- function(formula, clindat, mval){
  ## dataset use , female or male, lnpfoa ##
  gender1 = unlist(strsplit( deparse(substitute(clindat)), "_"))[2]
  gender2 = unlist(strsplit( deparse(substitute(mval)), "_"))[2]
  if(gender1 != gender2)
    stop("should use the same gender data")
  ## the type of pfoa
  chem = unlist(strsplit( deparse(substitute(formula)), "_"))[2]
  prefix_name = paste(gender1,"_", chem,sep = "")
  ### build design matrix ###
  str(model <- model.frame(formula, clindat))
  design = model.matrix(formula, model)
  
  ######## DMR analysis ########
  ## Build cpg anno object
  cpganno = cpg.annotate("array",
                         ## female 
                         # with unique Illumina probe IDs as rownames and unique
                         # sample IDs as column names o
                         mval,
                         what = "M",
                         arraytype = "450K",
                         analysis.type = "differential",
                         ## female
                         design = design,
                         ## only have the intercept and the pfoa concentration
                         ##  it will still be coef = 2. 
                         ## We are telling DMRcate that we are interested in the log transformed pfoa.
                         coef = 2
                         # Highly recommended as the primary thresholding parameter for calling DMRs 
                         # fdr = 0.05
  )
  ## get results ##
  dmrcoutput = dmrcate(cpganno, 
                       lambda=1000, 
                       c=2, 
                       p.adjust.method="BH", 
                       consec=FALSE, 
                       pcutoff = 0.05,
                       mc.cores = 4)
  #######################################
  # Get list of significant DMRs and CpGs
  # associated with those DMRs
  #######################################
  #get annotation;
  anno = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  
  #get probes in DMR regions;
  getDMRprobes = function(a, anno){
    chr = sapply(strsplit(as.character(a[1]), split=":", fixed=TRUE), "[[", 1)
    start =  as.numeric(sapply(strsplit(sapply(strsplit(as.character(a[1]), split=":", fixed=TRUE), "[[", 2), split="-", fixed=TRUE), "[[", 1))
    stop =  as.numeric(sapply(strsplit(as.character(a[1]), split="-", fixed=TRUE), "[[", 2))
    
    anno.chr = anno[which(anno$chr==chr),] 
    anno2 = as.data.frame(anno.chr[which(anno.chr$pos >= start & anno.chr$pos <= stop),])
    
    return(anno2) 
  }
  dmr.results = dmrcoutput$results
  ## save dmrcate results
  write.csv(dmrcoutput$results, 
            row.names = F,
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate/",  
                  Sys.Date(),"_", 
                  prefix_name,"__DMRcate_DMR", ".csv", sep= ""))
  ## save genomewide results
  write.csv(dmrcoutput$input, 
            row.names = F,
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate_genome/",  
                  Sys.Date(),"_", 
                  prefix_name,"__CpGs_withChem", ".csv", sep= ""))
  ########### get detailed DMR level results ############### 
  ##Build list of annotation data for probes in each DMR
  DMRprobes = list()
  for(i in 1:nrow(dmr.results)){
    DMRprobes[[i]] = getDMRprobes(dmr.results[i,], anno)
  }
  
  names(DMRprobes) = dmr.results$coord
  names(DMRprobes) = gsub(":", "_", names(DMRprobes))
  # Write results
  readme_sheet = data_frame(
    Columns = c(
      "Annotations for CpGs assocaited with DMR's unadjusted p-value < 0.05",
      "", colnames(DMRprobes[[1]])
    ))
  
  readme_sheet = list(README = readme_sheet)
  names(readme_sheet) = "README"
  openxlsx::write.xlsx(c(readme_sheet, 
                         DMRprobes),
                       paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate/",  
                             Sys.Date(),"_", 
                             prefix_name,"_DMR", ".xlsx", sep= ""))
  
  ############# get individual cpgs level results data frame ################
  n_cpg = sum(dmrcoutput$results$no.cpgs)
  cpg_sum = (data.frame(matrix(NA, 0,4)))
  ## generate summary table
  for(i in 1:nrow(dmr.results)){
    cpg_sum  = rbind(cpg_sum, 
                     cbind(dmr.results[i,1], 
                           dmrcoutput$input [dmrcoutput$input$ID %in% getDMRprobes(dmr.results[i,], anno)[,4], c(1,3,8)]) )
  }
  colnames(cpg_sum)[c(1,4)] = c("DMR_Identifier", "raw_p")
  if(nrow(cpg_sum) != n_cpg)
    stop("wrong total number of cpgs")
  cpg_min = cpg_sum %>% group_by(DMR_Identifier) %>% filter(raw_p == min(raw_p))
  if(nrow(cpg_min) != nrow(dmr.results))
    stop("wrong number of top1 cpgs in each dmr")
  write.csv(cpg_sum, 
            row.names = F,
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate/",  
                  Sys.Date(),"_", 
                  prefix_name,"_DMR_allCpGs", ".csv", sep= ""))
  write.csv(cpg_min, 
            row.names = F,
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/dmrcate/",  
                  Sys.Date(),"_", 
                  prefix_name,"_DMR_top1_CpG", ".csv", sep= ""))
}

#################### RUN #########################
## formula
formula_pfdea <-  ~as.numeric(lnpfdea) + Race +  Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC
formula_pfna <-  ~as.numeric(lnpfna) + Race +  Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC


#### all female
# all clinical data, chemical conc, obesity outcomes
clinchem_f <- read.csv("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/clin_chem_f_nada.csv")
dim(clinchem_f)
# mval
dtall_f <- fread( "~/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/dt_all_f.csv", header = T)
dim(dtall_f)
# dtall_f <- data.frame(dtall_f)
cpgs <- colnames(dtall_f)[-c(1:26)]
# all pids were aligned 
sum(dtall_f$pid != clinchem_f$pid)

##
tdtall_f <- t(dtall_f[,-c(1:26)])
## check pid
sum(tdtall_f[11, ] != dtall_f[,37])

########## run #########3
DMRcate_wrapper(formula_pfdea, clinchem_f, tdtall_f)

DMRcate_wrapper(formula_pfna, clinchem_f, tdtall_f)

######### all male 
clinchem_m <- read.csv("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/clin_chem_m_nada.csv")
dim(clinchem_m)
# all clinchem
clinchem <- read.csv( "/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/clin_chem.csv",
                      header = T)
t_mval <- fread( "~/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/t_mval.csv",
                 header = T)
dt_all <- merge(clinchem, t_mval, by = "pid")

dim(dt_all)
dt_all <- data.table(dt_all)
dtall_m <-  dt_all[infant_sex == "Male"]
dim(dtall_m)

##
##
tdtall_m <- t(dtall_m[,-c(1:26)])
## check pid
sum(tdtall_m[11, ] != dtall_m[,37])
##
colnames(tdtall_m) <- clinchem_m$pid
########## run #########3
DMRcate_wrapper(formula_pfdea, clinchem_m, tdtall_m)

DMRcate_wrapper(formula_pfna, clinchem_m, tdtall_m)
