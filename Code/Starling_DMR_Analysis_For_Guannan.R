#############################################################
#ON YOUR LOCAL MACHINE
# scp smiharry@lead.ucdenver.pvt:/home/smiharry/forGuannan/dmr.tar.gz
# ssh key passphrase: Be@trice9002
# NOTES: You will need to change the directories, and you will
# need to filter/use the reduced mval dataset.
#############################################################

library(DMRcate)
library(data.table)

## working dir
setwd("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/dmr/")
getwd()

## Read in clinical data (log2 transformed PFOA)
clindat <- fread(file = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/dmr/forGuannan/clindat.csv", 
                 header = T)
head(clindat)
dim(clindat)

## read in male female list
female <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/11_05_pfas_femalecpg.csv", header = T)
femaleid <- female$pid
length(femaleid)
male <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/11_05_pfas_malecpg.csv", header = T)
maleid <- male$pid
length(maleid)
# read in Filter CpGs based on data from Weiming(N = 433360)
filtered_cpg <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/MeanBeta_ExtremePFOA_M/MeanBeta_ExtremePFOA_M.txt",
                      header = F)
dim(filtered_cpg)
f_cpg <- filtered_cpg$V1 

############## build design matrix
formula <- ~as.numeric(log2pfoa) + Race
# formula
###### for female
clindat_f <- clindat[which(clindat$pid %in% femaleid) ,]
str(model <- model.frame(formula, clindat_f))
design_f <- model.matrix(formula, model)
dim(design_f)
##### for male 
clindat_m <- clindat[which(clindat$pid %in% maleid) ,]
str(model <- model.frame(formula, clindat_m))
design_m <- model.matrix(formula, model)
dim(design_m)
## read in methylation data (This should be cpgs as rows and subjects as columns)

## myMs <- logit2(myBetas)

mval <- fread(file = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/dmr/forGuannan/HS_450K_CB_Mval_normbatch_ForDMR-02-07-19.csv", header = T)

# cpg name as row name
rnames <- mval$V1
mval <- as.matrix(mval[, -1])
rownames(mval) <- rnames
## Filter subjects based on clinical data (N = 583)
mval <- mval[, which(colnames(mval) %in% clindat$pid)]
## Filter CpGs based on data from Weiming(N = 433360)
## GUANNAN: I DID NOT DO THIS HERE, SO DO NOT FORGET THIS STEP ##
mval <- mval[which(rownames(mval) %in% f_cpg ), ]
dim(mval)
mval_f <- mval[, which(colnames(mval) %in% femaleid )]
mval_m <- mval[, which(colnames(mval) %in% maleid )]

########### for female ##############
mval_f <- mval[, which(colnames(mval) %in% femaleid )]
mval_m <- mval[, which(colnames(mval) %in% maleid )]
## Build cpg anno object
cpganno <- cpg.annotate("array",
                        ## female 
                        mval_f,
                        what = "M",
                        arraytype = "450K",
                        analysis.type = "differential",
                        ## female
                        design = design_f,
                        ## only have the intercept and the pfoa concentration
                        ##  it will still be coef = 2. 
                        ## We are telling DMRcate that we are interested in the log transformed pfoa.
                        coef = 2)

# Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
# Your contrast returned no individually significant probes. Try increasing the fdr.
# Alternatively, set pcutoff manually in dmrcate() to return DMRs, 
# but be warned there is an increased risk of Type I errors.

dmrcoutput <- dmrcate(cpganno, 
                      lambda=1000, 
                      c=2, 
                      p.adjust.method="BH", 
                      consec=FALSE, 
                      pcutoff=0.05)

#######################################
# Get list of significant DMRs and CpGs
# associated with those DMRs
#######################################

dmr.results = dmrcoutput$results

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

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

##Build list of annotation data for probes in each DMR
DMRprobes = list()
for(i in 1:nrow(dmr.results)){
  DMRprobes[[i]] = getDMRprobes(dmr.results[i,], anno)
}

names(DMRprobes) <- dmr.results$coord
## explore dmr results
ls(dmrcoutput)

## find the DMR identifier
dmrcoutput$results$coord == names(DMRprobes)
length(names(DMRprobes))

## find FDR for individual cpg and merge data
dmr.results

## format name for naming the sheet
names(DMRprobes) <- gsub(":", "_", names(DMRprobes))



############################################
# Write results
############################################
library(tibble)
readme_sheet <- data_frame(
  Columns = c(
    "Annotations for CpGs assocaited with DMR's unadjusted p-value < 0.05",
    "", colnames(DMRprobes[[1]])
  ))

readme_sheet <- list(README = readme_sheet)
names(readme_sheet) <- "README"

openxlsx::write.xlsx(c(readme_sheet, 
                       DMRprobes),
                     "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/dmr/forGuannan/DMR_Results_female_loged_02-11-19.xlsx")


################## for male #####################
mval_m <- mval[, which(colnames(mval) %in% maleid )]
## Build cpg anno object
cpganno <- cpg.annotate("array",
                        mval_m,
                        what = "M",
                        arraytype = "450K",
                        analysis.type = "differential",
                        design = design_m,
                        ## only have the intercept and the pfoa concentration
                        coef = 2)

# Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
# Your contrast returned no individually significant probes. Try increasing the fdr.
# Alternatively, set pcutoff manually in dmrcate() to return DMRs, 
# but be warned there is an increased risk of Type I errors.

dmrcoutput <- dmrcate(cpganno, 
                      lambda=1000, 
                      c=2, 
                      p.adjust.method="BH", 
                      consec=FALSE, 
                      pcutoff=0.05)

#######################################
# Get list of significant DMRs and CpGs
# associated with those DMRs
#######################################

dmr.results = dmrcoutput$results

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

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

##Build list of annotation data for probes in each DMR
DMRprobes = list()
for(i in 1:nrow(dmr.results)){
  DMRprobes[[i]] = getDMRprobes(dmr.results[i,], anno)
}

names(DMRprobes) <- dmr.results$coord
names(DMRprobes) <- gsub(":", "_", names(DMRprobes))

############################################
# Write results
############################################
library(tibble)
readme_sheet <- data_frame(
  Columns = c(
    "Annotations for CpGs assocaited with DMR's unadjusted p-value < 0.05",
    "", colnames(DMRprobes[[1]])
  ))

readme_sheet <- list(README = readme_sheet)
names(readme_sheet) <- "README"

openxlsx::write.xlsx(c(readme_sheet, 
                       DMRprobes),
                     "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/dmr/forGuannan/DMR_Results_male_loged_02-11-19.xlsx")

## Your contrast returned no individually significant probes. 
## Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, 
## but be warned there is an increased risk of Type I errors.

######## transfer a variable name to a string ######
unlist(strsplit( deparse(substitute(clindat_f)), "_"))[2]
####### for this analysis Data ###########
## need 
c("clindat_f", "clindat_m", "which has log2pfoa, lnpfoa, Race",
  "mval_f", "mval_m",
  "formula_log2", "formula_ln")
## must put pfoa at the first place 
formula_log2 <-  ~as.numeric(log2pfoa) + Race
formula_ln <-  ~as.numeric(lnpfoa) + Race
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tibble)

###################################################
######################################################
####################################################################
################ function wrapper for DMR analysis ###################
DMRcate_wrapper <- function(formula, clindat, mval){
  ## dataset use , female or male, log2pfoa or lnpfoa ##
  gender1 = unlist(strsplit( deparse(substitute(clindat)), "_"))[2]
  gender2 = unlist(strsplit( deparse(substitute(mval)), "_"))[2]
  if(gender1 != gender2)
    stop("should use the same gender data")
  ## the type of pfoa
  log = unlist(strsplit( deparse(substitute(formula)), "_"))[2]
  prefix_name = paste(gender1,"_",log, "pfoa",sep = "")
  ### build design matrix ###
  str(model <- model.frame(formula, clindat))
  design = model.matrix(formula, model)
  
  ######## DMR analysis ########
  ## Build cpg anno object
  cpganno = cpg.annotate("array",
                          ## female 
                          mval,
                          what = "M",
                          arraytype = "450K",
                          analysis.type = "differential",
                          ## female
                          design = design,
                          ## only have the intercept and the pfoa concentration
                          ##  it will still be coef = 2. 
                          ## We are telling DMRcate that we are interested in the log transformed pfoa.
                          coef = 2)
  ## get results ##
  dmrcoutput = dmrcate(cpganno, 
                        lambda=1000, 
                        c=2, 
                        p.adjust.method="BH", 
                        consec=FALSE, 
                        pcutoff=0.05)
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
                       paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/dmr/",  prefix_name,
                             Sys.Date(),"_DMR", ".xlsx", sep= ""))
  
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
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/dmr/",  prefix_name,"_", 
                  Sys.Date(),"_DMR_allCpGs", ".csv", sep= ""))
  write.csv(cpg_min, 
            row.names = F,
            paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/dmr/",  prefix_name,"_", 
                  Sys.Date(),"_DMR_top1_CpG", ".csv", sep= ""))
}

DMRcate_wrapper(formula_log2, clindat_f, mval_f)
DMRcate_wrapper(formula_ln, clindat_f, mval_f)
DMRcate_wrapper(formula_log2, clindat_m, mval_m)
DMRcate_wrapper(formula_ln, clindat_m, mval_m)

library(data.table)
## fwrite from data.table 
fwrite(clindat_f, 
          row.names = F,
          paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/",  
                "clindat_f", ".csv", sep= ""))
fwrite(clindat_m, 
          row.names = F,
          paste("/home/guanshim/Documents/gitlab/ECCHO_github/DataProcessed/for_obesity/",  
                "clindat_m", ".csv", sep= ""))

