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
clindat <- fread(file = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/dmr/forGuannan/pid583_pfoa_log.csv", 
                 header = T)
head(clindat)

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
formula <- ~as.numeric(log2pfoa)
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

########### for female ##############
mval_f <- mval[, which(colnames(mval) %in% femaleid )]
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






