## set up workspace
rm(list = ls())
library(data.table)
library(tidyverse)
library(magrittr)
library(stats)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tibble)
library(grid) # low-level grid functions are required
options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
## not in function
'%nin%' <- Negate('%in%')


######### read and sort ###########
dir2 <- "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/combp_DMR/"