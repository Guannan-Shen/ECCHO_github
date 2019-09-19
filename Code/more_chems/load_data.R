# clear the memory
rm(list = ls())
library(tidyverse)
library(magrittr)

options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
# small functions, %nin%, %||%, isnothing
source( paste0(dir, "Code/small_fun.R") )

######### load the more chemicals data ############

dir = "~/Documents/gitlab/ECCHO_github/"
phthal_di2<- read.csv(paste0(dir,"DataRaw/more_chems/phthal_di2_aug2019.csv"), 
                            header=TRUE)

dim(phthal_di2)
colnames(phthal_di2)
