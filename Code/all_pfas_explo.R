rm(list = ls())
library(data.table)
library(knitr)
library(tidyverse)
library(magrittr)
library(stats)
library(DMRcate)
library(wesanderson)
library(outliers)
options(stringsAsFactors = F)
options(dplyr.width = Inf)
setwd("~/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/")
getwd()
## not in function
'%nin%' <- Negate('%in%')


# PFNA and PFDA data
pfna_da <- fread("pfas_methyl_di4.csv") %>% as.data.frame()
# data by gender
clinchem_f <- read.csv("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/clin_chem_f_nada.csv")
dim(clinchem_f)
clinchem_m <- read.csv("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/clin_chem_m_nada.csv")
dim(clinchem_m)
## pid, pfoa concentration and cpg M values, and clinical gender information 
pid583 <- fread("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/pid583.csv")
pid583 <- as.vector(unlist(pid583))
# check pid 
sum(pfna_da$pid != pid583)
new_chem <- c("PFDEA_ng_ml", "PFNA2_ng_ml")

## missingness check below the detection limit
apply(pfna_da[, new_chem], 2, function(x){
  # 0.1 is the detection limits 
  sum(x < 0.1)
})
apply(pfna_da[pfna_da$infant_sex == 1, new_chem], 2, function(x){
  # 0.1 is the detection limits 
  sum(x < 0.1)
})
apply(pfna_da[pfna_da$infant_sex == 2, new_chem], 2, function(x){
  # 0.1 is the detection limits 
  sum(x < 0.1)
})
pfna_da[pfna_da$infant_sex == 1, new_chem] %>% dim
########## more covariates ############
pfas_co <- fread("pfas_methyl_di5.csv") %>% as.data.frame()
# infant sex == 1, female, == 2, male
######### 
########### education 5 levels and 2 levels ########
bar_infantSex(pfas_co, pfas_co$epi5, "Education", legendname = "Maternal Education Level")
bar_infantSex(pfas_co, pfas_co$edu2, "Education_2", legendname = "Race")

######### Race 4 levels and 2 levels ########
bar_infantSex(pfas_co, pfas_co$race_4, "Race", legendname = "Race")
bar_infantSex(pfas_co, pfas_co$race2, "Race_2", legendname = "Race")

########## gestsmoking #############
bar_infantSex(pfas_co, pfas_co$gestsmoking, "Smoke", legendname = "Smoke")

#######33 previous pregnancies ##########
bar_infantSex(pfas_co, pfas_co$prev_preg, "Previous Pregnancy", legendname = "Previous Pregnancy")

############3 continuous ###########
#### maternal age #333
density_bygender(pfas_co, pfas_co$maternal_age, "Maternal Age")
density_bygender(pfas_co, pfas_co$pre_preg_bmi, "Maternal BMI")
density_bygender(pfas_co, log(pfas_co$pre_preg_bmi , base = exp(1)), "Maternal BMI (ln)")
grubbs.test(pfas_co$pre_preg_bmi, 10)

grubbs.test(log(pfas_co$pre_preg_bmi , base = exp(1)), 10)

########### generate gender straitified datasets #########
# get female data
pfas_f <- pfas_co %>% dplyr::filter(pid %in% pid583) %>% dplyr::filter(infant_sex == 1)
dim(pfas_f)
sum(pfas_f$pid != clinchem_f$pid)
clin_pfas_f <- clinchem_f %>% dplyr::mutate(edu = pfas_f$edu2 ,
                                            prev_preg = pfas_f$prev_preg,
                                            prev_preg_bmi_ln = log(pfas_f$pre_preg_bmi , base = exp(1)),
                                            smoke = pfas_f$gestsmoking)
# male 
pfas_m <- pfas_co %>% dplyr::filter(pid %in% pid583) %>% dplyr::filter(infant_sex == 2)
dim(pfas_m)
sum(pfas_m$pid != clinchem_m$pid)
clin_pfas_m <- clinchem_m %>% dplyr::mutate(edu = pfas_m$edu2 ,
                                            prev_preg = pfas_m$prev_preg,
                                            prev_preg_bmi_ln = log(pfas_m$pre_preg_bmi , base = exp(1)),
                                            smoke = pfas_m$gestsmoking)
################### save #############
## datasets with more covariates 
write.csv(clin_pfas_f, row.names = F,
          "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/clin_pfas_f_co.csv")
write.csv(clin_pfas_m, row.names = F,
          "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/clin_pfas_m_co.csv")

###################33 univariate test against chemical concentration and clinical outcomes #######3
#### obesity outcomes and chemicals
colnames(clin_pfas_f[c(5:13, 24:28)])
##### 6 covariates 
colnames(clin_pfas_f)[c(3:4, 29:32)]
#######3 univariate lm get p-value of beta1
get_uni_b<- function(x, y, df){
  a_lm = lm(y ~ x, data = df)
  sum_lm = summary(a_lm)
  beta = sum_lm$coefficients[2,1]
  return(beta )
}
get_uni_p <- function(x, y, df){
  a_lm = lm(y ~ x, data = df)
  sum_lm = summary(a_lm)
  p = sum_lm$coefficients[2,4]
  return(p )
}

######### test of birth weight and race 2 
a_lm <- lm(clin_pfas_f$birth_weight ~ clin_pfas_f$Race, data = clin_pfas_f)
sum_lm <- summary(a_lm)
sum_lm$coefficients[2,c(1,4)]

######### run #############
uni_test <- function(data, x, gender){
  ####33333 get coef ######3
  df1 = apply(data[, c(5:13, 24:28)], 2, function(y){
    round(get_uni_b(data[,x] , y, data), 4)
  })  %>% as.matrix()
  ####### get p #########
  df2 = apply(data[, c(5:13, 24:28)], 2, function(y){
    format.pval( get_uni_p(data[,x] , y, data), digits = 4)
  })  %>% as.matrix() 
  ######3 combine ##########
  df = cbind(df1, df2) %>% as.data.frame() %>% set_colnames(paste0(c("Coef. ", "p value "), x) ) 
  write.csv(df, row.names = T,
            paste0("/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/unitest/", 
                   x, gender, ".csv")
  )
  return(df)

}
########### test against 6 covariates #######
colnames(clin_pfas_f)[c(3:4, 29:32)]
for ( i in colnames(clin_pfas_f)[c(3:4, 29:32)]){
  uni_test(clin_pfas_f, i, "female")
}
##### male ##########
for ( i in colnames(clin_pfas_m)[c(3:4, 29:32)]){
  uni_test(clin_pfas_m, i, "male")
}

# colnames(clin_pfas_m)[c(3:4, 29:32)]
# colnames(clin_pfas_m)[c(24:28)]
dim(clin_pfas_f)
dim(clin_pfas_m)
# get_uni_p <- function(x, y, df){
#   a_lm = lm(y ~ x, data = df)
#   sum_lm = summary(a_lm)
#   beta = sum_lm$coefficients[2,1]
#   p = sum_lm$coefficients[2,4]
#   return(c(beta, p) )
# }