# Import Dataset
# load the data 
library(readr)
library(magrittr)
library(tidyverse)
# the data set 
dental <- read_csv("C:/Users/hithr/Documents/Stats/6624AdvStatsAnalysis/bios6624-Guannan-Shen/Project0/DataRaw/Project0_dental_data.csv")
#########################################
######### Data Preprocess
#########################################
# brief of the data 
str(dental)
# with missing value, create categorical variable and change
dental_change_t <- dental %>% mutate(trtgroup = factor(trtgroup, levels = c(2,1,3,4,5), 
                                              labels= c("no treatment", "placebo", "low", "medium", "high")))%>% 
  mutate(race = factor(race, levels= c(5,4,2,1), labels = c("White", "Asian", "African American", "Native American")))%>% 
  mutate(sex = factor(sex, levels = c(1,2), labels = c("Male", "Female"))) %>% 
  mutate(smoker = factor(smoker, levels = c(1,0), labels = c("Yes", "No"))) %>% 
  mutate(attach_change = attach1year - attachbase, pd_change = pd1year - pdbase)
str(dental_change_t)

dental_change_t <- data.frame(dental_change_t)
save(dental_change_t, file = "dental_change_t")
apply(dental_change_t, 2, is.numeric)
apply(dental, 2, is.numeric)
########################################
########## Univariate model loop of linear model
########################################
## standard
# test the categorical variable 
change_trt <- lm(attach_change ~ trtgroup, data = dental_change_t)
summary(change_trt)
# test the continuous variable 
change_age <- lm(attach_change ~ sex, data = dental_change_t)
summary(change_age)

####################
#########  apply lapply sapply
#####################
## lm apply
## in the apply, the lm(), the levels of categorical variable were forced to be alphabetical
# &&&&&& all the continuous variables were treated as categorical variables  
################################################################
uni_1 <- apply(dental_change_t[, c(2:8,10)], 2,
      function(col){
         lm = lm(dental_change_t$attach_change ~ col)
        coef = round(summary(lm)$coefficients[2,4],3)
        return(coef)
      })
uni_1
# lm sapply 
uni_sapply_lm <- sapply( c(2:8,10), function(col){
  lm = lm(dental_change_t$attach_change ~ dental_change_t[ , col])
  coef = round(summary(lm)$coefficients[2,4],3)
  return(coef)
})
names(uni_sapply_lm) <- colnames(dental_change_t)[ c(2:8,10)]
uni_sapply_lm

# lm lapply 
uni_lapply_lm <- unlist(lapply( c(2:8,10), function(col){
  lm = lm(dental_change_t$attach_change ~ dental_change_t[ , col])
  coef = round(summary(lm)$coefficients[2,4],3)
  return(coef)
}))
names(uni_lapply_lm) <- colnames(dental_change_t)[ c(2:8,10)]
uni_lapply_lm

#################################
# return a data.frame?  for table generation 
##############################
uni_lapply <- data.frame(lapply( c(2:8,10), function(col){
  lm = lm(dental_change_t$attach_change ~ dental_change_t[ , col])
  coef = round(summary(lm)$coefficients[2,],3)
  return(coef)
}))
uni_lapply
## lm apply
uni_2 <- apply(dental_change_t[, c(2:8,10)], 2,
               function(col){
                 lm = lm(dental_change_t$attach_change ~ col)
                 coef = round(summary(lm)$coefficients[2,],3)
                 return(coef)
               })
uni_2

# sapply 
# anova
uni_sapply <- sapply(c(2:4,7,8,10), function(col){
  base_lm = lm(dental_change_t$attach_change ~ 1 )
  lm = lm(dental_change_t$attach_change ~ dental_change_t[ , col])
  round(anova(base_lm,lm)$`Pr(>F)`[2],4)
})
names(uni_sapply) <- colnames(dental_change_t)[c(2:4,7,8,10)]
uni_sapply



###################################
## glm, 
####################################
uni_2 <- apply(dental_change_t[, c(2:8,10)], 2,
               function(col){
                 glm = glm(dental_change_t$attach_change ~ col, family = "gaussian")
                 coef = round(summary(glm)$coefficients[2,],3)
                 return(coef)
               })
uni_2
