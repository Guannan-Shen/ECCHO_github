options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(tidyverse)
library(magrittr)
library(tools)
getwd()

#### summary table ####

sum_table <- function(df){
  sapply(colnames(df), function(x){
    c(summary(df[, x]), quantile(df[, x], 0.975))
  }) %>% as.data.frame() %>% t %>% as.data.frame()
}

