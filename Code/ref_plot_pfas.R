############# diagnostic plots ###########
## histogram and boxplot with stat_summary ##
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(reshape2)
library(readxl)
library(tidyverse)
library(magrittr)
library(tools)
library(extrafont)
library(wesanderson)
setwd("~/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/")
getwd()


############# barplot ############

bar_infantSex <- function(pfas_co, fill, title, legendname){
  ########### by sex #############
  p1 = ggplot(pfas_co, aes(x = as.character(infant_sex), 
                      fill = as.character(fill))) +
    geom_bar(position="dodge", stat="count") +
    scale_x_discrete( labels = c("1" = "Female","2" = "Male"  ), name = "Infant Sex" ) + 
    scale_y_continuous( name = "Count" ) +
    scale_fill_grey(name = legendname) +
    # scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1"), name = "Maternal Education Level") +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal" ) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial")) 
  print(p1)
  ggsave(filename =  paste0( title, "_by_gender.tiff"), device = NULL,
         path = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/",
         dpi = 300, compression = "lzw")
  ################# all ###############
  p2 = ggplot(pfas_co, aes(x = "All", 
                      fill = as.character(fill) )) +
    geom_bar(position="dodge", stat="count")  +
    scale_x_discrete(  name = "" ) + 
    scale_y_continuous( name = "Count" ) +
    scale_fill_grey(name = legendname) +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal" ) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial")) 
  print(p2)
  ggsave(filename =  paste0( title, "_all.tiff"), device = NULL,
         path = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/",
         dpi = 300, compression = "lzw")
}

density_bygender <- function(data, x, title){
  p1 = ggplot(data, aes(x = x)) +
    geom_dotplot(method="histodot",  dotsize = 0.25) +
    scale_x_continuous(  name = paste(title, "All") ) + 
    scale_y_continuous( name = "Count" )  +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial")) 
  print(p1)
  ggsave(filename =  paste0( title, "_all.tiff"), device = NULL,
         path = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/",
         dpi = 300, compression = "lzw")
  
  p2 = ggplot(data, aes(x = x, fill = as.character(infant_sex)) ) +
    geom_density(alpha = 0.5 )  + 
    scale_x_continuous(  name = paste(title, "by Gender") ) + 
    scale_y_continuous( name = "Density" )  +
    scale_fill_discrete(name = "Gender", labels = c("1" = "Female", "2" = "Male")) +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal" ) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial")) 
  print(p2)
  ggsave(filename =  paste0( title, "_by_gender.tiff"), device = NULL,
         path = "/home/guanshim/Documents/gitlab/ECCHO_github/DataRaw/more_pfas/",
         dpi = 300, compression = "lzw")
  
  
}


