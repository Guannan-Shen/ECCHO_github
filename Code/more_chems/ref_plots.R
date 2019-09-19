############# diagnostic plots ###########
options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(tidyverse)
library(magrittr)
library(tools)
library(reshape2)
getwd()

## dot plots by chemicals ##
dot_groupby <- function(data, y, groupby, dotsize, xlab, ylab){
  p = ggplot(data, aes(x=groupby, y=y)) + 
    # dot plot
    geom_dotplot(binaxis='y', stackdir='center', dotsize = dotsize)+
    scale_colour_grey(start = 0.2, end = 0.8) +
    theme_bw() +
    # mean and 2 std
    stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max) +
    labs(x = xlab, y = ylab) +
    theme(axis.text.x = element_text(angle = 90))
  print(p)
}

## density plot for skewness ##
density_values_group <- function(data, title){
  print("A column called group and a column called values")
  p = ggplot(data, aes(x = values, color = ind, fill = ind)) + 
    # alpha controls the transparency 
    geom_density(alpha = 0.6) +
    # gray scale color
    # grey scale plot
    scale_colour_grey() +
    # start = 0.2, end = 0.8
    scale_fill_grey() +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal") +
    labs(caption = title) 
  print(p)
}

## heatmap of correlation 
cor_heatmap <- function(df, method, reorder, hclust_method, text_size){
  print("Use pairwise.complete.obs, and methods from pearson, kendall, spearman")
  cormat = cor(df, use = "pairwise.complete.obs", method = method)
  ###### Reordered correlation data visualization #####
  # using correlation as dist
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    # default settings of 
    hc <- hclust(dd, method = hclust_method)
    cormat <-cormat[hc$order, hc$order]
  }
  if(reorder){
    cormat <- reorder_cormat(cormat)
  }else{}
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] = NA
    return(cormat)
  }
  upper_tri = get_upper_tri(cormat)
  
  # reshape data 
  melted_cormat = melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name= method) +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed() +
    labs(x = "", y = "") + 
    geom_text(aes(Var2, Var1, label = round(value,2) ), color = "black", size = text_size)  +
    theme(
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  # Print the heatmap
  print(ggheatmap)
  return(cormat)
}



