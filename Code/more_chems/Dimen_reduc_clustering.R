## rescale data 
rescale0_1 <- function(df){
  data = apply(df, 2, scale)
  return(data)
}


## correlation distance ##
dist_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd = as.dist((1-cormat)/2)
  return(dd)
}

## Hierarchical Clustering
my_h_clust <- function(dist, method, title){
 hc = hclust(dist, method = method)
 plot(hc,
      main = title)
 return(hc)
}

######### PCA analysis by prcomp is preferred #######
### prcomp is preferred for numerical accuracy #########

## 'princomp' can only be used with more units than variables
wrapper_prcomp <- function(data){
  print("Using prcomp, standard observation, variable data frame")
  # observation id 
  pid = rownames(data)
  # variable names
  var_names = colnames(data)
  pca = prcomp(data, center = TRUE, scale. = TRUE, retx = TRUE)
  # pc names 
  pc = colnames(pca$x)
  #### make variance by pc plot ####
  sd = ( pca$sdev/sum(pca$sdev) )* 100
  # ggplot2 to check the % of variance pca analysis 
  # Horizontal bar plot
  pc.sd = data.frame(PCs = pc,
                    Variance = sd)
  
  p_sd = ggplot(data=pc.sd, aes(x=1:length(pc), y=Variance)) +
    geom_bar(stat="identity", fill = "Black") +
    ylab("% Variance") +
    scale_x_continuous(name = "Principal components", breaks = 1:length(pc),
                       labels = pc)  +
    # Horizontal bar plot
    coord_flip() +
    theme_bw()
  print(p_sd)
  #### make pc1 vs pc2 plot ####
  pca12 = ggplot(data = data.frame(pca$x), aes(x = PC1, y = PC2)) +
    geom_point(shape = 19) +
    xlab(paste("PC1 (", round(pc.sd$Variance[1],2) ,"%", ")", 
               sep = "") ) +
    ylab(paste("PC2 (", round(pc.sd$Variance[2],2) ,"%", ")", 
               sep = "") ) +
    theme_bw()
  print(pca12)
  
  return(list(pc_byfeature = pca$rotation, pc_var = pca$sdev, obs_bypc = pca$x))
}


## $rotation is the decomposition of ecah pc components at the feature level

## kmeans clustering 
elbow_kmeans <- function(data, maxk, iter.max){
  k_c = 1:maxk
  k_sws = NULL ; k_wsm = NULL
  for(i in k_c){
    km <- kmeans(data, i, iter.max = iter.max)
    k_sws[i] <- sum(km$withinss)
    if(i == 1){
      k_wsm <- c(km$withinss)
    }else{
      k_wsm <- c(k_wsm, km$withinss)
    }
  }
  plot(k_c, k_sws, type = "o", 
       xlab = "Number of Clusters", ylab = "Sum of Within cluster sum of squares")
}


