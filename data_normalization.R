
#------------------------------------------ #### 1) Function for normalization of epigenetic features and gene expression values #### -----------------------------------------------------

normalize.features <- function(df, method = 'log') {
  if (method == 'log') {
    
    normalized_features <- log.normalize.features(df)
    
  } else if (method == 'minmax') {
    
    normalized_features <- min.max.normalize.features(df)
    
  } else if (method == 'none') {
    
    normalized_features <- df
    
  } else {
    
    possible_methods <- c('log', 'minmax', 'none')
    stop(paste("Invalid method '", method, "'. Possible arguments are:", paste(possible_methods, collapse = ", ")))
  }
  
  return(normalized_features)
}

normalize.gex <- function(df, method = 'log') {
  if (method == 'log') {
    
    normalized_gex <- log.normalize.gex(df)
    
  } else if (method == 'minmax') {
    
    normalized_gex <- min.max.normalize.gex(df)
    
  } else if (method == 'none') {
    
    normalized_gex <- df
    
  } else {
    
    possible_methods <- c('log', 'minmax', 'none')
    stop(paste("Invalid method '", method, "'. Possible arguments are:", paste(possible_methods, collapse = ", ")))
  }
  
  return(normalized_gex)
}

# feature value for each bin of each gene per marker is normalized by the maximum value of the marker over all bins of all genes
min.max.normalize.features <- function(df){
  df.norm <- df
  maximums <- vector("numeric")
  for (i in 1:ncol(df)) {
    if (colnames(df[i]) %in% markers){
      column_max <- max(unlist(lapply(df[[i]], max)))
      maximums <- c(maximums, column_max)
      for(j in 1:nrow(df)){
        df.norm[j,i][[1]] <- lapply(df[j,i],function(x) x/column_max)
      }
    }
  }
  return(df.norm)
}

log.normalize.features <- function(df){
  df.norm <- df
  for (i in 1:ncol(df)) {
    if (colnames(df[i]) %in% markers){
      for(j in 1:nrow(df)){
        df.norm[j,i][[1]] <- lapply(df[j,i],function(x) log(x + 1e-5))
      }
    }
  }
  return(df.norm)
}


log.normalize.gex <- function(df){
  df.norm <- df
  df.norm$gex <- log(df.norm$gex + 1e-5)
  return(df.norm)
}

min.max.normalize.gex <- function(df){
  df.norm <- df
  max.gex <- max(df$gex)
  df.norm$gex <- df.norm$gex/max.gex
  return(df.norm)
}
