#------------------------------------------ #### 1) Functions for loading gene info and gene expression data for the train, test and validation data  #### -----------------------------------------------------

# Cell lines X1 and X2 are used for training and validation. The same genes (on the same chromosomes) from X1 and X2 are used for training and the same genes for validation.
# Cell line X3 is used for testing.
# The final training data set will be a data frame with the gene names, the gene info (TSS start site, chromosome etc.) and the gene expression values. 
# The validation set is structured in the same way as the train set, as is the test set (without the gene expression column).
# Each gene in the validation and the train set is present as a duplicate (once from X1, once from X2).

gene.loader <- function(path, data.set.names){
  train.info.df <- data.frame()
  train.y.df <- data.frame()
  val.info.df <- data.frame()
  val.y.df <- data.frame()
  test.info.df <- data.frame()
  
  for (data.set.name in data.set.names){
    file_list <- list.files(path = path, pattern = paste0("^",data.set.name))
    for (file in file_list){
      my_parts <- str_split(file, pattern = "_")
      role <- my_parts[[1]][2]
      information <- my_parts[[1]][3]
      if (role == "train"){
        if (information == "y.tsv"){
          df <- read_tsv(paste0(path,file))
          df$dataset <- data.set.name
          train.y.df <- rbind(train.y.df, df)
        } else if (information == "info.tsv"){
          df <- read_tsv(paste0(path,file))
          df$dataset <- data.set.name
          train.info.df <- rbind(train.info.df, df)
        }
        
      } else if (role == "val"){
        if (information == "y.tsv"){
          df <- read_tsv(paste0(path,file))
          df$dataset <- data.set.name
          val.y.df <- rbind(val.y.df, df)
        } else if (information == "info.tsv") {         
          df <- read_tsv(paste0(path,file))
          df$dataset <- data.set.name
          val.info.df <- rbind(val.info.df, df)
        }
      } else if (role == "test"){
        df <- read_tsv(paste0(path,file))
        df$dataset <- data.set.name
        test.info.df <- rbind(test.info.df, df)
      }
    }
  }
  train.df <- merge(x=train.info.df,y=train.y.df, by=c("gene_name","dataset"))
  val.df <- merge(x=val.info.df,y=val.y.df, by=c("gene_name","dataset"))
  test.df <- test.info.df
  
  return(list(train.df, val.df, test.df))
}

#------------------------------------------ #### 2) Epigenetic features are extracted from the bigwig files #### -----------------------------------------------------

# To keep the correct orientation of the genetic information encoded in the genes, the bins for genes on the '-' strand are reversed.

epigenetic.marker.feature.extraction <- function(markers.path, markers, data.frame, n.bins, window.size){
  feature.df <- data.frame(data.frame$gene_name, data.frame$dataset, data.frame$chr, data.frame$strand)
  feature.df[markers] <- 0
  for (marker in markers){
    print(paste0(marker,": starting extraction"))
    for (i in c(1:nrow(data.frame))){
      if (i%%1000 == 0){
        print(paste0(i," genes extracted from marker: ", marker))
      }
      feature.df[[marker]][i] <- I(list(rep(0, n.bins)))
      path <- paste0(markers.path, marker, "-bigwig/", as.character(data.frame$dataset[i]),".bigwig")
      if (!file.exists(path)) {
        path <- paste0(markers.path, marker, "-bigwig/", as.character(data.frame$dataset[i]),".bw")
      }
      
      TSS <- (data.frame$TSS_start[i] + data.frame$TSS_end[i])/2
      start.bins <- TSS - window.size/2
      end.bins <- TSS + window.size/2
      bin.borders <- round(seq(from = start.bins, to = end.bins, length.out = (n.bins + 1)), digits = 0)
      
      bin.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("chr", "start", "end"))
      for (j in c(1:(length(bin.borders) - 1))){
        bin.df <- rbind(bin.df, data.frame(chr = data.frame$chr[i], start = bin.borders[j], end = bin.borders[j + 1]))
      }
      
      tryCatch({
        features <- PopSV::bin.bw(bw.file = path, bin.df = bin.df, outfile.prefix = NULL, appendIndex.outfile = FALSE,
                                  chunk.size = n.bins, check.chr.name = TRUE, no.checks = FALSE,
                                  read.length = 100, fromSummaries = TRUE)
        feature.df[[marker]][i] <- I(list(features$bc$bc))
        if (feature.df$data.frame.strand[i] == "-"){
          feature.df[[marker]][[i]] <- rev(feature.df[[marker]][[i]])
        }
      }, error = function(e) {
        print(paste0("no reads found from: ", marker," in gene:", as.character(data.frame$gene_name[i])))
        #feature.df[[marker]][i] <- I(list(features$bc$bc)) #I(list(rep(0, n.bins)))
      })
      
    }
  }
  return(feature.df)
}


flatten.bins <- function(bin.df, markers, n.bins){
  for (marker in markers){
    print(marker)
    for (i in c(1:n.bins)){
      new.col <- paste0(marker,"_bin_",i)
      bin.df[[new.col]] <- lapply(bin.df[[marker]], `[[`, i)
    }
  }
  return(bin.df)
}
