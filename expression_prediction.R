###############################################################################################################################\
##
# This program was originally written for the project work of the course 'Machine Learning for Genomics' at ETH Zurich.
# The approach was developed together with my buddy Paul Doucet and coded in Python originally.
# As an exercise for myself, I optimized and streamlined parts of the pipeline and implemented it in R.
##
###############################################################################################################################\






#------------------------------------------ #### 0) Loading packages #### -----------------------------------

library(tidyverse)
library(stringr)
library(PopSV)
library(ranger)
library(Hmisc)
library(latticeExtra)

#------------------------------------------ #### 1) Loading functions #### -----------------------------------------------------

source('data_loading.R')
source('data_normalization.R')

#------------------------------------------ #### 2) Path #### ------------------------------------------------------------------

file.path <- "./CAGE-train/"
data.sets <- c('X1', 'X2', 'X3')

#------------------------------------------ #### 3) Loading gene expression data #### ------------------------------------------

dfs <- gene.loader(path = file.path, data.set.names = data.sets)

X.train <- dfs[[1]]
X.val <- dfs[[2]]
X.test <- dfs[[3]]

#------------------------------------------ #### 4) Feature and Model specifications #### --------------------------------------

markers <- c('H3K4me3',"H3K4me1",'H3K36me3','DNase','H3K27me3') #epigenetic features used for prediction

window      <- 10000                     # size of the window around the TSS of each gene in which epigenetic features are used
bin.number  <- 7                         # number of bins, the window around the TSS of each gene will be partitioned into

# specify hyperparameter grid for random forest regressor
rf.min.trees <- 50                       # minimal number of trees of the random forest regressor that should be tried during model fitting
rf.max.trees <- 400                      # maximal number of trees of the random forest regressor that should be tried during model fitting
rf.min.depth <- 0                        # minimal depth of the trees used in the random forest regressor (set to 0 for unlimited depth)
rf.max.depth <- 5                        # maximal depth of the trees used in the random forest regressor

#------------------------------------------ #### 5) Epigenetic feature extraction #### -----------------------------------------

train.features <- epigenetic.marker.feature.extraction(markers.path = "./", markers = markers, 
                                                       data.frame = X.train, n.bins = bin.number, window.size = window)
val.features <- epigenetic.marker.feature.extraction(markers.path = "./", markers = markers, 
                                                       data.frame = X.val, n.bins = bin.number, window.size = window)
test.features <- epigenetic.marker.feature.extraction(markers.path = "./", markers = markers, 
                                                       data.frame = X.test, n.bins = bin.number, window.size = window)

#------------------------------------------ #### 6) Epigenetic feature and gene expression normalization #### ------------------

# Gene expression values have already been normalized by library size
# normalize.gex and normalize.features normalize gene expression and epigenetic feature read counts. 
# method-argument: 'log' (log normalization), 'minmax' (min-max normalization), 'none' (no further normalization)

X.train.norm <- normalize.gex(X.train, method = 'log')
X.val.norm <- normalize.gex(X.val, method = 'log')

train.features.norm <- normalize.features(train.features, method = 'log')
test.features.norm <- normalize.features(test.features, method = 'log')
val.features.norm <- normalize.features(val.features, method = 'log')

#------------------------------------------ #### 7) Epigenetic feature reshaping #### ------------------------------------------

# The epigenetic features for each marker and each gene are stored in a data frame, each entry being a vector of read counts per bin.
# The data structure with shape (#genes x #markers x #bins) is flattened to 2D (#genes x (#markers x #bins)), the columns being the separate bins for each marker.

train.features.flattened <- flatten.bins(train.features.norm, markers = markers, bin.number)
test.features.flattened <- flatten.bins(test.features.norm, markers = markers, bin.number)
val.features.flattened <- flatten.bins(val.features.norm, markers = markers, bin.number)

y_train <- X.train.norm$gex
y_val <- X.val.norm$gex
X_train <- train.features.flattened[,(ncol(train.features) + 1):ncol(train.features.flattened)]
X_val <- val.features.flattened[,(ncol(val.features) + 1):ncol(val.features.flattened)]
X_test <- test.features.flattened[,(ncol(test.features) + 1):ncol(test.features.flattened)]

#------------------------------------------ #### 8) Random forest regressor for gene expression prediction #### -----------------

# Grid search over the hyperparameter ranges specified in 4. 
# For the final model, the combination of hyperparamters will be used that reaches the highest average Spearman correlation during testing.
# Each set of hyperparameters is tested by 10 - fold cross validation on the training data (train and test set drawn from training data in ordered fashion).
# The fold indices are not randomly drawn to ensure that no gene is present in both the train and the test set.
#(each gene is in the training twice, once for cell line X1 and once for cell line X2).

# The validation data is only used to test the performance of the final model which is fitted again on the whole training data set.

trees.vec <- seq(from = rf.min.trees, to = rf.max.trees, by = as.integer((rf.max.trees - rf.min.trees)/10))
depth.vec <- c(rf.min.depth:rf.max.depth)

correlations <- expand.grid(trees = trees.vec, depth = depth.vec)
correlations$corr <- 0

folds <- cut(seq(1,nrow(X_train)),breaks=10,labels=FALSE)

counter <- 1
for (trees in trees.vec){
  for(depth in depth.vec){
    correlation <- rep(0,10)
    for(i in 1:10){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      
      testData <- X_train[testIndexes, ]
      trainData <- X_train[-testIndexes, ]
      testFeatures <- y_train[testIndexes]
      trainFeatures <- y_train[-testIndexes]
    
      clf <- ranger(trainFeatures ~ ., data = trainData, num.trees = trees, max.depth = depth)
      pred <- predict(clf, data = testData)$predictions
      correlation[i] <- rcorr(testFeatures, pred, type = "spearman")$r[2]
    }
    avg.correlation <- mean(correlation)
    row_index <- which(correlations$trees == trees & correlations$depth == depth)
    correlations[row_index,]$corr <- avg.correlation
  }
  
  print(counter/length(trees.vec))
  counter <- counter + 1
}


n_trees <- correlations[which.max(correlations$corr),]$trees
max.depth <- correlations[which.max(correlations$corr),]$depth

# final model for validation
clf <- ranger(y_train ~ ., data = X_train, num.trees = n_trees, max.depth = max.depth)
pred <- predict(clf, data = X_val)$predictions
corr <- rcorr(y_val, pred, type = "spearman")

corr$r[2]

# prediction of test set
pred.test <- predict(clf, data = X_test)$predictions

X.test$gex_pred <- pred.test

write.csv(X.test, "Test_gene_expression_predicted.csv", row.names=FALSE)

