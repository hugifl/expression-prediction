# expression-prediction

## Project description
The goal of this project was to predict gene expression from epigenetic features of the DNA, represented by multiple histone mark ChIP-seq datasets and one chromatin accessibility DNase-seq dataset.

## Data

The epigenetic data was provided as bigWig enrichment/read coverage files over the whole genome for three cell types X1/2/3.
The gene expression data was provided as CPM-normalized counts from two cell types X1/X2. 
For all genes, the chromosome, the location of the gene and its transcription start site on the chromosome and the direction of the gene on the chromosome were provided.

X1 and X2 were used to train and validate the model and X3 was used to test the model.


## How to use
The whole project can be run from the 'expression_prediction.R' script.
The bigWig files of the epigenetic tracks are not provided due to their size. The R workspace 'allfeaturesloaded.RData' contains the pre-loaded epigenetic features (7 bins, windowsize of 10'000) as well as the pre-loaded gene expression data (corresponding to the result of running all data-loading and -processing steps in 'expression_prediction.R' up to and including step 5 ( "5) Epigenetic feature extraction")).
