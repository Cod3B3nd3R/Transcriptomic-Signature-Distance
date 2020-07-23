
rm(list = ls())
# install.packages('Metrics')
# install.packages('ggplot2')
# install.packages('matrixStats')
# install.packages('glasso')
# install.packages('huge')
# install.packages('matlib')
# install.packages('dplyr')
# install.packages('wCorr')
# install.packages('philentropy')

library(Metrics)
library(ggplot2)
library(matrixStats)
library(glasso)
library(huge)
#library(matlib)
library(dplyr)
library(wCorr)
library(philentropy)


# add the working directory
setwd('....')


# load the reference tissue file
Reference_Tissue_temp = read.csv(
  file = './InputFiles/Reference_healthy.csv',
  row.names = 1,
  header = TRUE,
  sep = ','
)
Reference_Tissue <-
  data.frame(Reference_Tissue_temp[, 1:dim(Reference_Tissue_temp)[2]], row.names = rownames(Reference_Tissue_temp)) # keep the reference samples


# load the file with the tissue samples that will be compared to the reference samples
Compared_Tissue_temp = read.csv(
  file = './InputFiles/Compare_disease.csv',
  row.names = 1,
  header = TRUE,
  sep = ','
)
Compared_Tissue <-
  data.frame(Compared_Tissue_temp[, 1:dim(Compared_Tissue_temp)[2]], row.names = rownames(Compared_Tissue_temp)) # keep the compared samples




# Load the signature genes
AtlasGenes_temp = read.csv(
  file = './InputFiles/HPAgenes.csv',
  row.names = NULL,
  header = TRUE,
  sep = ','
)


######################################################################################################################
######################################################################################################################
######################################################################################################################
#filter out genes with zero expression across samples.

common_genes <-
  intersect(rownames(Reference_Tissue), rownames(Compared_Tissue))


keepRows1 <- c()
keepRows2 <- c()
for (i in 1:length(common_genes)) {
  keep1 <- which(rownames(Reference_Tissue) == common_genes[i])
  keepRows1 <- rbind(keepRows1, keep1)
  
  keep2 <- which(rownames(Compared_Tissue) == common_genes[i])
  keepRows2 <- rbind(keepRows2, keep2)
}


Reference_Tissue <- Reference_Tissue[keepRows1, , drop = FALSE]

Compared_Tissue <- Compared_Tissue[keepRows2, , drop = FALSE]



combinedDatasets <- cbind(Reference_Tissue, Compared_Tissue)

# keep only the genes that do not have zero expression across all samples
keep <-
  which(rowSums(combinedDatasets) != 0)
Reference_Tissue <- Reference_Tissue[keep, , drop = FALSE]
Compared_Tissue <- Compared_Tissue[keep, , drop = FALSE]
combinedDatasets <- combinedDatasets[keep, , drop = FALSE]

######################################################################################################################
######################################################################################################################
######################################################################################################################


# Find the Atlas signature genes in the Reference and Compared tissue
kept_index1 <- c()
for (i in 1:length(rownames(combinedDatasets))) {
  index <-
    which(AtlasGenes_temp$Ensembl == rownames(combinedDatasets)[i])
  if (length(index) > 0) {
    kept_index1 <- cbind(kept_index1, index)
  }
}

AtlasGenes = AtlasGenes_temp$Ensembl[kept_index1]#Final tissue Atlas signature genes



######################################################################################################################
######################################################################################################################
######################################################################################################################

# constructing the Atlas dataset
keeprow <- c()
for (i in 1:length(AtlasGenes)) {
  index <- which(rownames(combinedDatasets) == AtlasGenes[i])
  if (length(index) > 0) {
    keeprow <- rbind(keeprow, index)
  }
}


#Reference Tissue!
Reference_Tissue_Atlas_DataSet <-
  Reference_Tissue[keeprow, , drop = FALSE]
Reference_Tissue_Atlas_DataSet <-
  t(Reference_Tissue_Atlas_DataSet) # we transpose the matrix: the rows correspond to samples and columns to genes


#Compared Tissue
Compared_Tissue_Atlas_DataSet <-
  Compared_Tissue[keeprow, , drop = FALSE]

Compared_Tissue_Atlas_DataSet <-
  t(Compared_Tissue_Atlas_DataSet) # we transpose the matrix: the rows correspond to samples and columns to genes




# constructing the Atlas probability vector for Reference Tissue!
Prob_Reference_Tissue_Atlas_DataSet <-
  Reference_Tissue_Atlas_DataSet

for (i in 1:dim(Reference_Tissue_Atlas_DataSet)[1]) {
  # for each tissue
  for (j in 1:dim(Reference_Tissue_Atlas_DataSet)[2]) {
    # for each gene
    Prob_Reference_Tissue_Atlas_DataSet[i, j] <-
      Reference_Tissue_Atlas_DataSet[i, j] / sum(Reference_Tissue_Atlas_DataSet[i,])
    
  }
}



# constructing the Atlas probability vector for Compared Tissue!
Prob_Compared_Tissue_Atlas_DataSet <- Compared_Tissue_Atlas_DataSet
for (i in 1:dim(Compared_Tissue_Atlas_DataSet)[1]) {
  # for each tissue
  for (j in 1:dim(Compared_Tissue_Atlas_DataSet)[2]) {
    # for each gene
    Prob_Compared_Tissue_Atlas_DataSet[i, j] <-
      Compared_Tissue_Atlas_DataSet[i, j] / sum(Compared_Tissue_Atlas_DataSet[i,])
    
  }
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################


# constructing the Ranking matrix of the Atlas signature genes for the Reference Tissue!

#initialization
Rankings_Atlas_Reference_Tissue <-
  matrix(
    data = 0,
    nrow = length(AtlasGenes),
    ncol = dim(Reference_Tissue)[2]
  )


for (i in 1:dim(Reference_Tissue)[2]) {
  # for each sample
  #we order the matrix of expressions in decreasing order
  keepRow <- c()
  Reference_Tissue_Rankings <-
    Reference_Tissue[order(Reference_Tissue[, i], decreasing = TRUE), , drop = FALSE]
  for (j in 1:length(AtlasGenes)) {
    # for each atlas gene
    index <-
      which(rownames(Reference_Tissue_Rankings) == AtlasGenes[j])
    Rankings_Atlas_Reference_Tissue[j, i] <- index
  }
  
}


Rankings_Atlas_Reference_Tissue <-
  t(Rankings_Atlas_Reference_Tissue) # transpose: the columns correspond to Atlas genes




# constructing the Ranking matrix of the atlas signaturegenes for Compared Tissue!

Rankings_Atlas_Compared_Tissue <-
  matrix(
    data = 0,
    nrow = length(AtlasGenes),
    ncol = dim(Compared_Tissue)[2]
  )

for (i in 1:dim(Compared_Tissue)[2]) {
  # for each sample
  #we order the matrix of expressions in decreasing order
  keepRow <- c()
  Compared_Tissue_Rankings <-
    Compared_Tissue[order(Compared_Tissue[, i], decreasing = TRUE), , drop = FALSE]
  for (j in 1:length(AtlasGenes)) {
    # for each atlas gene
    index <-
      which(rownames(Compared_Tissue_Rankings) == AtlasGenes[j])
    Rankings_Atlas_Compared_Tissue[j, i] <- index
  }
  
}



Rankings_Atlas_Compared_Tissue <-
  t(Rankings_Atlas_Compared_Tissue) # transpose: the columns correspond to Atlas genes


# square root JSD distances between the tissue groups of samples (Reference and Compared)
SRJSD <- c() # square root Jenshen-Shannon Divergence
# RCD distances between the tissue groups of samples (Reference and Compared)
RCD <- c() # Ranking correlation distance

for (i in 1:dim(Prob_Reference_Tissue_Atlas_DataSet)[1]) {
  # for each tissue in Reference group
  for (j in 1:dim(Prob_Compared_Tissue_Atlas_DataSet)[1]) {
    # for each tissue in comparison
    
    
    w1 = 0.5
    w2 = 0.5
    
    
    No_rescal = rbind(Prob_Reference_Tissue_Atlas_DataSet[i,],
                      Prob_Compared_Tissue_Atlas_DataSet[j,])
    No_rescal_JSD <-
      gJSD(
        No_rescal,
        unit = "log2",
        weights = c(w1, w2),
        est.prob = NULL
      )
    
    
    SRJSD <- cbind(SRJSD, sqrt(No_rescal_JSD))
    
    
    #note: Ranking correlation distance
    temp3 <-
      sqrt(1 - max(
        0,
        cor(
          Rankings_Atlas_Reference_Tissue[i,],
          Rankings_Atlas_Compared_Tissue[j,],
          method =  "pearson"
        )
      ))
    RCD <- cbind(RCD, temp3)
    
    
  }
}

# Transcriptomic signature distance
TSD <- (0.5 * SRJSD) + (0.5 * RCD)


write.csv(TSD, file = "./Distance_Results/TSD.csv", row.names = FALSE)
write.csv(RCD, file = "./Distance_Results/RCD.csv", row.names = FALSE)
write.csv(SRJSD, file = "./Distance_Results/SRJSD.csv", row.names = FALSE)
