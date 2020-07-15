


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
library(matlib)
library(dplyr)
library(wCorr)
library(philentropy)


setwd('.....')

Disease_temp = read.csv(
  file = './InputFiles/Disease.csv',
  row.names = 1,
  header = TRUE,
  sep = ','
)
Disease <-
  Disease_temp[, 2:dim(Disease_temp)[2]] # keep the Disease samples



Healthy_temp = read.csv(
  file = './InputFiles/Healthy.csv',
  row.names = 1,
  header = TRUE,
  sep = ','
)
Healthy <-
  Healthy_temp[, 2:dim(Healthy_temp)[2]] # keep the Healthy samples




######################################################################################################################
######################################################################################################################
######################################################################################################################
#filter out genes with many zeros accross samples.

combinedDatasets <- cbind(Healthy, Disease)
keep <-
  which(rowSums(combinedDatasets == 0) <= round(0.8 * dim(combinedDatasets)[2]))
Healthy <- Healthy[keep, ]
Disease <- Disease[keep, ]


# Atlas Elevated genes
AtlasGenes_temp = read.csv(
  file = './InputFiles/HPAgenes.csv',
  row.names = NULL,
  header = TRUE,
  sep = ','
)




# Find the Atlas genes that exist in the Healthy dataset
kept_index1 <- c()
for (i in 1:length(rownames(Healthy))) {
  index <- which(AtlasGenes_temp$Ensembl == rownames(Healthy)[i])
  if (length(index) > 0) {
    kept_index1 <- cbind(kept_index1, index)
  }
}


# Find the Atlas genes that exist in the Disease dataset
kept_index2 <- c()
for (i in 1:length(rownames(Disease))) {
  index <- which(AtlasGenes_temp$Ensembl == rownames(Disease)[i])
  if (length(index) > 0) {
    kept_index2 <- cbind(kept_index2, index)
  }
}


#common atlas genes in Healthy and Disease datasets
common_indices <- intersect(kept_index1, kept_index2)

AtlasGenes = AtlasGenes_temp$Ensembl[common_indices]#Final tissue Atlas genes


# The two groups of tissues that we want to compare
group1 = Healthy #Reference tissue samples
group2 =  Disease # The tissue samples that we want to measure its distance from the reference (e.g. Healthy Tissue)



######################################################################################################################
######################################################################################################################
######################################################################################################################

# constructing the Atlas dataset for Healthy!
keeprow <- c()
for (i in 1:length(AtlasGenes)) {
  index <- which(rownames(group1) == AtlasGenes[i])
  if (length(index) > 0) {
    keeprow <- rbind(keeprow, index)
  }
}

group1_Atlas_DataSet <- group1[keeprow, ]
group1_Atlas_DataSet <-
  t(group1_Atlas_DataSet) # we transpose the matrix such as rows being the samples and columns the genes


# constructing the Atlas probability vector for Healthy!
Prob_group1_Atlas_DataSet <- group1_Atlas_DataSet

for (i in 1:dim(group1_Atlas_DataSet)[1]) {
  # for each tissue
  for (j in 1:dim(group1_Atlas_DataSet)[2]) {
    # for each gene
    Prob_group1_Atlas_DataSet[i, j] <-
      group1_Atlas_DataSet[i, j] / sum(group1_Atlas_DataSet[i, ])
    
  }
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
# constructing the Atlas dataset for Disease!
keeprow <- c()
for (i in 1:length(AtlasGenes)) {
  index <- which(rownames(group2) == AtlasGenes[i])
  if (length(index) > 0) {
    keeprow <- rbind(keeprow, index)
  }
}

group2_Atlas_DataSet <- group2[keeprow, ]
group2_Atlas_DataSet <-
  t(group2_Atlas_DataSet) # we transpose the matrix such as rows being the samples and columns the genes


# constructing the Atlas probability vector for Disease!
Prob_group2_Atlas_DataSet <- group2_Atlas_DataSet
for (i in 1:dim(group2_Atlas_DataSet)[1]) {
  # for each tissue
  for (j in 1:dim(group2_Atlas_DataSet)[2]) {
    # for each gene
    Prob_group2_Atlas_DataSet[i, j] <-
      group2_Atlas_DataSet[i, j] / sum(group2_Atlas_DataSet[i, ])
    
  }
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################


# constructing the Ranking matrix of the atlas genes for Healthy!
Rankings_Atlas_group1 <-
  matrix(
    data = 0,
    nrow = length(AtlasGenes),
    ncol = dim(group1)[2]
  ) # The uncertainties of the genes

for (i in 1:dim(group1)[2]) {
  # for each sample
  #we order the matrix of expressions in decreasing order
  keepRow <- c()
  group1_Rankings <- group1[order(group1[, i], decreasing = TRUE), ]
  for (j in 1:length(AtlasGenes)) {
    # for each atlas gene
    index <-
      which(rownames(group1_Rankings) == AtlasGenes[j])
    Rankings_Atlas_group1[j, i] <- index
  }
  
}



#For each tissue calculate the standarized expression metrix of the atlas genes
Rankings_Atlas_group1 <-
  t(Rankings_Atlas_group1) # transpose matrix such as the columns to correspond to Atlas genes




# constructing the Ranking matrix of the atlas genes for Disease!

Rankings_Atlas_group2 <-
  matrix(
    data = 0,
    nrow = length(AtlasGenes),
    ncol = dim(group2)[2]
  ) # The uncertainties of the genes

for (i in 1:dim(group2)[2]) {
  # for each sample
  #we order the matrix of expressions in decreasing order
  keepRow <- c()
  group2_Rankings <- group2[order(group2[, i], decreasing = TRUE), ]
  for (j in 1:length(AtlasGenes)) {
    # for each atlas gene
    index <-
      which(rownames(group2_Rankings) == AtlasGenes[j])
    Rankings_Atlas_group2[j, i] <- index
  }
  
}


#For each tissue calculate the standarized expression metrix of the atlas genes
Rankings_Atlas_group2 <-
  t(Rankings_Atlas_group2) # transpose matrix such as the columns to correspond to Atlas genes


# calculate pairwise the JSD distances between the tissue groups of sampels (healthy and Disease)
SRJSD <- c()
# calculate pairwise the RCD distances between the tissue groups of sampels (healthy and Disease)
RCD <- c() # Ranking correlation distance

for (i in 1:dim(Prob_group1_Atlas_DataSet)[1]) {
  # for each tissue in Healthy
  for (j in 1:dim(Prob_group2_Atlas_DataSet)[1]) {
    # for each tissue in Disease
    
    
    w1 = 0.5
    w2 = 0.5
    
    
    No_rescal = rbind(Prob_group1_Atlas_DataSet[i, ], Prob_group2_Atlas_DataSet[j, ])
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
        cor(Rankings_Atlas_group1[i, ], Rankings_Atlas_group2[j, ], method =  "pearson")
      ))
    RCD <- cbind(RCD, temp3)
    
    
  }
}

# Transcriptomic signature distance
TSDM <- (0.5 * SRJSD) + (0.5 * RCD)


write.csv(TSDM, file = "Healthy_vs_Disease_TSD.csv", row.names = FALSE)
write.csv(RCD, file = "Healthy_vs_Disease_RCD.csv", row.names = FALSE)
write.csv(SRJSD, file = "Healthy_vs_Disease_SRJSD.csv", row.names = FALSE)
