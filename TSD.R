

rm(list=ls())
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
library('proxy') # Library of similarity/dissimilarity measures for 'dist()'


setwd('/Users/dimitris.manatakis/Desktop/quarantineDesktop/BioinformaticsRevision/Additional_Experiments_Revision/CompareToOtherDistances_Experiment2_AllGenes/LiverCancer/')


for (times in 1:3) {
  
  rm(list=setdiff(ls(), "times"))  
  
  if(times == 1){
    Cancer_Liver_Data_temp = read.csv(file = './InputFiles/HCC.csv', row.names = 1, header = TRUE,sep = ',')
    
  }
  
  else if (times ==2) {
    Cancer_Liver_Data_temp = read.csv(file = './InputFiles/CHC.csv', row.names = 1, header = TRUE,sep = ',')
    
    
  }  
  
  else if(times == 3){
    Cancer_Liver_Data_temp = read.csv(file = './InputFiles/CC.csv', row.names = 1, header = TRUE,sep = ',')
    
  }
  

  
  
  Control_temp = read.csv(file = './InputFiles/Control.csv', row.names = 1, header = TRUE,sep = ',')
  Control <- Control_temp[,2:dim(Control_temp)[2]] # keep the control samples
  
  
  Cancer_Liver_Data <- Cancer_Liver_Data_temp[,2:dim(Cancer_Liver_Data_temp)[2]] # keep the cancer liver samples
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  #filter out genes with many zeros accross samples.
  
  combinedDatasets <- cbind(Control,Cancer_Liver_Data)
  keep <- which(rowSums(combinedDatasets==0) <= round(0.5*dim(combinedDatasets)[2]))
  Control<-Control[keep,]
  Cancer_Liver_Data<-Cancer_Liver_Data[keep,]

  
  # The two groups of tissues that we want to compare
  Prob_group1_Atlas_DataSet = t(Control) #Reference tissue
  Atlas_Gold_Standard_Genes = rownames(Cancer_Liver_Data)#Reference tissue Atlas genes
  Prob_group2_Atlas_DataSet =  t(Cancer_Liver_Data) # The tissue 
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  # # constructing the Ranking matrix of the atlas genes for group1!
  # 
  # Rankings_Atlas_group1 <- matrix(data = 0, nrow = length(T1_atlas_gene_names), ncol =dim(group1)[2]) # The uncertainties of the genes
  # 
  # for(i in 1 : dim(group1)[2]){ # for each sample
  #   #we order the matrix of expressions in decreasing order
  #   keepRow<-c()
  #   group1_Rankings <- group1[order(group1[,i],decreasing = TRUE),]
  #   for(j in 1 : length(T1_atlas_gene_names)){ # for each atlas gene
  #     index <- which(rownames(group1_Rankings) == T1_atlas_gene_names[j])
  #     Rankings_Atlas_group1[j,i] <- index  
  #   }
  #   
  # }
  # 
  # 
  # 
  # #step 1 for each tissue calculate the standarized expression metrix of the atlas genes 
  # Rankings_Atlas_group1 <- t(Rankings_Atlas_group1) # transpose matrix such as the columns to correspond to Atlas genes
  # 
  # 
  # 
  # 
  # # constructing the Ranking matrix of the atlas genes for group2!
  # 
  # Rankings_Atlas_group2 <- matrix(data = 0, nrow = length(T2_atlas_gene_names), ncol =dim(group2)[2]) # The uncertainties of the genes
  # 
  # for(i in 1 : dim(group2)[2]){ # for each sample
  #   #we order the matrix of expressions in decreasing order
  #   keepRow<-c()
  #   group2_Rankings <- group2[order(group2[,i],decreasing = TRUE),]
  #   for(j in 1 : length(T2_atlas_gene_names)){ # for each atlas gene
  #     index <- which(rownames(group2_Rankings) == T2_atlas_gene_names[j])
  #     Rankings_Atlas_group2[j,i] <- index  
  #   }
  #   
  # }
  # 
  # 
  # #step 1 for each tissue calculate the standarized expression metrix of the atlas genes 
  # Rankings_Atlas_group2 <- t(Rankings_Atlas_group2) # transpose matrix such as the columns to correspond to Atlas genes
  # 
  
  # calculate pairwise the JSD distances between the tissue groups
  euclid <-c()
  # calculate pairwise the RCD distances between the tissue groups
  manhat <-c() # Ranking correlation distance
  
  Pearson <- c()
  Kendall <- c()
  cosine <- c()
  
  for(i in 1: dim(Prob_group1_Atlas_DataSet)[1]){ # for each tissue in group1
    for(j in 1: dim(Prob_group2_Atlas_DataSet)[1]){ # for each tissue in group2
   
      euclid <- cbind(euclid,dist(rbind(Prob_group1_Atlas_DataSet[i,],Prob_group2_Atlas_DataSet[j,]),method = 'euclidean'))
      manhat <- cbind(manhat,dist(rbind(Prob_group1_Atlas_DataSet[i,],Prob_group2_Atlas_DataSet[j,]),method = 'manhattan' ))
      Pearson <- cbind(Pearson,cor(Prob_group1_Atlas_DataSet[i,], Prob_group2_Atlas_DataSet[j,],  method = "pearson"))
      Kendall <- cbind(Kendall, cor(Prob_group1_Atlas_DataSet[i,], Prob_group2_Atlas_DataSet[j,],  method = "kendall"))
      cosine <- cbind(cosine, dist(rbind(Prob_group1_Atlas_DataSet[i,], Prob_group2_Atlas_DataSet[j,]),  method = "cosine"))
      
      
    }
  }
  
  # Transcriptomic signature distance

  
  if(times == 1){
    write.csv(euclid, file = "Control_vs_HCC_Euclid.csv",row.names=FALSE)
    write.csv(manhat, file = "Control_vs_HCC_Manhat.csv",row.names=FALSE)
    write.csv(Pearson, file = "Control_vs_HCC_Pearson.csv",row.names=FALSE)
    write.csv(Kendall, file = "Control_vs_HCC_Kendall.csv",row.names=FALSE)
    write.csv(cosine, file = "Control_vs_HCC_Cosine.csv",row.names=FALSE)
    
  }
  
  else if (times == 2) {
    write.csv(euclid, file = "Control_vs_CHC_Euclid.csv",row.names=FALSE)
    write.csv(manhat, file = "Control_vs_CHC_Manhat.csv",row.names=FALSE)
    write.csv(Pearson, file = "Control_vs_CHC_Pearson.csv",row.names=FALSE)
    write.csv(Kendall, file = "Control_vs_CHC_Kendall.csv",row.names=FALSE)
    write.csv(cosine, file = "Control_vs_CHC_Cosine.csv",row.names=FALSE)
    

  }  
  
  else if (times == 3) {
    write.csv(euclid, file = "Control_vs_CC_Euclid.csv",row.names=FALSE)
    write.csv(manhat, file = "Control_vs_CC_Manhat.csv",row.names=FALSE)
    write.csv(Pearson, file = "Control_vs_CC_Pearson.csv",row.names=FALSE)
    write.csv(Kendall, file = "Control_vs_CC_Kendall.csv",row.names=FALSE)
    write.csv(cosine, file = "Control_vs_CC_Cosine.csv",row.names=FALSE)
    
  }

}



