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


for (counter1 in 14:16){
  
  rm(list= ls()[!(ls() %in% c("counter1"))])  
  
  library(Metrics)
  library(ggplot2)
  library(matrixStats)
  library(glasso)
  library(huge)
  library(matlib)
  library(wCorr)
  library(philentropy)
  
  
  
  setwd('C:/Users/dimitris.manatakis/Desktop/3.1Evaluating_wTSD_UsingDataFromDifferentOrganTissues/TSD_forThe13_OrganTissues/Bladder')
  
  
  # Loading Datasets
  
  Liver_Data = read.csv(file = '../../DatasetsofOrganTissues/Liver.csv', row.names = 1, header = TRUE,sep = ',')
  
  Kidney_Data = read.csv(file = '../../DatasetsofOrganTissues/Kidney.csv', row.names = 1, header = TRUE,sep = ',')
  
  Small_Intestine_Data = read.csv(file = '../../DatasetsofOrganTissues/Small_intestine.csv', row.names = 1, header = TRUE,sep = ',')
  
  Stomach_Data = read.csv(file = '../../DatasetsofOrganTissues/Stomach.csv', row.names = 1, header = TRUE,sep = ',')
  
  Skin_Data = read.csv(file = '../../DatasetsofOrganTissues/Skin.csv', row.names = 1, header = TRUE,sep = ',')
  
  Brain_Data = read.csv(file = '../../DatasetsofOrganTissues/Brain.csv', row.names = 1, header = TRUE,sep = ',')
  
  Pancreas_Data = read.csv(file = '../../DatasetsofOrganTissues/Pancreas.csv', row.names = 1, header = TRUE,sep = ',')
  
  Lung_Data = read.csv(file = '../../DatasetsofOrganTissues/Lung.csv', row.names = 1, header = TRUE,sep = ',')
  
  Ovary_Data = read.csv(file = '../../DatasetsofOrganTissues/Ovary.csv', row.names = 1, header = TRUE,sep = ',')
  
  Prostate_Data = read.csv(file = '../../DatasetsofOrganTissues/Prostate.csv', row.names = 1, header = TRUE,sep = ',')
  
  Bladder_Data = read.csv(file = '../../DatasetsofOrganTissues/Bladder.csv', row.names = 1, header = TRUE,sep = ',')
  
  Esophagus_Data = read.csv(file = '../../DatasetsofOrganTissues/Esophagus.csv', row.names = 1, header = TRUE,sep = ',')
  
  Thyroid_Data = read.csv(file = '../../DatasetsofOrganTissues/Thyroid.csv', row.names = 1, header = TRUE,sep = ',')
  
  Adrenal_Data = read.csv(file = '../../DatasetsofOrganTissues/Adrenal_gland.csv', row.names = 1, header = TRUE,sep = ',')
  
  Cervix_Data = read.csv(file = '../../DatasetsofOrganTissues/Cervix.csv', row.names = 1, header = TRUE,sep = ',')
  
  Skeletal_Data = read.csv(file = '../../DatasetsofOrganTissues/Skeletal_muscle.csv', row.names = 1, header = TRUE,sep = ',')
  
  
  
  # Atlas Elevated genes
  
  Liver_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Liver_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Kidney_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Kidney_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Intestine_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/SmallIntestine_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Stomach_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Stomach_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Skin_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Skin_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Brain_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Brain_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Pancreas_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Pancreas_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Lung_AtlasGenes = read.csv(file = '../../HPA_SignatureGenes/Lung_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Ovary_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Ovary_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Prostate_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Prostate_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Bladder_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Bladder_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Esophagus_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Esophagus_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Thyroid_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Thyroid_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Adrenal_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Adrenal_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Cervix_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Cervix_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  Skeletal_Atlasgenes = read.csv(file = '../../HPA_SignatureGenes/Skeletal_Atlasgenes.csv', row.names = 1, header = TRUE,sep = ',')
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  # The two groups of tissues that we want to compare
  
  group1 = Bladder_Data #Reference tissue
  
  
  if(counter1 == 1){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Liver_Data))#Reference tissue Atlas genes
    group2 =  Liver_Data # Tissue compared to the Reference  
  }
  else if(counter1 == 2){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Kidney_Data))#Reference tissue Atlas genes
    group2 =  Kidney_Data # Tissue compared to the Reference 
  }
  else if(counter1 == 3){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Small_Intestine_Data))#Reference tissue Atlas genes
    group2 =  Small_Intestine_Data # Tissue compared to the Reference  
  }
  else if(counter1 == 4){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Stomach_Data))#Reference tissue Atlas genes
    group2 =  Stomach_Data # Tissue compared to the Reference  
  }  
  
  else if(counter1 == 5){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Skin_Data))#Reference tissue Atlas genes
    group2 =  Skin_Data # Tissue compared to the Reference  
  }  
  
  else if(counter1 == 6){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Brain_Data))#Reference tissue Atlas genes
    group2 =  Brain_Data # Tissue compared to the Reference  
  }  
  
  
  else if(counter1 == 7){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Pancreas_Data))#Reference tissue Atlas genes
    group2 =  Pancreas_Data # Tissue which will be compared to the Reference  
  }    
  
  else if(counter1 == 8){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Lung_Data))#Reference tissue Atlas genes
    group2 =  Lung_Data # Tissue compared to the Reference  
  }
  else if(counter1 == 9){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Ovary_Data))#Reference tissue Atlas genes
    group2 =  Ovary_Data # Tissue compared to the Reference  
  }
  else if(counter1 == 10){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Prostate_Data))#Reference tissue Atlas genes
    group2 =  Prostate_Data # Tissue compared to the Reference  
  }
  else if(counter1 == 11){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Bladder_Data))#Reference tissue Atlas genes
    group2 =  Bladder_Data # Tissue compared to the Reference  
  }  
  
  else if(counter1 == 12){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Esophagus_Data))#Reference tissue Atlas genes
    group2 =  Esophagus_Data # Tissue compared to the Reference  
  }  
  
  else if(counter1 == 13){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Thyroid_Data))#Reference tissue Atlas genes
    group2 =  Thyroid_Data # Tissue  compared to the Reference  
  }    
  
  
  else if(counter1 == 14){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Adrenal_Data))#Reference tissue Atlas genes
    group2 =  Adrenal_Data # Tissue  compared to the Reference  
  } 
  
  else if(counter1 == 15){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Cervix_Data))#Reference tissue Atlas genes
    group2 =  Cervix_Data # Tissue  compared to the Reference  
  } 
  
  else if(counter1 == 16){
    Atlas_Gold_Standard_Genes = intersect(Bladder_Atlasgenes$Ensembl, rownames(Skeletal_Data))#Reference tissue Atlas genes
    group2 =  Skeletal_Data # Tissue  compared to the Reference  
  } 
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  #filter out genes with many zeros accross samples.
  
  combinedDatasets <- cbind(group1,group2)
  keep <- which(rowSums(combinedDatasets==0) <= round(0.6*dim(combinedDatasets)[2]))
  group1<-group1[keep,]
  group2<-group2[keep,]
  Atlas_Gold_Standard_Genes = intersect(Atlas_Gold_Standard_Genes, rownames(group1))
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  # constructing the Atlas dataset for group1!
  
  
  keeprow <- c()
  for(i in 1 : length(Atlas_Gold_Standard_Genes)){
    index <- which(rownames(group1) == Atlas_Gold_Standard_Genes[i])
    keeprow <- rbind(keeprow,index)
  }
  
  group1_Atlas_DataSet <- group1[keeprow,]+1
  group1_Atlas_DataSet <- t(group1_Atlas_DataSet) # we transpose the matrix such as rows being the samples and columns the genes
  
  # constructing the Atlas probability vector for group1!
  
  Prob_group1_Atlas_DataSet <- group1_Atlas_DataSet
  for (i in 1 : dim(group1_Atlas_DataSet)[1]){ # for each tissue
    for (j in 1 : dim(group1_Atlas_DataSet)[2]){ # for each gene
      Prob_group1_Atlas_DataSet[i,j] <- group1_Atlas_DataSet[i,j]/sum(group1_Atlas_DataSet[i,])
    }
  }
  
  
  V <- dim(Prob_group1_Atlas_DataSet)[1] # number of samples in group1
  mean_group1 <- colMeans(Prob_group1_Atlas_DataSet)  # the means of the genes using probability distribution matrix
  variance_group1 <- colVars(Prob_group1_Atlas_DataSet) # the variances of the genes using probability distribution matrix
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  # calculate the likelihood (confidence) about the probability of appearace of each Atlas gene
  
  #calculate the uncertainty
  phi_group1 <- matrix(data = NA, nrow = dim(Prob_group1_Atlas_DataSet)[1], ncol = dim(Prob_group1_Atlas_DataSet)[2]) # The uncertainties of the genes
  for(i in 1:dim(phi_group1)[1]){ # for each tissue
    for(j in 1:dim(phi_group1)[2]){ # for each gene
      phi_group1[i,j] <- log10((1/sqrt(2*pi*variance_group1[j]))*(exp(-((Prob_group1_Atlas_DataSet[i,j]-mean_group1[j])^2)/(2*variance_group1[j])))+1)
    }
  }
  
  
  #normalize the weights of each tissue such as to sum to 1 accross genes
  for (i in 1 : dim(phi_group1)[1]){
    phi_group1[i,] <- phi_group1[i,]/sum(phi_group1[i,])
  }
  
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  
  
  
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  # constructing the Atlas dataset for group2!
  
  keeprow <- c()
  for(i in 1 : length(Atlas_Gold_Standard_Genes)){
    index <- which(rownames(group2) == Atlas_Gold_Standard_Genes[i])
    keeprow <- rbind(keeprow,index)
  }
  
  group2_Atlas_DataSet <- group2[keeprow,]+1
  group2_Atlas_DataSet <- t(group2_Atlas_DataSet) # we transpose the matrix such as rows being the samples and columns the genes
  
  
  # constructing the Atlas probability vector for group2!
  Prob_group2_Atlas_DataSet <- group2_Atlas_DataSet
  for (i in 1 : dim(group2_Atlas_DataSet)[1]){ # for each tissue
    for (j in 1 : dim(group2_Atlas_DataSet)[2]){ # for each gene
      Prob_group2_Atlas_DataSet[i,j] <- group2_Atlas_DataSet[i,j]/sum(group2_Atlas_DataSet[i,])
    }
  }
  
  
  U <- dim(Prob_group2_Atlas_DataSet)[1] # number of samples in group2
  
  mean_group2 <- colMeans(Prob_group2_Atlas_DataSet)  # the means of the genes using probability distribution matrix
  variance_group2 <- colVars(Prob_group2_Atlas_DataSet) # the variances of the genes using probability distribution matrix
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  # calculate the likelihood (confidence) about the probability of appearace of each Atlas gene
  
  #calculate the uncertainty
  phi_group2 <- matrix(data = NA, nrow = dim(Prob_group2_Atlas_DataSet)[1], ncol = dim(Prob_group2_Atlas_DataSet)[2]) # The uncertainties of the genes
  
  for(i in 1:dim(phi_group2)[1]){ # for each tissue
    for(j in 1:dim(phi_group2)[2]){ # for each gene
      phi_group2[i,j] <- log10((1/sqrt(2*pi*variance_group2[j]))*(exp(-((Prob_group2_Atlas_DataSet[i,j]-mean_group2[j])^2)/(2*variance_group2[j])))+1)
      
    }
  }
  
  
  #normalize the weights of each tissue such as to sum to 1 accross genes
  for (i in 1 : dim(phi_group2)[1]){
    phi_group2[i,] <- phi_group2[i,]/sum(phi_group2[i,])
  }
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  
  
  
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  # calculating the confidences that will be used for the calculation of the weights w_i and w_j for the JSD
  
  #step 1 calculate the standarized expression metrices of the atlas genes (zscore)
  K1<-scale(as.matrix(t(group1_Atlas_DataSet)), center = TRUE, scale = TRUE)
  x1 <- t(K1) # transpose the matrix
  #Here we fixed the optimal value of lambda as calculated by STARS in order to speed up the execution of the algorithm.
  #If you want to run StARs, comment the following lines and uncomment the lines below them.
  out.glasso = huge(x1, method = "glasso",lambda =  0.3 ,cov.output = TRUE) 
  inv_cov <- out.glasso$icov[[1]]
  inv_cov <- matrix(inv_cov,nrow = dim(K1)[1],ncol = dim(K1)[1])
  cov_Matr <- out.glasso$cov[[1]]
  cov_Matr <- matrix(cov_Matr,nrow = dim(K1)[1],ncol = dim(K1)[1])
  
  
  # estimate the confidence for each tissue to be sampled from a normal distribution with parameters N(0,Sigma).
  zeta_group1 <- c()
  for(i in 1:dim(K1)[2]){ # for each tissue in group_1 calculate the confidence (log-likelihood) 
    temp <- log10(exp(-0.5*sqrt(as.vector(K1[,i])%*%inv_cov%*%as.vector(K1[,i])))/(((2*pi)^(1/dim(K1)[1]))*sqrt(det(cov_Matr)))+1)
    zeta_group1 <- cbind(zeta_group1,temp)
  }
  
  
  
  
  
  
  # ## This part uses StARS to select the optimal sparsity parameter. This makes the execution slower!
  # out.glasso = huge(x1, method = "glasso",lambda = seq(from = 0.3, to = 0.4, by = 0.1),cov.output = TRUE)
  # ##model selection using stars
  # out.select = huge.select(out.glasso, criterion = "stars", stars.thresh = 0.05,rep.num= 5)
  # inv_cov <- out.select$opt.icov
  # index = which(out.select$lambda == out.select$opt.lambda)
  
  
  
  ## estimate the confidence for each tissue to be sampled from a normal distribution with parameters N(0,Sigma).
  #zeta_group1 <- c()
  #for(i in 1:dim(K1)[2]){ # for each tissue in group_1 calculate the confidence (log-likelihood) 
  #temp <- log(exp(-0.5*sqrt(as.vector(K1[,i])%*%inv_cov%*%as.vector(K1[,i])))/(((2*pi)^(1/dim(K1)[1]))*sqrt(det(out.select$cov[[index]]))))
  #zeta_group1 <- cbind(zeta_group1,temp)
  #}
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  # calculating the confidences that will be used for the calculation of the weights w_i and w_j for the JSD
  
  #step 1 calculate the standarized expression metrices of the atlas genes (zscore)
  K2<-scale(as.matrix(t(group2_Atlas_DataSet)), center = TRUE, scale = TRUE)
  x2 <- t(K2) # transpose the matrix
  #Here we fixed the optimal value of lambda as calculated by STARS in order to speed up the execution of the algorithm.
  #If you want to run StARs, comment the following lines and uncomment the lines below them.
  out.glasso = huge(x2, method = "glasso",lambda =  0.3 ,cov.output = TRUE)
  inv_cov <- out.glasso$icov[[1]]
  inv_cov <- matrix(inv_cov,nrow = dim(K2)[1],ncol = dim(K2)[1])
  cov_Matr <- out.glasso$cov[[1]]
  cov_Matr <- matrix(cov_Matr,nrow = dim(K2)[1],ncol = dim(K2)[1])
  
  
  # estimate the confidence for each tissue to be sampled from a normal distribution with parameters N(0,Sigma).
  zeta_group2 <- c()
  for(i in 1:dim(K2)[2]){ # for each tissue in group_1 calculate the confidence (log10-likelihood) 
    temp <- log10(exp(-0.5*sqrt(as.vector(K2[,i])%*%inv_cov%*%as.vector(K2[,i])))/(((2*pi)^(1/dim(K2)[1]))*sqrt(det(cov_Matr)))+1)
    zeta_group2 <- cbind(zeta_group2,temp)
  }
  
  
  
  ## This part uses StARS to select the optimal sparsity parameter. This makes the execution slower!
  #out.glasso = huge(x1, method = "glasso",lambda = seq(from = 0.3, to = 0.4, by = 0.1),cov.output = TRUE)
  ##model selection using stars
  #out.select = huge.select(out.glasso, criterion = "stars", stars.thresh = 0.05,rep.num= 5)
  #inv_cov <- out.select$opt.icov
  #index = which(out.select$lambda == out.select$opt.lambda)
  
  
  
  ## estimate the confidence for each tissue to be sampled from a normal distribution with parameters N(0,Sigma).
  #zeta_group2 <- c()
  #for(i in 1:dim(K2)[2]){ # for each tissue in group_2 calculate the confidence (log-likelihood) 
  #temp <- log(exp(-0.5*sqrt(as.vector(K2[,i])%*%inv_cov%*%as.vector(K2[,i])))/(((2*pi)^(1/dim(K2)[1]))*sqrt(det(out.select$cov[[index]]))))
  #zeta_group2 <- cbind(zeta_group2,temp)
  #}
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  # Weighted ranking Correlation
  
  
  # constructing the Ranking matrix of the atlas genes for group1!
  
  Rankings_Atlas_group1 <- matrix(data = 0, nrow = length(Atlas_Gold_Standard_Genes), ncol =dim(group1)[2]) # The uncertainties of the genes
  
  for(i in 1 : dim(group1)[2]){ # for each sample
    #we order the matrix of expressions in decreasing order
    keepRow<-c()
    group1_Rankings <- group1[order(group1[,i],decreasing = TRUE),]
    for(j in 1 : length(Atlas_Gold_Standard_Genes)){ # for each atlas gene
      index <- which(rownames(group1_Rankings) == Atlas_Gold_Standard_Genes[j])
      Rankings_Atlas_group1[j,i] <- index  
    }
    
  }
  
  
  
  Rankings_Atlas_group1 <- t(Rankings_Atlas_group1) # transpose matrix such as the columns to correspond to Atlas genes
  
  
  ######################################################################################################################
  # calculate the likelihood (confidence) about the probability of appearace of each Atlas gene
  
  
  meanRank_group1 <- colMeans(Rankings_Atlas_group1)  # mean rankings of the Atlas genes accross tissues in group 1
  varianceRank_group1 <- colVars(Rankings_Atlas_group1) # the ranking variances of the Atlas genes accross tissues in group 1
  
  #calculate the uncertainty
  tau_group1 <- matrix(data = NA, nrow = dim(Rankings_Atlas_group1)[1], ncol = dim(Rankings_Atlas_group1)[2]) # The uncertainties of the genes
  
  for(i in 1:dim(tau_group1)[1]){ # for each tissue
    for(j in 1:dim(tau_group1)[2]){ # for each gene
      tau_group1[i,j] <- log10((1/sqrt(2*pi*varianceRank_group1[j]))*(exp(-((Rankings_Atlas_group1[i,j]-meanRank_group1[j])^2)/(2*varianceRank_group1[j])))+1)
      
    }
  }
  
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  
  
  # constructing the Ranking matrix of the atlas genes for group2!
  Rankings_Atlas_group2 <- matrix(data = 0, nrow = length(Atlas_Gold_Standard_Genes), ncol =dim(group2)[2]) # The uncertainties of the genes
  
  for(i in 1 : dim(group2)[2]){ # for each sample
    #we order the matrix of expressions in decreasing order
    keepRow<-c()
    group2_Rankings <- group2[order(group2[,i],decreasing = TRUE),]
    for(j in 1 : length(Atlas_Gold_Standard_Genes)){ # for each atlas gene
      index <- which(rownames(group2_Rankings) == Atlas_Gold_Standard_Genes[j])
      Rankings_Atlas_group2[j,i] <- index  
    }
  }
  
  Rankings_Atlas_group2 <- t(Rankings_Atlas_group2) # transpose matrix such as the columns to correspond to Atlas genes
  
  
  
  ######################################################################################################################
  # calculate the likelihood (confidence) about the probability of appearace of each Atlas gene
  
  
  meanRank_group2 <- colMeans(Rankings_Atlas_group2)  # mean rankings of the Atlas genes accross tissues in group 2
  varianceRank_group2 <- colVars(Rankings_Atlas_group2) # the ranking variances of the Atlas genes accross tissues in group 2
  
  #calculate the uncertainty
  tau_group2 <- matrix(data = NA, nrow = dim(Rankings_Atlas_group2)[1], ncol = dim(Rankings_Atlas_group2)[2]) # The uncertainties of the genes
  
  for(i in 1:dim(tau_group2)[1]){ # for each tissue
    for(j in 1:dim(tau_group2)[2]){ # for each gene
      tau_group2[i,j] <- log10((1/sqrt(2*pi*varianceRank_group2[j]))*(exp(-((Rankings_Atlas_group2[i,j]-meanRank_group2[j])^2)/(2*varianceRank_group2[j])))+1)
      
    }
  }
  
  
  
  
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  ######################################################################################################################
  
  
  
  
  
  # calculate pairwise the JSD distances between the tissue groups
  SRwJSD <-c()
  # calculate pairwise the wRCD distances between the tissue groups
  wRCD <-c() #weighted ranking correlation distance
  
  
  
  
  if(counter1 != 11){
    
    for(i in 1: dim(Prob_group1_Atlas_DataSet)[1]){ # for each tissue in group1
      for(j in 1: dim(Prob_group2_Atlas_DataSet)[1]){ # for each tissue in group2
        
        
        w1 = (zeta_group1[i]+log10(dim(K1)[2]))/((zeta_group1[i]+log10(dim(K1)[2])) +(zeta_group2[j]+log10(dim(K2)[2])))
        
        w2 = (zeta_group2[j]+log10(dim(K2)[2]))/((zeta_group1[i]+log10(dim(K1)[2])) +(zeta_group2[j]+log10(dim(K2)[2])))
        
        
        # rescaling probability vectors
        Prob_group1_Atlas_DataSet_reschaled <- c()
        for(q in 1 :dim(Prob_group1_Atlas_DataSet)[2]){# for each gene
          temp <- (phi_group1[i,q]*Prob_group1_Atlas_DataSet[i,q])/sum((phi_group1[i,]*Prob_group1_Atlas_DataSet[i,]))
          Prob_group1_Atlas_DataSet_reschaled <- cbind(Prob_group1_Atlas_DataSet_reschaled,temp)
          
        }
        
        
        
        # rescaling probability vectors
        Prob_group2_Atlas_DataSet_reschaled <- c()
        for(q in 1 :dim(Prob_group2_Atlas_DataSet)[2]){
          temp <- (phi_group2[j,q]*Prob_group2_Atlas_DataSet[j,q])/sum((phi_group2[j,]*Prob_group2_Atlas_DataSet[j,]))
          Prob_group2_Atlas_DataSet_reschaled <- cbind(Prob_group2_Atlas_DataSet_reschaled,temp)
          
        }
        
        
        
        rescal =rbind(Prob_group1_Atlas_DataSet_reschaled,Prob_group2_Atlas_DataSet_reschaled)
        rescal_JSD <- gJSD(rescal, unit = "log2", weights = c(w1,w2), est.prob = NULL)
        
        
        
        SRwJSD <- cbind(SRwJSD,sqrt(rescal_JSD))
        
        
        
        # calculate ranking correlation distance      
        beta <- c()
        for(q in 1:dim(Rankings_Atlas_group2)[2]){ # for all Atlas genes
          tempa <- tau_group1[i,q] + tau_group2[j,q]
          beta <- cbind(beta,tempa)
        }
        
        
        #note: ranking correlation distance
        temp3 <- sqrt(1 - max(0,weightedCorr(Rankings_Atlas_group1[i,], Rankings_Atlas_group2[j,], method =  "Pearson",weights = beta)))
        wRCD <- cbind(wRCD,temp3)
        
      }
    }
    
  }
  
  
  else if(counter1 == 11){
    
    for(i in 1: (dim(Prob_group1_Atlas_DataSet)[1]-1)){ # for each tissue in group1
      for(j in (i+1) : dim(Prob_group2_Atlas_DataSet)[1]){ # for each tissue in group2
        
        
        w1 = (zeta_group1[i]+log10(dim(K1)[2]))/((zeta_group1[i]+log10(dim(K1)[2])) +(zeta_group2[j]+log10(dim(K2)[2])))
        w2 = (zeta_group2[j]+log10(dim(K2)[2]))/((zeta_group1[i]+log10(dim(K1)[2])) +(zeta_group2[j]+log10(dim(K2)[2])))
        
        
        # rescaling probability vectors
        Prob_group1_Atlas_DataSet_reschaled <- c()
        for(q in 1 :dim(Prob_group1_Atlas_DataSet)[2]){# for each gene
          temp <- (phi_group1[i,q]*Prob_group1_Atlas_DataSet[i,q])/sum((phi_group1[i,]*Prob_group1_Atlas_DataSet[i,]))
          Prob_group1_Atlas_DataSet_reschaled <- cbind(Prob_group1_Atlas_DataSet_reschaled,temp)
          
        }
        
        
        
        # rescaling probability vectors
        Prob_group2_Atlas_DataSet_reschaled <- c()
        for(q in 1 :dim(Prob_group2_Atlas_DataSet)[2]){
          temp <- (phi_group2[j,q]*Prob_group2_Atlas_DataSet[j,q])/sum((phi_group2[j,]*Prob_group2_Atlas_DataSet[j,]))
          Prob_group2_Atlas_DataSet_reschaled <- cbind(Prob_group2_Atlas_DataSet_reschaled,temp)
          
        }
        
        
        
        
        rescal =rbind(Prob_group1_Atlas_DataSet_reschaled,Prob_group2_Atlas_DataSet_reschaled)
        rescal_JSD <- gJSD(rescal, unit = "log2", weights = c(w1,w2), est.prob = NULL)
        
        
        
        SRwJSD <- cbind(SRwJSD,sqrt(rescal_JSD))
        
        
        
        # calculate ranking correlation distance
        beta <- c()
        for(q in 1:dim(Rankings_Atlas_group2)[2]){ # for all Atlas genes
          tempa <- tau_group1[i,q] + tau_group2[j,q]
          beta <- cbind(beta,tempa)
        }
        
        
        #note: ranking correlation distance
        temp3 <- sqrt(1 - max(0,weightedCorr(Rankings_Atlas_group1[i,], Rankings_Atlas_group2[j,], method =  "Pearson",weights = beta)))
        wRCD <- cbind(wRCD,temp3)
        
      }
    }
    
  }
  
  
  
  # Transcriptomic signature distance
  wTSDM <- (0.5*SRwJSD)+(0.5*wRCD)
  
  
  if(counter1 == 1){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_liver_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_liver_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_liver_SRwJSD.csv",row.names=FALSE)
  }
  else if(counter1 == 2){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_Kidney_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_Kidney_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_Kidney_SRwJSD.csv",row.names=FALSE) 
  }
  else if(counter1 == 3){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_smallIntestine_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_smallIntestine_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_smallIntestine_SRwJSD.csv",row.names=FALSE) 
  }
  else if(counter1 == 4){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_stomach_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_stomach_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_stomach_SRwJSD.csv",row.names=FALSE) 
  }  
  
  else if(counter1 == 5){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_skin_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_skin_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_skin_SRwJSD.csv",row.names=FALSE) 
  }  
  
  
  else if(counter1 == 6){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_brain_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_brain_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_brain_SRwJSD.csv",row.names=FALSE) 
  }   

    
  else if(counter1 == 7){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_pancreas_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_pancreas_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_pancreas_SRwJSD.csv",row.names=FALSE) 
  }  
  
  
  else if(counter1 == 8){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_lung_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_lung_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_lung_SRwJSD.csv",row.names=FALSE)
  } 
  
  else if(counter1 == 9){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_ovary_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_ovary_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_ovary_SRwJSD.csv",row.names=FALSE) 
  }
  else if(counter1 == 10){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_prostate_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_prostate_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_prostate_SRwJSD.csv",row.names=FALSE) 
  }
  else if(counter1 == 11){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_bladder_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_bladder_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_bladder_SRwJSD.csv",row.names=FALSE) 
  }  
  
  else if(counter1 == 12){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_esophagus_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_esophagus_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_esophagus_SRwJSD.csv",row.names=FALSE) 
  }  
  
  
  else if(counter1 == 13){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_thyroid_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_thyroid_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_thyroid_SRwJSD.csv",row.names=FALSE) 
  }    
  
  else if(counter1 == 14){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_Adrenal_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_Adrenal_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_Adrenal_SRwJSD.csv",row.names=FALSE) 
  }   
  else if(counter1 == 15){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_Cervix_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_Cervix_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_Cervix_SRwJSD.csv",row.names=FALSE) 
  }   
  else if(counter1 == 16){
    write.csv(wTSDM, file = "../../Results/Heatmap_wTSDM/Bladder_vs_Skeletal_wTSDM.csv",row.names=FALSE)
    write.csv(wRCD, file = "../../Results/Heatmap_wRCD/Bladder_vs_Skeletal_wRCD.csv",row.names=FALSE)
    write.csv(SRwJSD, file = "../../Results/Heatmap_SRwJSD/Bladder_vs_Skeletal_SRwJSD.csv",row.names=FALSE) 
  }  
  
  print(counter1)
  
}


######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################




