
## according to the predictors to their relevance for the classification
## wi=sqrt(sum_j=1^J sum_k=1^Mi [P(Y=j|Xi=k)-P(Y=j)]^2)
library(infotheo)
library(fpc)

LF_train_releweight <- function(data,train_index){
  result=list()
  synornot_train=c(rep(1,831),rep(-1,831))[train_index]
  result$synornot_train = synornot_train
  data_all =NULL
  train_data=data
  for(i in 1:length(data)){
    data[[i]]<- data[[i]][,c(1,train_index[train_index<=831]+1)]
    train_data[[i]] <- data[[i]][train_index,]
    data_all <- cbind(data_all,as.numeric(apply(data[[i]][,-1],1,max)))
    data_t <- as.matrix(data_all[train_index,])
  }
  result$data=data
  result$train_data = train_data
  result$data_all=data_all
  #colnames(data_t)<-c('synornot','F1','F2','F3','F4','F5','F6','F7')
  train_pos = data_t[which(synornot_train==1),]
  train_neg = data_t[which(synornot_train==-1),] 

  result$train_pos=train_pos
  result$train_neg=train_neg


  weight=as.numeric(ncol(data_t))
  mi=as.numeric(ncol(data_t))
  for(i in 1:ncol(data_t)){
    if(length(unique(data_t[,i]))==2){
      dis_data=data_t[,i]
    }else{
      max_k  = 30
      silhouette_scores <- numeric(max_k)
      hclust_fit <- hclust(dist(data_t[,i]), method='average')
      for (k in 2:max_k) {
        hclust_result <- cutree(hclust_fit,k)
       silhouette_scores[k] <- cluster.stats(dist(data_t[,i]), hclust_result)$avg.silwidth
      }

      optimal_k <- max(which.max(silhouette_scores),2)
      hclust_result = cutree(hclust_fit,optimal_k)

      dis_data=hclust_result
    }

    contingency = as.matrix(table(dis_data,synornot_train)) 

    row_sums_broadcasted <- matrix(rowSums(contingency), nrow = nrow(contingency), ncol = ncol(contingency), byrow = FALSE)
    col_sums_broadcasted <- matrix(colSums(contingency), nrow = nrow(contingency), ncol = ncol(contingency), byrow = TRUE)

    com_ele <- contingency / row_sums_broadcasted - col_sums_broadcasted/sum(contingency)
  
    squared_df <- sum(apply(com_ele, MARGIN = c(1, 2), FUN = function(x) x^2))
    weight[i] <- squared_df/optimal_k
    mi[i]=mutinformation(dis_data, synornot_train)
  }

  result$weight=weight
  result$data_all=data_all
  result$train_pos=train_pos
  result$train_neg=train_neg
  data_t=data_t
  LR_train=matrix(NA,nrow(data_t),ncol(data_t))
  h1=h2=vector()
 
  h1 <- apply(train_pos, 2, function(col) density(col)$bw)
  h2 <- apply(train_neg, 2, function(col) density(col)$bw)
  #h1=h2=apply(data_t,2,function(col) density(col)$bw)
  
  result$h1=h1
  result$h2=h2
  LR_train_a=apply(LR_train,1,prod)
  result$normLR=LR_train_a
  

  return(result)
}



##### testing 
LF_test_releweight<-function(data_all,test_index,train_pos,train_neg,h1,h2,LR_rand,weight){
  #p=0
  result=list()
  synornot = c(rep(1,831),rep(-1,831))
  data_test = as.matrix(data_all[test_index,])
  synornot_test <- synornot[test_index]
  result$synornot_test = synornot_test
  LR_test=matrix(NA,nrow(data_test),ncol(data_test))
  
  
  for (k in 1:ncol(data_test)) {
    P1 <- sapply(data.frame(data_test)[, k], function(x) {
      sum(exp(-(x - data.frame(train_pos)[, k])^2 / (2 * h1[k]^2))) / (nrow(data.frame(train_pos)) * h1[k] * sqrt(2 * pi))+10^(-8)
    })
    P2 <- sapply(data.frame(data_test)[, k], function(x) {
      sum(exp(-(x - data.frame(train_neg)[, k])^2 / (2 * h2[k]^2))) / (nrow(data.frame(train_pos)) * h2[k] * sqrt(2 * pi))+10^(-8)
    })
    LR_test[, k] <- (P1 / P2) ^ weight[k]
  }

  
  LR_test_a=apply(LR_test,1,prod)
  
  result$LR_test_a = LR_test_a

  return(result)
}

auc(as.numeric(synornot[test_index]),as.numeric(LR_test_a))


library(pROC)
require('caret')
library(PRROC)
set.seed(1234)

folds1<-createFolds(c(rep(1,831),rep(-1,831)), k = 10, list = TRUE, returnTrain = FALSE) 





  # input cut
  auc_con=aupr_con=kappa=accuracy=sensitivity=f1_score=precision=recall=vector()
  result=NULL
  y_true_all = NULL
  y_score_all = NULL
  auc_con=aupr_con=NULL
  weight=NULL
  #data=list(F1_CV)
  data=list(F1_CV,F2_CV,F3_CV,F4_CV,F5_CV,F6_CV,F7_CV)
  roc_data =NULL
  roc_curves=list()
  for(i in 1:10){

    train_index <- c(1:1662)[-folds[[i]]]
    test_index <- folds[[i]]
   
    train=LF_train_releweight(data,train_index)
    test=LF_test_releweight(train$data_all,test_index,train$train_pos,train$train_neg,train$h1,train$h2,train$LR_rand,train$weight)
    
    y_true <- test$synornot_test
    y_score <- test$LR_test_a
    y_score <- 2/(1+exp(-y_score))-1
    weight <- rbind(weight,train$weight)
    
    y_true[which(y_true==-1)]=0
    
    auc_con=rbind(auc_con,auc(y_true,y_score))
    roc_data=c(roc_data,roc(y_true,y_score))
    aupr_con <- rbind(aupr_con,pr.curve(scores.class0 = y_score, weights.class0 = y_true)$auc.integral)
    
    roc_curve=roc(y_true,y_score)
    roc_curves[[i]]=roc_curve

    y_score_bin <- ifelse(y_score>cut, 1, 0)
 
    confusion_mat <- confusionMatrix(factor(y_score_bin), factor(y_true))
    kappa <- rbind(kappa,confusion_mat$overall["Kappa"])
    
    accuracy <- rbind(accuracy,confusion_mat$overall["Accuracy"])
    
    precision <- rbind(precision,confusion_mat$byClass["Precision"])
    recall <- rbind(recall,confusion_mat$byClass["Sensitivity"])
    
    f1_score <- rbind(f1_score,2 * (precision * recall) / (precision + recall))
    
    result=rbind(result,cbind(test_index,y_true,y_score))
  }
 
  mean=c(mean(auc_con),mean(aupr_con),mean(accuracy),mean(f1_score),mean(precision),mean(recall),mean(kappa))
  sd = c(sd(auc_con),sd(aupr_con),sd(accuracy),sd(f1_score),sd(precision),sd(recall),sd(kappa))
  print(rbind(mean,sd))

  
  

####################################################################################
#####################################################################################
## mutual information
library(infotheo)
LF_train_MIweight <- function(data,train_index){
  result=list()
  synornot_train=c(rep(1,831),rep(-1,831))[train_index]
  result$synornot_train = synornot_train
  data_all =NULL
  train_data=data
  for(i in 1:length(data)){
    data[[i]]<- data[[i]][,c(1,train_index[train_index<831]+1)]
    train_data[[i]] <- data[[i]][train_index,]
    data_all <- cbind(data_all,as.numeric(apply(data[[i]][,-1],1,max)))
    data_t <- data_all[train_index,]
  }

  mi=numeric(ncol(data_t))
  for(i in 1:ncol(data_t)){
    if(length(unique(data_t[,i]))==2){
      dis_data=data_t[,i]
    }else{
      max_k <- 30
      silhouette_scores <- numeric(max_k)
      hclust_fit <- hclust(dist(data_t[,i]), method='average')
      for (k in 2:max_k) {
        hclust_result <- cutree(hclust_fit,k)
        silhouette_scores[k] <- cluster.stats(dist(data_t[,i]), hclust_result)$avg.silwidth
      }

      optimal_k <- max(which.max(silhouette_scores),2)
      hclust_result = cutree(hclust_fit,optimal_k)

      dis_data=hclust_result
    }
    mi[i]=mutinformation(dis_data, synornot_train)
  }
 

 weight=mi
 result$weight=weight
 result$data_all=data_all
 train_pos = data_t[which(synornot_train==1),]
 train_neg = data_t[which(synornot_train==-1),] 
 result$train_pos=train_pos
 result$train_neg=train_neg
 data_t=data_t
 LR_train=matrix(NA,nrow(data_t),ncol(data_t))
 h1=h2=vector()
  
  LR_train=matrix(NA,nrow(data_t),ncol(data_t))
  h1=h2=vector()
  
  for(k in 1:ncol(data_t)){
  
    h1[k]=density(train_pos[,k])$bw
    h2[k]=density(train_neg[,k])$bw
    
    for(i in 1:nrow(data_t)){
      if(length(unique(train_pos[,k]))==2){
        P1=sum(data_t[i,k]==train_pos[,k])
        P2=sum(data_t[i,k]==train_neg[,k])
      }else{
        P1<-density2(train_pos[,k],h1[k],data_t[i,k])+10^(-8)  
        P2<-density2(train_neg[,k],h2[k],data_t[i,k])+10^(-8)  

      }
      LR_train[i,k]<-(P1/P2)**weight[k]
    }
  }
  
 
  result$h1=h1
  result$h2=h2
  LR_train_a=apply(LR_train,1,prod)
  result$normLR=LR_train_a
 
  
  return(result)
}



##### testing 
LF_test_MIweight<-function(data_all,test_index,train_pos,train_neg,h1,h2,LR_rand,weight){
 
  result=list()
  synornot = c(rep(1,831),rep(-1,831))
  data_test = data_all[test_index,]
  synornot_test <- synornot[test_index]
  result$synornot_test = synornot_test
  LR_test=matrix(NA,nrow(data_test),ncol(data_test))
  
  for(k in 1:ncol(data_test)){
   
    for(i in 1:nrow(data_test)){
      if(length(unique(train_pos[,k]))==2){
        P1=sum(data_test[i,k]==train_pos[,k])/nrow(train_pos)
        P2=sum(data_test[i,k]==train_neg[,k])/nrow(train_neg)
      }else{
        P1<-density2(train_pos[,k],h1[k],data_test[i,k])+10^(-8)  #posterior probability given positive
        P2<-density2(train_neg[,k],h2[k],data_test[i,k])+10^(-8)  #posterior probability given false
        
      }
      LR_test[i,k]<-(P1/P2)**weight[k] 
    }
  }
  
  
  
  LR_test_a=apply(LR_test,1,prod)
  
  result$LR_test_a = LR_test_a
 
  return(result)
}





#input cut
  auc_con=aupr_con=kappa=accuracy=sensitivity=f1_score=precision=recall=vector()
  result=NULL
  y_true_all = NULL
  y_score_all = NULL
  auc_con=aupr_con=NULL
  data=list(F1_CV,F2_CV,F3_CV,F4_CV,F5_CV,F6_CV,F7_CV)

  for(i in 1:10){
   
    train_index <- c(1:1662)[-folds[[i]]]
    test_index <- folds[[i]]
    
    train=LF_train_MIweight(data,train_index)
    test=LF_test_MIweight(train$data_all,test_index,train$train_pos,train$train_neg,train$h1,train$h2,train$LR_rand,train$weight)
    
    y_true <- test$synornot_test
    y_score <- test$LR_test_a
   
   
    y_true[which(y_true==-1)]=0
  
    auc_con=rbind(auc_con,auc(y_true,y_score))
    roc_data=c(roc_data,roc(y_true,y_score))
    aupr_con <- rbind(aupr_con,pr.curve(scores.class0 = y_score, weights.class0 = y_true)$auc.integral)
    
    roc_curve=roc(y_true,y_score)
    roc_curves[[i]]=roc_curve
   
    y_score_bin <- ifelse(y_score>cut, 1, 0)
   
    confusion_mat <- confusionMatrix(factor(y_score_bin), factor(y_true))
    kappa <- rbind(kappa,confusion_mat$overall["Kappa"])
    
   
    accuracy <- rbind(accuracy,confusion_mat$overall["Accuracy"])
    
    
    precision <- rbind(precision,confusion_mat$byClass["Precision"])
    recall <- rbind(recall,confusion_mat$byClass["Sensitivity"])
    
    f1_score <- rbind(f1_score,2 * (precision * recall) / (precision + recall))
    
    result=rbind(result,cbind(test_index,y_true,y_score))
  }
  
  
  mean=c(mean(auc_con),mean(aupr_con),mean(accuracy),mean(f1_score),mean(precision),mean(recall),mean(kappa))
  sd = c(sd(auc_con),sd(aupr_con),sd(accuracy),sd(f1_score),sd(precision),sd(recall),sd(kappa))
  print(rbind(mean,sd))
 