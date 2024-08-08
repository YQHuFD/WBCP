
# input  S1,S2,...,S7 (list format)
# input known drug combinations "valiadation"


Fn<-function(Sn,Snlist,validation){
  Fn<-matrix(NA,nrow(Snlist),nrow(validation))
  
  for(i in 1:nrow(Snlist)){
    for(j in 1:nrow(validation)){
      a=sqrt(Sn[which(rownames(Sn)==Snlist[i,1]),
                which(colnames(Sn)==validation[j,1])]
             *Sn[which(rownames(Sn)==Snlist[i,2]),
                 which(colnames(Sn)==validation[j,2])])##
      b=sqrt(Sn[which(rownames(Sn)==Snlist[i,1]),
                which(colnames(Sn)==validation[j,2])]
             *Sn[which(rownames(Sn)==Snlist[i,2]),
                 which(colnames(Sn)==validation[j,1])])
      Fn[i,j]<-max(a,b)
    }
    print(i)
  }
  rownames(Fn)<-paste(Snlist[,1],Snlist[,2],sep="+")
  colnames(Fn)<-paste(validation[,1],validation[,2],sep="+")
  return(Fn)
}



Fn_optimized <- function(Sn, Snlist, validation) {
  n_Snlist <- nrow(Snlist)
  n_validation <- nrow(validation)
  
  Fn <- matrix(NA, n_Snlist, n_validation)
  
  for (i in 1:n_Snlist) {
    idx1 <- which(rownames(Sn) == Snlist[i, 1])
    idx2 <- which(rownames(Sn) == Snlist[i, 2])
    
    for (j in 1:n_validation) {
  
      idx3 <- which(colnames(Sn) == validation[j, 1])
      idx4 <- which(colnames(Sn) == validation[j, 2])
      
      
      a <- sqrt(Sn[idx1, idx3] * Sn[idx2, idx4])
      b <- sqrt(Sn[idx1, idx4] * Sn[idx2, idx3])
      
     
      Fn[i, j] <- max(a, b)
    }
    print(i)
  }
  
  
  rownames(Fn) <- paste(Snlist[, 1], Snlist[, 2], sep = "+")
  colnames(Fn) <- paste(validation[, 1], validation[, 2], sep = "+")
  
  return(Fn)
}


start_time <- Sys.time()
F7<-Fn(S7,S7list,validation)
write.csv(F7,"F7.csv")
end_time <- Sys.time()





## for example : F7
colnames(F7)<-gsub('\\.','\\+',colnames(F7))
name1=colnames(F7)
name2<-sapply(1:length(name1),function(i){paste(strsplit(name1,split="\\+")[[i]][2],strsplit(name1,split="\\+")[[i]][1],sep="+")})
colnames(F7)[which(name2%in%rownames(F7))]=rownames(F7)[which(rownames(F7)%in%name2)]
common_columns <- intersect(colnames(F1), colnames(F7))

F7 <- F7[, common_columns, drop = FALSE]

F7 <- cbind(F1[,1:2],F7)
write.csv(F7, 'F7.csv')

