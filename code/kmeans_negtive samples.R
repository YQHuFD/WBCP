##### data preprocession
### input F1-F7


diag(F7[1:831,])=0

feature_max=cbind(as.numeric(apply(F1,1,max)),
                  as.numeric(apply(F2,1,max)),
                  as.numeric(apply(F3,1,max)),
                  as.numeric(apply(F4,1,max)),
                  as.numeric(apply(F5,1,max)),
                  as.numeric(apply(F6,1,max)),
                  as.numeric(apply(F7,1,max)))

feature_max
write.csv(feature_max,'feature_all_max.csv')
######KMEANS
##Calinsky Criterion
require(vegan)
fit <- cascadeKM(scale(feature_max[-c(1:831),],center=TRUE,scale=TRUE),1,40,iter=20)
#plot(fit,sortg=TRUE,grpmts.plot=TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,])) ##4

print(fit$results[2,]) 
plot(1:30,as.numeric(fit$results[2,]))

## k=4 for negtive samples
library(cluster)
km <- kmeans(feature_max[-c(1:831),],831,nstart=24)
cluster=data.frame(km$cluster)
library('plyr')
cluster[,2]=rownames(cluster) 

choose_neg <- ddply(cluster,.(km.cluster),function(x) x[sample(nrow(x),1),]) ##832
choose_neg <- ddply(cluster,.(km.cluster),function(x) x[sample(nrow(x),208),]) ##832

randomrow=sample(nrow(choose_neg),1)
choose_neg <- choose_neg[-randomrow,] ##delete a random row 
synornot<-c(rep(1,831),rep(-1,831))

index_new_0403=c(1:831,floor(as.numeric(choose_neg[,2])+831))

tiqu_index=floor(as.numeric(choose_neg[,2])+831)
tiqu_index = c(1:831,unlist(tiqu_index))


lapply(tiqu_index, as.numeric)

F1_CV <- cbind(synornot,F1[unlist(tiqu_index),])
F2_CV <- cbind(synornot,F2[unlist(tiqu_index),])
F3_CV <- cbind(synornot,F3[unlist(tiqu_index),])
F4_CV <- cbind(synornot,F4[unlist(tiqu_index),])
F5_CV <- cbind(synornot,F5[unlist(tiqu_index),])
F6_CV <- cbind(synornot,F6[unlist(tiqu_index),])
F7_CV <- cbind(synornot,F7[unlist(tiqu_index),])


