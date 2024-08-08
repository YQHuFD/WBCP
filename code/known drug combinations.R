
install.packages("readr")
install.packages("dplyr")

library(readr)
library(dplyr)

## preprocessing DCDB
#input DCDB, DC_USAGE,DC_TO_DCU

DC_efficacious <- DC_USAGE[which(DC_USAGE$EFFICACY=='Efficacious'),]

DC_efficacious <- merge(DC_USAGE, DC_TO_DCU[, c("DC_ID", "DCU_ID")], by.x = "DCU_ID", by.y = "DCU_ID", all.x = TRUE)
DCDB_eff <- DCDB[DCDB$DrugCombination_ID%in%DC_efficacious$DC_ID,]   
write.csv(DCDB_eff,'DCDB_eff.csv')

your_data_split <- cbind(DCDB_eff[,1], do.call(rbind, strsplit(as.character(DCDB_eff[,2]), ";")))
your_data_split <- data.frame(ID = rep(DCDB_eff[,1], sapply(strsplit(as.character(DCDB_eff[,2]), ";"), length)),
                              Value = unlist(strsplit(as.character(DCDB_eff[,2]), ";")))

your_data=your_data_split
library(reshape2)

your_data$row_num <- ave(your_data$ID, your_data$ID, FUN = seq_along)
result_data <- dcast(your_data, ID ~ paste0("Value", row_num), value.var = "Value")[, -1]


your_data$row_num <- ave(your_data$ID, your_data$ID, FUN = function(x) rep(1:ceiling(length(x)/2), each = 2, length.out = length(x)))
result_data <- dcast(your_data, ID ~ paste0("Value", row_num), value.var = "Value")[, -1]



result_data <- your_data_split %>%
  inner_join(your_data, by = "ID") %>%
  filter(Value.x < Value.y) %>%
  select(ID, Drug1 = Value.x, Drug2 = Value.y)


write.csv(result_data,'DCDB_eff_combinations.csv')


#######################################################################


# data1 
DCDB_name <- DCDB[,c('Drug1','Drug2')]
colnames(DCDB_name) <- c('name1','name2')

library(readxl)
ASDCD_name <- ASDCD[,c('name1','name2')]

library('dplyr')
name_union <- union(DCDB_name,ASDCD_name)
name_union<-bind_rows(DCDB_name,ASDCD_name)
newpos_name <- name_union
colnames(newpos_name)=c('X1','X2')

library(stringr)
newpos_name$X1 <- str_trim(newpos_name$X1)
newpos_name$X2 <- str_trim(newpos_name$X2)
newpos_name$name1 <- str_to_lower(newpos_name$X1)
newpos_name$name2 <- str_to_lower(newpos_name$X2)
newpos_name$name1 <-  word(newpos_name$name1, 1)
newpos_name$name2 <-  word(newpos_name$name2, 1)

newpos_pairs <-  paste(newpos_name$name1, newpos_name$name2, sep = "+")
newpos_pairs2 <- paste(newpos_name$name2, newpos_name$name1, sep = "+")

newpos<- c(newpos_pairs,newpos_pairs2)
newpos<-unique(newpos)
sum(rownames(FSNF)%in%newpos)
sum(newpos%in%rownames(FSNF))

positive_575<-intersect(newpos,rownames(F1))
write.csv(positive_575,'positive_575.csv')


positive_575_s <- data.frame(t(apply(positive_575, 1, function(row) sort(row))))


positive_2019NC<-data.frame(colnames(F1_413))
write.csv(positive_2019NC,'positive_2019NC.csv')

positive_831<-union(positive_575,positive_2019NC)
write.csv(positive_831,"positive_831.csv")

result<- strsplit(as.character(positive_831[,1]),'\\+')

po_831 = colnames(positive_831)

pos_831=matrix(NA,831,2)
for(i in 1:nrow(pos_831)){
  pos_831[i,1]<-result[[i]][1]
  pos_831[i,2]<-result[[i]][2]
}
# 对每一行的字符串按照顺序排列
positive_831_s <- data.frame(t(apply(pos_831, 1, function(row) sort(row))))
positive_831_s <- paste(positive_831_s[,1],positive_831_s[,2],sep='+')
unique(positive_831_s)

colnames(F1_413)%in%positive_831_s

write.csv(positive_831_s,'positive_831_s.csv')

sum(positive_831_s[,1]%in%rownames(FSNF))

#####把413的列名重新规划
result<- strsplit(colnames(F1_413[,-c(1,2)]),'\\+')

pos_413 = colnames(positive_413)

pos_413=matrix(NA,413,2)
for(i in 1:nrow(pos_413)){
  pos_413[i,1]<-result[[i]][1]
  pos_413[i,2]<-result[[i]][2]
}
# 对每一行的字符串按照顺序排列
positive_413_s <- data.frame(t(apply(pos_413, 1, function(row) sort(row))))
positive_413_s <- paste(positive_413_s[,1],positive_413_s[,2],sep='+')
unique(positive_413_s)
positive_413_s <- c('synornot','antiornot',positive_413_s)
names(F1_413)<-positive_413_s
names(F2_413)<-positive_413_s
names(F3_413)<-positive_413_s
names(F4_413)<-positive_413_s
names(F5_413)<-positive_413_s
names(F6_413)<-positive_413_s
names(F7_413)<-positive_413_s

rownames(F1_413)<-rownames(F1_418)
rownames(F2_413)<-rownames(F2_418)
rownames(F3_413)<-rownames(F3_418)
rownames(F4_413)<-rownames(F4_418)
rownames(F5_413)<-rownames(F5_418)
rownames(F6_413)<-rownames(F6_418)
rownames(F7_413)<-rownames(F7_418)

write.csv(F1_413,'F1_413.csv')
write.csv(F2_413,'F2_413.csv')
write.csv(F3_413,'F3_413.csv')
write.csv(F4_413,'F4_413.csv')
write.csv(F5_413,'F5_413.csv')
write.csv(F6_413,'F6_413.csv')
write.csv(F7_413,'F7_413.csv')

positive_417_s=positive_831_s[-which(positive_831_s%in%positive_413_s)]
unique(positive_417_s)

colnames(F1_413)%in%positive_831_s
