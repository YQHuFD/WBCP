#通路相似性计算

setwd("E:\\药物联用项目\\操作\\相似性计算\\kegg/")

drug_target = read.table("drug_target.txt",sep="\t",header=T,check.names=F)

kegg <- read.table("KEGG_NEED.txt",sep="\t",header=T,check.names=F)

samples <- intersect(kegg$name,drug_target$gene_name)

kegg <- kegg[which(kegg$name %in% samples),]
drug_target <- drug_target[which(drug_target$gene_name %in% samples),]
unique(drug_target$drug_name)

  drug_kegg <- NULL
for (i in 1:nrow(drug_target)) {
  
  kegg_need <- kegg[which(kegg$name == drug_target[i,"gene_name"]),]
  kegg_need$name <- drug_target[i,"drug_name"]
  
  drug_kegg <- rbind(kegg_need,drug_kegg)
}

write.table(drug_kegg, "drug_kegg_input.txt", sep="\t", quote=F, row.names=T) 
  
  
  