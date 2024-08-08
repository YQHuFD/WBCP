library(data.table)
library(dplyr)
library(lsa)
setwd("E:\\药物联用项目\\操作\\相似性计算\\kegg/")
# Read drugbank_id -> side_effects file
db2se <- read.table("drug_kegg_input.txt",sep="\t",header=T,check.names=F)

# Create empty matrix
sim <- matrix(data = 0, nrow = length(unique(db2se$drug)), ncol = length(unique(db2se$event)))
rownames(sim) <- unique(db2se$drug)
colnames(sim) <- unique(db2se$event)

# Count 1s and 0s
for (i in 1:nrow(sim)) {
  sel <- db2se %>% filter(drug == rownames(sim)[i]) %>% select(event)
  sim[i, ] <- as.numeric(colnames(sim) %in% sel$event)
}

# Compute IDF
tf <- sim
idf <- log(nrow(sim) / colSums(sim))
tfidf <- sim

# Compute term TF * IDF
for(word in names(idf)){
  tfidf[, word] <- tf[, word] * idf[word]
}

# Normalize matrix
tfidf_norm <- tfidf / sqrt(rowSums(tfidf^2))

# Compute similarity matrix
cos <- cosine(as.matrix(t(tfidf_norm)))

# Write table
write.table(cos, "kegg_similarity.txt", sep="\t", quote=F, row.names=T)

