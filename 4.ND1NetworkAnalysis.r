# ND1 Network Analysis

library(GENIE3)
library(parallel)

# ChIP Peak & RNA coorelation
peak <- "ND1_ChIP.peaks.anno"
possibleTargets <- peak$`Gene Name`

exprMatr <- as.matrix(sc.vitro@assays$RNA@data)
possibleTargets <- intersect(rownames(exprMatr), possibleTargets)
exprMatr <- exprMatr[possibleTargets,]
exprMatr <- exprMatr[which(rowSums(exprMatr) > 1),]

DEG <- read.csv("VitroMarker.csv")
tfdb <- read.delim("Rattus_norvegicus_TF.txt")
DEG <- subset(DEG, gene %in% tfdb$Symbol)
NeuDEG.df <- subset(NeuDEG.df, p_val_adj < 0.0001 & abs(avg_log2FC) > 0.5)
potential_regulatores <- c(unique(NeuDEG.df$gene), "Meis2")
potential_regulatores <- potential_regulatores[which(potential_regulatores %in% possibleTargets)]

regulators <- potential_regulatores
weightMat <- GENIE3(exprMatr, regulators=regulators, nCores = 15)

linkList <- getLinkList(weightMat)

for (i in 1:nrow(linkList)) {
  gene1 <- linkList[i, "regulatoryGene"]
  gene2 <- linkList[i, "targetGene"]
  res <- cor.test(exprMatr[gene1,], exprMatr[gene2,])
  linkList[i,"cor"] <- res$estimate
  linkList[i,"pvalue"] <- res$p.value
  print(i)
}
linkList$padj <- p.adjust(linkList$pvalue, method = "BH")

write.csv(linkList, file = file.path(FIG4, "ND1_Gene_NetWork_20240413_addcor.csv"))