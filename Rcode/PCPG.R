install.packages("BiocManager")

# 使用 BiocManager 安装 GenomicFeatures 和 AnnotationDbi 包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 使用 BiocManager 安装 GenomicFeatures 和 AnnotationDbi 包
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# 安装 EnhancedVolcano 包
BiocManager::install("EnhancedVolcano")

BiocManager::install("DESeq2")
library(DESeq2)

cts <- read.csv("/Users/xiaoqiluo/Desktop/生信实验大作业/TCGA-PCPG.htseq_counts.tsv",sep="\t",row.names="Ensembl_ID", header = TRUE)
cts <- 2^cts - 1
cts <- round(cts)

thecolnames <- colnames(cts)
thelabel <- vector()
for (var in thecolnames){
  if (substr(var,14,14) == "0"){
    thelabel <- append(thelabel,"tumor")
  } else {
    thelabel <- append(thelabel,"normal")
  }
}
trymatrix <- cbind(thecolnames,thelabel)
rownames(trymatrix)=trymatrix[,1]  #取出第一列
trymatrix=trymatrix[,-1]          #将第一列删除
trymatrix_new <- as.matrix(trymatrix)
colnames(trymatrix_new)<-c("condition")
coldata <- trymatrix_new

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("normal","tumor"))
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","tumor","normal")) 

resOrdered <- res[order(res$pvalue),]
summary(res)
write.csv(as.data.frame(resOrdered), file = "PCPG_condition_treated_results.csv")

## 火山图绘制
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 安装 EnhancedVolcano 包
install.packages("ggthemes")
library(ggthemes)
EnhancedVolcano(resOrdered,
                x = 'log2FoldChange',
                y = 'pvalue',
                lab = rownames(resOrdered),
                title = 'PCPG Cancer versus Normal',
                xlim=c(-10,15),
                ylim=c(-10,80),
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize=1.5,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 5,
                drawConnectors=TRUE,
                widthConnectors=0.5)
ggsave("PCPG_new.png",width = 10, height = 6, dpi = 300)