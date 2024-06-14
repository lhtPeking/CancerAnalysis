library(tidyverse)
library(DESeq2)
library("BiocParallel")
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ggthemes)

#清空变量
rm(list=ls(all=TRUE))
options(stringsAsFactors = F)

# setwd("/Users/xiaoqiluo/Desktop/生信实验大作业")

##读取
GBM <-read.csv('bioinfoHW/DataFiles/TCGA-GBM.htseq_counts.tsv', header = T, sep = "") 
rownames(GBM) <- GBM[,1] 
GBM <- GBM[,-1]

PCPG <-read.csv('bioinfoHW/DataFiles/TCGA-PCPG.htseq_counts.tsv', header = T, sep = "") 
rownames(PCPG) <- PCPG[,1] 
PCPG <- PCPG[,-1]


##保留矩阵中癌细胞的样本
GBM <- GBM[,as.numeric(str_sub(colnames(GBM), 14, 15)) < 10]
PCPG <- PCPG[,as.numeric(str_sub(colnames(PCPG), 14, 15)) < 10]


rawcount <- cbind(GBM,PCPG) %>% as.matrix() #合并转化为表达矩阵
rawcount <- (2**rawcount)-1 #获得原始reads数


#设置样品信息矩阵
condition <- factor(c(rep("GBM",172),rep("PCPG",179)),levels = c("GBM","PCPG"))
colData <- data.frame(row.names = colnames(rawcount),condition)

#用DESeq2进行差异表达分析
dds <- DESeqDataSetFromMatrix(round(rawcount),colData,design = ~condition)
dds_remove <- DESeq(dds[rowSums(counts(dds)) >= 10 , ]) #保留表达量加起来大于等于10的行

#两两比较
register(SnowParam(4))#4线程加快运行速度
find_sig_DEGs = function(contra){
  res = results(dds_remove,contrast = contra)
  DEG = res[order(res$padj),] %>% na.omit()
  sig_DEGs = DEG[DEG$padj < 0.01  ,] 
  return(sig_DEGs)
} ##函数，根据p-value 升序，去除NA行，选出pvalue小于0.01的基因

GBMvsPCPG_contrast<- c("condition","GBM","PCPG")
GBMvsPCPG_DEG <- find_sig_DEGs(GBMvsPCPG_contrast)

#PCA
vsd <- vst(dds, blind = TRUE)
# 创建 PCA 图形对象
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  geom_point(size = 0.001) +
  ggtitle(label = "Principal Component Analysis (PCA)")
# 保存 PCA 图形为高分辨率 PNG 文件
ggsave("PCA_Plot.png", plot = pca_plot, width = 10, height = 6, dpi = 300)


# 进行K-means聚类
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
# 提取前两个主成分
pcaData <- as.data.frame(pcaData)
pcaData <- pcaData[, 1:2] # 提取前两个主成分
set.seed(123) # 设置随机种子以获得可重复的结果
kmeans_result <- kmeans(pcaData, centers = 2) # 聚成2类（根据需求调整centers的值）
pcaData$cluster <- as.factor(kmeans_result$cluster)
# 创建 PCA 图形对象，并根据K-means聚类结果进行标注
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = rownames(pcaData)), size = 3) +
  ggtitle(label = "PCA with K-means Clustering")
ggsave("PCA_Kmeans_Clustering.png", plot = pca_plot, width = 10, height = 6, dpi = 300)


#统计差异表达基因的信息
get_statistics = function(contra){
  filter_up <- subset(contra,log2FoldChange >1)#上调基因
  filter_down <- subset(contra,log2FoldChange < -1)#下调基因
  print(paste("差异上调基因数量：",nrow(filter_up)))#打印上调基因数量
  print(paste("差异下调基因数量：",nrow(filter_down)))#打印下调基因数量
  a <- substring(deparse(substitute(contra)),1,10)
  #保存结果
  write.table(contra, file = paste(a,"differential_gene.txt",sep = "_"),quote = F)
  write.table(filter_up,file = paste(a,"filter_up_gene.txt",sep = "_"),quote = F)
  write.table(filter_down,file = paste(a,"filter_down_gene.txt",sep = "_"),quote = F)
}
get_statistics(GBMvsPCPG_DEG)
head(GBMvsPCPG_DEG)

library(ggplot2)
library(ggrepel)

# 基因表达分类
GBMvsPCPG_DEG <- as.data.frame(GBMvsPCPG_DEG)
GBMvsPCPG_DEG$group <- ifelse(GBMvsPCPG_DEG$log2FoldChange >= 1 & GBMvsPCPG_DEG$padj <= 0.01, "Up",
                              ifelse(GBMvsPCPG_DEG$log2FoldChange <= -1 & GBMvsPCPG_DEG$padj <= 0.01, "Down", "Not sig"))

# 查看分类情况
table(GBMvsPCPG_DEG$group)

# 挑选一些明显表达上调或者下调的基因
GBMvsPCPG_DEG$label <- ifelse((GBMvsPCPG_DEG$padj < 1.6e-60 & GBMvsPCPG_DEG$log2FoldChange <= -12) | 
                                (GBMvsPCPG_DEG$padj < 1.6e-90 & GBMvsPCPG_DEG$log2FoldChange >= 6.5),
                              as.character(row.names(GBMvsPCPG_DEG)), '')

# 绘制火山图
ggplot(GBMvsPCPG_DEG, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = group), size = 0.5) +  # 调整点的大小，size=1 表示点的大小较小
  scale_color_manual(values = c("#eb7979", "grey", "#4a86e8"), limits = c("Up", "Not sig", "Down")) +
  theme_bw(base_size = 20) + 
  ggtitle("GBMvsPCPG Volcano Plot") +
  theme(plot.title = element_text(size = 25, hjust = 0.5)) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0, 400)) +
  geom_hline(yintercept = 2, linetype = 2) +   # 添加虚线
  geom_vline(xintercept = c(-1, 1), linetype = 2)  # 添加虚线
ggsave("GBMvsPCPG Volcano Plot_new.png",width = 10, height = 6, dpi = 300)
