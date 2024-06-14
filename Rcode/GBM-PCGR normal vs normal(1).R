# 加载所需的包
library(DESeq2)
library(ggplot2)
library(ggthemes)

# 读取两组癌症细胞表达矩阵数据
expr_matrix1 <- read.table("/Users/xiaoqiluo/Desktop/生信实验大作业/TCGA-GBM.htseq_counts.tsv", header=TRUE,sep="")
expr_matrix2 <- read.table("/Users/xiaoqiluo/Desktop/生信实验大作业/TCGA-PCPG.htseq_counts.tsv", header=TRUE,sep="")

# 筛选正常细胞
normal_1 <- grepl("11A", colnames(expr_matrix1))
normal_2 <- grepl("11A", colnames(expr_matrix2))
normal_matrix1 <- expr_matrix1[, normal_1]
normal_matrix2 <- expr_matrix2[, normal_2]

#数据处理
normal_matrix1<-(2**normal_matrix1)-1
normal_matrix1<-round(normal_matrix1)
normal_matrix2<-(2**normal_matrix2)-1
normal_matrix2<-round(normal_matrix2)
nA <- ncol(normal_matrix1)
nB <- ncol(normal_matrix2)

#创建用于差异分析的矩阵
con12<-cbind(normal_matrix1,normal_matrix2)
group12 <- factor(c(rep(1, nA), rep(2, nB)))
colData12 <- data.frame(Sample=colnames(con12), Group=group12)

# 差异分析
dds12 <- DESeqDataSetFromMatrix(countData=con12, colData=colData12, design=~Group)
keep <- rowSums(counts(dds12)) >= 10
dds12 <- dds12[keep,]
dds12$Group <- factor(dds12$Group)
dds12 <- DESeq(dds12)
res12 <- results(dds12)
res12$padj[is.na(res12$padj)] <- 1

# 保存结果
res_ordered12 <- res12[order(res12$padj),]
write.csv(res_ordered12, "GBM-PCPG normal.csv")
head(res12)

#作火山图
volcano12 <- data.frame('log2FC' = res12$log2FoldChange, '-log10Padj' = log10(res12$padj),'pvalue'=res12$pvalue)
volcano12$Sig <- ifelse(volcano12$pvalue < 0.0001&abs(volcano12$log2FC)>=2,ifelse(volcano12$log2FC >2,'Up','Down'),'None')
ggplot(volcano12, aes(x = log2FC, y = -X.log10Padj,colour = Sig)) +
  geom_point(alpha = 0.5) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value")
ggsave("GBM-PCPG normal.png",width = 10, height = 6, dpi = 300)



# 合并两个矩阵
normal_matrix <- cbind(normal_matrix1, normal_matrix2)
# 创建分组信息
group <- factor(c(rep(1, nA), rep(2, nB)))
colData <- data.frame(Sample = colnames(normal_matrix), Group = group)
# 创建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = normal_matrix, colData = colData, design = ~ Group)
# 过滤低表达基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# 设置分组因子
dds$Group <- factor(dds$Group)
# 运行 DESeq 分析
dds <- DESeq(dds)
res <- results(dds)
res$padj[is.na(res$padj)] <- 1
# 进行变换
vsd <- vst(dds, blind = TRUE)
# 创建 PCA 图形对象
pca_plot <- plotPCA(vsd, intgroup = "Group") +
  geom_point(size = 0.001)
  ggtitle(label = "Principal Component Analysis (PCA)")
# 保存 PCA 图形为高分辨率 PNG 文件
ggsave("PCA_Plot normal vs normal.png", plot = pca_plot, width = 10, height = 6, dpi = 300)

library(ggrepel)
# 进行K-means聚类
vsd <- vst(dds12, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
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
ggsave("PCA_Kmeans_Clustering_normal.png", plot = pca_plot, width = 10, height = 6, dpi = 300)







