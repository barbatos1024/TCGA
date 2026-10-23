
rm(list = ls())

###设置路径
getwd()
setwd("/Users/naturexv/Datasets/TCGA/Isoform_TCGA")

###加载包  这个包是用来加快数据读取的
library(data.table)
library(dplyr)
library(tidyverse)

### 读取基因信息
geneid = fread("probeMap-gencode.v23.annotation.transcript.probemap", header = T, sep = '\t', data.table = F)

### 读取表达rnaseq的数据
mem.maxVSize(vsize = 1e19)
tpm = fread("TcgaTargetGtex_rsem_isoform_tpm.gz", header = T, sep = '\t', data.table = F)
tpm <- column_to_rownames(tpm, var = "sample")
tpm[1:5, 1:5]

### 读取临床信息
phe = fread("TcgaTargetGTEX_phenotype.txt.gz", header = T, sep = '\t', data.table = F, check.names = F)
colnames(phe) <- c("sample", "detailedcategory", "primarydiseasetissue", "primarysite", "sampletype", "gender", "study")
skin <- filter(.data = phe, primarysite == "Skin")
table(skin$primarydiseasetissue)
skin.1 <- skin[c(2:470,727:1282),]   # 去掉癌旁组织 和 细胞系
table(skin.1$primarydiseasetissue)
table(skin.1$sampletype)
skin.1$ST <- c(rep("T", 102), rep("M", 367), rep("N", 556))

skcm <- intersect(skin.1$sample, colnames(tpm))
skcm.data <- tpm[,skcm]

rm(tpm)

library(limma)
normalized.data = normalizeBetweenArrays(skcm.data)
normalized.data[1:5, 1:5]

group_list = c(rep("T", 102), rep("M", 367), rep("N", 556))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "N")

#去除批次效应后
boxplot(normalized.data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-20,18))
#去除批次效应前
boxplot(skcm.data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-20,13))

id.igsf3 <- filter(.data = geneid, gene == "IGSF3")
id.igsf3$id

##########################################################################
"ENST00000369483.5" "ENST00000369486.7" "ENST00000318837.6" "ENST00000481589.1"

################################################# 箱线图
library(ggpubr)
normalized.data <- as.data.frame(normalized.data)
GENE = normalized.data["ENST00000481589.1", ]
z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框

x = c(rep("T", 102), rep("M", 367), rep("N", 556))

z$Type = x

z$Type <- factor(z$Type, levels = c("N", "T", "M"))

colnames(z)

p1 = ggboxplot(z, x = "Type", y = "ENST00000481589.1",
               width = 0.6, fill = "Type",
               notch = F, palette = c("grey", "red", "blue", "yellow"),
               add = "jitter", shape = "Type")

p1
###加个p值

p2 = ggboxplot(z, x = "Type", y = "ENST00000481589.1",
               width = 0.6, fill = "Type",
               notch = F, palette = c("grey", "red", "blue", "black", "green", "yellow"),
               add = "jitter", shape = "Type") + 
  labs(x = "", y = "Gene Expression log2 (TPM + 0.001)", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 0, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c("N", "T"), c("N", "M"), c("T", "M")), 
                     label.x = 1.7, label.y = c(8.5,10,11.5), method = "wilcox.test") +   ### t.test wilcox.test
  ggtitle("IGSF3") 

p2

ggsave("igsf3.ENST00000481589.1.png", width = 4, height = 5, dpi = 300)

################# isoform之间比较 ###############
normalized.data <- as.data.frame(normalized.data)
tpm.igsf3 <- normalized.data %>% filter(rownames(normalized.data) %in% id.igsf3$id)
z = t(tpm.igsf3) #转置
z = as.data.frame(z) #转化为数据框

x = c(rep("T", 102), rep("M", 367), rep("N", 556))

z$Type = x

z$Type <- factor(z$Type, levels = c("N", "T", "M"))

colnames(z)

library(tidyverse)

# 转换为long format，便于 ggboxplot 绘图
df_long <- z %>%
  pivot_longer(
    cols = starts_with("ENST"),  # 或者你手动指定四个列名
    names_to = "Transcript",
    values_to = "Expression"
  )

ggboxplot(df_long,
          x = "Transcript",     # x轴为转录本
          y = "Expression",     # y轴为表达值
          color = "Type",       # 颜色按样本类型
          add = "jitter",       # 加入散点
          palette = "lancet") +    # 颜色主题可换，如 "npg", "lancet", "jco"
  theme_classic2() +
  labs(x = "", y = "Gene Expression log2 (TPM + 0.001)", color = "") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "grey40") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "top") +
  ggtitle("IGSF3 isoform") 

ggsave("igsf3.isoform.png", width = 7, height = 7, dpi = 300)

######
# 创建一个新的数据框，仅包含两个转录本的表达值
df_plot <- data.frame(
  IGSF3_v1 = z$`ENST00000369483.5`,
  IGSF3_v2 = z$`ENST00000369486.7`,
  Type = z$Type
)

# 画散点图并加上相关系数
ggscatter(df_plot, 
          x = "IGSF3_v1", 
          y = "IGSF3_v2",
          color = "Type", 
          add = "reg.line",              # 添加回归线
          conf.int = TRUE,               # 添加置信区间
          cor.coef = TRUE,               # 显示相关系数
          cor.method = "pearson",        # 使用 Pearson 方法
          xlab = "ENST00000369483.5(IGSF3_v1)", 
          ylab = "ENST00000369486.7(IGSF3_v2)",
          palette = "npg",
          add.params = list(color = "red"))

ggsave("igsf3.correlation.png", width = 7, height = 7, dpi = 300)



