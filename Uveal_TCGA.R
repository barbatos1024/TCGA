
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/TCGA/Uveal_TCGA")

list.files()

library(data.table)

tpm <- fread("tcga_target_no_normal_rsem_gene_tpm.gz", header = T, sep = "\t")
phenotype <- read.table("TCGA_TARGET_phenotype", header = T, sep = "\t", check.names = F)
probe <- read.table("probeMap-gencode.v23.annotation.gene.probemap", header = T, sep = "\t", check.names = F)

tpm[1:5, 1:5]
probe[1:5, 1:5]

tpm$sample <- gsub("\\..*", "", tpm$sample)
probe$id <- gsub("\\..*", "", probe$id)

library(tidyverse)

tpm_with_symbol <- tpm %>%
  inner_join(probe %>% select(id, gene), by = c("sample" = "id")) %>%
  select(gene, everything(), -sample)
tpm_with_symbol <- distinct(.data = tpm_with_symbol, gene, .keep_all = T)
tpm_with_symbol <- column_to_rownames(.data = tpm_with_symbol, var = "gene")

table(phenotype$`_cohort`)

library(dplyr)

mela_pheno <- phenotype %>%
  filter(`_cohort` %in% c("TCGA Melanoma", "TCGA Ocular melanomas"))

mela_samples <- mela_pheno$sampleID
tpm_samples <- colnames(tpm_with_symbol)

common_samples <- intersect(mela_samples, tpm_samples)

mela_tpm <- tpm_with_symbol[, common_samples]
mela_pheno_final <- mela_pheno[match(common_samples, mela_pheno$sampleID), ]     # 更新临床信息表，保持样本一一对应且顺序一致
table(mela_pheno_final$`_primary_site`)

library(limma)

group_list = c(rep('Skin', 469), rep("Eye", 79))
group_list = factor(group_list)
group_list <- relevel(group_list, ref = "Skin")

data = normalizeBetweenArrays(mela_tpm)
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-13, 18))
boxplot(mela_tpm, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-13, 18))


library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

gene_name <- "IGSF3"

if (gene_name %in% rownames(data)) {
  plot_data <- data.frame(
    Sample = colnames(data),
    Expression = as.numeric(data[gene_name, ]),
    Group = c(rep('Skin', 469), rep("Eye", 79)) 
  )
} else {
  stop(paste0("找不到", gene_name))
}

plot_data$Group <- factor(plot_data$Group, levels = c("Skin", "Eye"))


my_colors <- c("Skin" = "#2E31F3", "Eye" = "#FF1E1E")

ggplot(plot_data, aes(x = Group, y = Expression, color = Group)) +
  # 1. 绘制背景箱线图 (虚线边框)
  geom_boxplot(aes(fill = Group), 
               outlier.shape = NA, 
               alpha = 0, 
               linetype = "dashed", 
               color = "black", 
               width = 0.5) +
  # 2. 叠加中位数实线
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.4, color = "black", size = 0.5) +
  # 3. 叠加散点 (因为点少，dodge.width 设为 0 居中显示)
  geom_jitter(size = 4, width = 0.1, alpha = 0.9) +
  # 4. 颜色与主题
  scale_color_manual(values = my_colors) +
  theme_bw() + 
  ylim(0,10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    axis.title = element_text(size = 13),
    legend.position = "none" # 点少，直接看横轴即可，不需要图例
  ) +
  stat_compare_means(method = "t.test", 
                     label = "p.format",      # 显示具体 P 值，如 p = 0.0012
                     label.x = 1.35,           # 标注显示的横坐标位置（1.5 位于两组中间）
                     label.y = max(plot_data$Expression) * 1.05, # 标注显示的纵坐标高度
                     comparisons = list(c("Skin", "Eye")), # 指定比较对象画连线
                     bracket.size = 0.4,
                     size = 5,                # 字体大小
                     color = "black") +       # 字体颜色
  labs(
    title = "TCGA_Melanoma",
    y = paste0(gene_name, " Expression log2 (TPM + 0.001)"),
    x = ""
  )

# ggsave(filename = "Melanoma_Skin_Eye.png", width = 3, height = 5.5)


ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
  # 1. 绘制小提琴图，draw_quantiles 会画出分位数线（0.5是中位数）
  geom_violin(trim = TRUE, 
              color = "black", 
              linewidth = 0.6,
              draw_quantiles = c(0.25, 0.5, 0.75), 
              linetype = "dashed") + 
  
  # 2. 覆盖中位数线，使其变为实线（对应图中效果）
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.2, color = "black", linewidth = 0.7) +
  
  # 3. 添加显著性连线 (比之前更接近图中风格)
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     comparisons = list(c("Skin", "Eye")), # 指定比较对象画连线
                     bracket.size = 0.4,
                     size = 5) +
  
  # 4. 颜色与坐标轴精修
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # 让底部从0开始，顶部留点空间
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    legend.position = "none" # 点少，直接看横轴即可，不需要图例
  ) +
  labs(
    y = paste0(gene_name, " Expression log2 (TPM + 0.001)"))
  )





