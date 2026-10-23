
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/TCGA/TCGA_SKCM_TPM")

list.files()

library(data.table)
library(tidyverse)

color1 = c("#53049d", "#7500aa", "#b17dd5", "#c587c7", 
           "#a0259f", "#d16ba5", "#dea7bc", "#f7b3ae")

color2 = c("#0f0a80", "#4b1475", "#6400a6", "#870fa1",
           "#aa2494", "#da596c", "#eb7954", "#f0bb3f")

color3 = c("#573299", "#6e3fc1", "#8a63cd", "#a78cd9", "#c4b3e8")

color4 = c("#1d136d", "#421b7b", "#600e86", "#9d2a7b",
           "#c44a60", "#ca5662", "#e5a544", "#fabdc0")

color5 = c("#360c6b", "#59187a", "#a3307a", "#d4466c", "#f46a57")


tpm <- fread("TCGA-SKCM.star_tpm.tsv.gz.tsv", header = T, sep = "\t", check.names = F)
clinical <- fread("TCGA-SKCM.clinical.tsv.gz.tsv", header = T, sep = "\t", check.names = F)
survival <- fread("TCGA-SKCM.survival.tsv.gz.tsv", header = T, sep = "\t", check.names = F)
genesymbol <- fread("gencode.v36.annotation.gtf.gene.probemap", header = T, sep = "\t", check.names = F)

head(genesymbol)
tpm[1:5, 1:5]

tpm_with_symbol <- merge(tpm, genesymbol[, .(id, gene)], 
                         by.x = "Ensembl_ID", 
                         by.y = "id", 
                         all.x = TRUE)
tpm_with_symbol <- distinct(.data = tpm_with_symbol, gene, .keep_all = T)     # 根据gene列去除重复基因
tpm_with_symbol <- column_to_rownames(.data = tpm_with_symbol, var = "gene")   
tpm_with_symbol <- tpm_with_symbol[,-1]

colnames(tpm_with_symbol)    # 473样本
clinical$sample    # 476样本
survival$sample    # 462样本
colnames(tpm_with_symbol) %in% clinical$sample 
survival$sample %in% colnames(tpm_with_symbol)

sample_recorded <- intersect(survival$sample, colnames(tpm_with_symbol))

tpm_recorded <- tpm_with_symbol[,sample_recorded]     # 有survival记录的sample



library(survival)
library(survminer)
library(dplyr)
library(ggplot2)


gene_name <- "IGSF3"
exp_vec <- tpm_recorded[gene_name, ]
exp_df <- data.frame(sample = names(exp_vec), 
                     Exp_TPM = as.numeric(exp_vec), 
                     stringsAsFactors = FALSE)

surv_df <- survival %>%
  left_join(exp_df, by = "sample") %>%
  filter(!is.na(OS.time) & !is.na(OS) & !is.na(Exp_TPM))

med_val <- median(surv_df$Exp_TPM)
surv_df$group <- ifelse(surv_df$Exp_TPM > med_val, "High", "Low")

n_high <- sum(surv_df$group == "High")
n_low <- sum(surv_df$group == "Low")

fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_df)

# --- 1. 设置统一的字体参数 ---
# 建议 14 或 12，这个值会同时作用于 Legend 和 P值
base_font_size <- 14  
font_family <- "Arial" 

# --- 2. 运行绘图逻辑 ---
splot <- ggsurvplot(
  fit,
  data = surv_df,
  palette = c("#f46a57", "#870fa1"), 
  linetype = "dashed",               
  censoring.shape = 3,                
  censor.size = 2,
  xlim = c(-500, 12000),                
  break.time.by = 2500,     
  ylim = c(-0.05, 1.05),
  legend = "none",  # 初始设为 none，后续通过 splot$plot 覆盖
  xlab = "Days",
  ylab = "Probability of survival",
  axes.offset = FALSE,               
  legend.labs = c(paste0(gene_name, " high (n=", n_high, ")"), 
                  paste0(gene_name, " low (n=", n_low, ")")),
  ggtheme = theme_bw() + 
    theme(
      # 设置全局字体和基础大小
      text = element_text(family = font_family, size = base_font_size),
      
      # 坐标轴与网格设置
      axis.line = element_line(linewidth = 0, color = "black"), # 建议加回 axis.line 更有出版质感
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length = unit(-0.15, "cm"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(linewidth = 0.8), # 保持方框风格
      
      axis.text.x = element_text(margin = margin(t = 10), family = font_family, size = base_font_size),
      axis.text.y = element_text(margin = margin(r = 10), family = font_family, size = base_font_size),
      
      # 图例属性设置 (位置在后续步骤强制覆盖)
      legend.text = element_text(family = font_family, size = base_font_size),
      legend.title = element_blank(),
      legend.background = element_blank(), 
      legend.key = element_blank()
    )
)

# --- 3. 准备 P 值文字 ---
logrank_p <- surv_pvalue(fit)$pval
cox_p <- summary(coxph(Surv(OS.time, OS) ~ group, data = surv_df))$waldtest["pvalue"]
p_text <- paste0("LogRank p = ", sprintf("%.3f", logrank_p), "\n",
                 "Cox p = ", sprintf("%.3f", cox_p))

# --- 4. 关键修正：统一文字大小并叠加 ---

# 计算 annotate 所需的 size (mm)
# ggplot2 的 annotate size 使用 mm，而 theme 使用 pt。 1 mm ≈ 2.845 pt
annotate_size <- base_font_size / 2.845

splot$plot <- splot$plot + 
  annotate("text", 
           x = 500, 
           y = 0.15, 
           label = p_text, 
           family = font_family, 
           size = annotate_size, # 物理大小现在与 Legend 一致
           hjust = 0, 
           lineheight = 1.1) +
  # 强制图例位置 (对应你代码里最后的设置)
  theme(
    legend.position = c(0.7, 0.9),
    legend.text = element_text(family = font_family, size = base_font_size)
  )

# 显示结果
print(splot)

ggsave("IGSF3_skcm_SKCM.svg", width = 4.25, height = 4.05, dpi = 300)


### PD-L1 IGSF3 correlate
igsf3 <- as.numeric(tpm_with_symbol["IGSF3", ])   # 提取两个基因的表达
cd274 <- as.numeric(tpm_with_symbol["CD274", ])
cor.test(igsf3, cd274, method = "spearman")     # Spearman 相关性

# 提取两个基因的表达
plot_df <- data.frame(
  sample = colnames(tpm_with_symbol),
  IGSF3 = as.numeric(tpm_with_symbol["IGSF3", ]),
  CD274 = as.numeric(tpm_with_symbol["CD274", ])
)

# 计算 Spearman 相关性
spearman_result <- cor.test(plot_df$IGSF3, plot_df$CD274, method = "spearman")
pearson_result <- cor.test(plot_df$IGSF3, plot_df$CD274, method = "pearson")

p <- ggplot(plot_df, aes(x = IGSF3, y = CD274)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = Inf, y = Inf,
           label = paste0(
             "Pearson R = ",
             round(pearson_result$estimate, 3),
             "\nP = ",
             signif(pearson_result$p.value, 3)
           ),
           hjust = 1.1, vjust = 1.5, size = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(x = "IGSF3 expression", y = "CD274 expression")
p
# ggsave("CD274_vs_IGSF3_skcm_tumortissues.png", height = 5, width = 5, dpi = 300)

