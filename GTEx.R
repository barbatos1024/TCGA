
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/TCGA/GTEx/")

list.files()

library(data.table)
library(dplyr)
library(tidyverse)
library(ggpubr)

color = c("#53049d", "#7500aa", "#b17dd5", "#c587c7", "#a0259f", "#d16ba5", "#dea7bc", "#f7b3ae",
          "#0f0a80", "#4b1475", "#6400a6", "#870fa1", "#aa2494", "#da596c", "#eb7954", "#f0bb3f",
          "#573299", "#6e3fc1", "#8a63cd", "#a78cd9", "#c4b3e8", "#1d136d", "#421b7b", "#600e86",
          "#9d2a7b", "#c44a60", "#ca5662", "#e5a544", "#fabdc0", "#360c6b", "#59187a", "#a3307a",
          "#d4466c", "#f46a57")


tpm <- fread("gtex_RSEM_gene_tpm.gz", header = T, sep = "\t")
annotation <- fread("probeMap-gencode.v23.annotation.gene.probemap")
metadata <- read.table("GTEX_phenotype.gz", header = T, sep = "\t", row.names = 1)

setnames(tpm, 1, "id")
tpm_with_gene <- merge(tpm, annotation[, .(id, gene)], by = "id", all.x = TRUE)

tpm_with_gene <- distinct(.data = tpm_with_gene, gene, .keep_all = T)
tpm_with_gene <- column_to_rownames(.data = tpm_with_gene, var = "gene")
tpm_with_gene <- tpm_with_gene[,-1]

metadata[1:5, 1:5]
tpm_with_gene[1:5, 1:5]


colnames(metadata)
setnames(metadata, 1, "body_site_detail")
setnames(metadata, 2, "primary_site")
setnames(metadata, 3, "gender")
setnames(metadata, 4, "patient")
setnames(metadata, 5, "cohort")

table(metadata$primary_site)
Adipose_Tissue = filter(metadata, primary_site == "Adipose Tissue")
Adrenal_Gland = filter(metadata, primary_site == "Adrenal Gland")
Bladder = filter(metadata, primary_site == "Bladder")
Blood = filter(metadata, primary_site == "Blood")
Blood_Vessel = filter(metadata, primary_site == "Blood Vessel")
Bone_Marrow = filter(metadata, primary_site == "Bone Marrow")
Brain = filter(metadata, primary_site == "Brain")
Breast = filter(metadata, primary_site == "Breast")
Cervix_Uteri = filter(metadata, primary_site == "Cervix Uteri")
Colon = filter(metadata, primary_site == "Colon")
Esophagus = filter(metadata, primary_site == "Esophagus")
Fallopian_Tube = filter(metadata, primary_site == "Fallopian Tube")
Heart = filter(metadata, primary_site == "Heart")
Kidney = filter(metadata, primary_site == "Kidney")
Liver = filter(metadata, primary_site == "Liver")
Lung = filter(metadata, primary_site == "Lung")
Muscle = filter(metadata, primary_site == "Muscle")
Nerve = filter(metadata, primary_site == "Nerve")
Ovary = filter(metadata, primary_site == "Ovary")
Pancreas = filter(metadata, primary_site == "Pancreas")
Pituitary = filter(metadata, primary_site == "Pituitary")
Prostate = filter(metadata, primary_site == "Prostate")
Salivary_Gland = filter(metadata, primary_site == "Salivary Gland")
Skin = filter(metadata, primary_site == "Skin")
Small_Intestine = filter(metadata, primary_site == "Small Intestine")
Spleen = filter(metadata, primary_site == "Spleen")
Stomach = filter(metadata, primary_site == "Stomach")
Testis = filter(metadata, primary_site == "Testis")
Thyroid = filter(metadata, primary_site == "Thyroid")
Uterus = filter(metadata, primary_site == "Uterus")
Vagina = filter(metadata, primary_site == "Vagina")

z1 = intersect(rownames(Adipose_Tissue), colnames(tpm_with_gene))
z2 = intersect(rownames(Adrenal_Gland), colnames(tpm_with_gene))
z3 = intersect(rownames(Bladder), colnames(tpm_with_gene))
z4 = intersect(rownames(Blood), colnames(tpm_with_gene))
z5 = intersect(rownames(Blood_Vessel), colnames(tpm_with_gene))
z6 = intersect(rownames(Bone_Marrow), colnames(tpm_with_gene))
z7 = intersect(rownames(Brain), colnames(tpm_with_gene))
z8 = intersect(rownames(Breast), colnames(tpm_with_gene))
z9 = intersect(rownames(Cervix_Uteri), colnames(tpm_with_gene))
z10 = intersect(rownames(Colon), colnames(tpm_with_gene))

z11 = intersect(rownames(Esophagus), colnames(tpm_with_gene))
z12 = intersect(rownames(Fallopian_Tube), colnames(tpm_with_gene))
z13 = intersect(rownames(Heart), colnames(tpm_with_gene))
z14 = intersect(rownames(Kidney), colnames(tpm_with_gene))
z15 = intersect(rownames(Liver), colnames(tpm_with_gene))
z16 = intersect(rownames(Lung), colnames(tpm_with_gene))
z17 = intersect(rownames(Muscle), colnames(tpm_with_gene))
z18 = intersect(rownames(Nerve), colnames(tpm_with_gene))
z19 = intersect(rownames(Ovary), colnames(tpm_with_gene))
z20 = intersect(rownames(Pancreas), colnames(tpm_with_gene))

z21 = intersect(rownames(Pituitary), colnames(tpm_with_gene))
z22 = intersect(rownames(Prostate), colnames(tpm_with_gene))
z23 = intersect(rownames(Salivary_Gland), colnames(tpm_with_gene))
z24 = intersect(rownames(Skin), colnames(tpm_with_gene))
z25 = intersect(rownames(Small_Intestine), colnames(tpm_with_gene))
z26 = intersect(rownames(Spleen), colnames(tpm_with_gene))
z27 = intersect(rownames(Stomach), colnames(tpm_with_gene))
z28 = intersect(rownames(Testis), colnames(tpm_with_gene))
z29 = intersect(rownames(Thyroid), colnames(tpm_with_gene))
z30 = intersect(rownames(Uterus), colnames(tpm_with_gene))
z31 = intersect(rownames(Vagina), colnames(tpm_with_gene))


Adipose_Tissue = tpm_with_gene[,z1]
Adrenal_Gland = tpm_with_gene[,z2]
Bladder = tpm_with_gene[,z3]
Blood = tpm_with_gene[,z4]
Blood_Vessel = tpm_with_gene[,z5]
Bone_Marrow = tpm_with_gene[,z6]
Brain = tpm_with_gene[,z7]
Breast = tpm_with_gene[,z8]
Cervix_Uteri = tpm_with_gene[,z9]
Colon = tpm_with_gene[,z10]
Esophagus = tpm_with_gene[,z11]
Fallopian_Tube = tpm_with_gene[,z12]
Heart = tpm_with_gene[,z13]
Kidney = tpm_with_gene[,z14]
Liver = tpm_with_gene[,z15]
Lung = tpm_with_gene[,z16]
Muscle = tpm_with_gene[,z17]
Nerve = tpm_with_gene[,z18]
Ovary = tpm_with_gene[,z19]
Pancreas = tpm_with_gene[,z20]
Pituitary = tpm_with_gene[,z21]
Prostate = tpm_with_gene[,z22]
Salivary_Gland = Prostate = tpm_with_gene[,z23]
Skin = Prostate = tpm_with_gene[,z24]
Small_Intestine = Prostate = tpm_with_gene[,z25]
Spleen = Prostate = tpm_with_gene[,z26]
Stomach = Prostate = tpm_with_gene[,z27]
Testis = Prostate = tpm_with_gene[,z28]
Thyroid = Prostate = tpm_with_gene[,z29]
Uterus = Prostate = tpm_with_gene[,z30]
Vagina = Prostate = tpm_with_gene[,z31]

table(metadata$primary_site)


all.data_GTEx = cbind(Adipose_Tissue, Adrenal_Gland, Bladder, Blood, Blood_Vessel,
                      Bone_Marrow, Brain, Breast, Cervix_Uteri, Colon, Esophagus,  
                      Fallopian_Tube, Heart, Kidney, Liver, Lung, Muscle, Nerve, 
                      Ovary, Pancreas, Pituitary, Prostate, Salivary_Gland, Skin,           
                      Small_Intestine, Spleen, Stomach, Testis, Thyroid, Uterus,
                      Vagina)

rm(Adipose_Tissue, Adrenal_Gland, Bladder, Blood, Blood_Vessel,
   Bone_Marrow, Brain, Breast, Cervix_Uteri, Colon, Esophagus,  
   Fallopian_Tube, Heart, Kidney, Liver, Lung, Muscle, Nerve, 
   Ovary, Pancreas, Pituitary, Prostate, Salivary_Gland, Skin,           
   Small_Intestine, Spleen, Stomach, Testis, Thyroid, Uterus,
   Vagina)

library(limma)
mem.maxVSize(vsize = 1e12)
normalized.data = normalizeBetweenArrays(all.data_GTEx)

x = c(rep("Adipose_Tissue",515), rep("Adrenal_Gland",128), rep("Bladder",9), 
      rep("Blood",444), rep("Blood_Vessel",606), rep("Bone_Marrow",70), 
      rep("Brain",1152), rep("Breast",179), rep("Cervix_Uteri",10), rep("Colon",308),
      rep("Esophagus",653), rep("Fallopian_Tube",5), rep("Heart",377), rep("Kidney",28),
      rep("Liver",110), rep("Lung",288), rep("Muscle",396), rep("Nerve",278), rep("Ovary",88),
      rep("Pancreas",167), rep("Pituitary",107), rep("Prostate",85), rep("Salivary_Gland",55),
      rep("Skin",812), rep("Small_Intestine",92), rep("Spleen",100), rep("Stomach",174),
      rep("Testis",165), rep("Thyroid",279), rep("Uterus",78), rep("Vagina", 85)
)


### BOXPLOT ###
normalized.data <- as.data.frame(normalized.data)
selected_gene <- "IGSF3"
GENE = normalized.data[selected_gene,]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Adipose_Tissue",515), rep("Adrenal_Gland",128), rep("Bladder",9), 
           rep("Blood",444), rep("Blood_Vessel",606), rep("Bone_Marrow",70), 
           rep("Brain",1152), rep("Breast",179), rep("Cervix_Uteri",10), rep("Colon",308),
           rep("Esophagus",653), rep("Fallopian_Tube",5), rep("Heart",377), rep("Kidney",28),
           rep("Liver",110), rep("Lung",288), rep("Muscle",396), rep("Nerve",278), rep("Ovary",88),
           rep("Pancreas",167), rep("Pituitary",107), rep("Prostate",85), rep("Salivary_Gland",55),
           rep("Skin",812), rep("Small_Intestine",92), rep("Spleen",100), rep("Stomach",174),
           rep("Testis",165), rep("Thyroid",279), rep("Uterus",78), rep("Vagina", 85)
)


colnames(z) # 
colnames(z) <- c("gene", "type")

p2 = ggboxplot(z, x = "type", y = "gene",
              fill = "type", palette = color,
              width = 0.75,               # 适当调宽小提琴
              add = "boxplot",           # 在小提琴内部添加箱线图
              add.params = list(fill = "white", width = 0.15), # 设置内部箱线图为白色且变窄
              error.plot = "crossbar") + 
  theme_bw() +
  labs(x = "", y = "Gene Expression log2 (TPM + 1)", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), 
        axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 13), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  ggtitle(paste0(selected_gene, " levels (GTEx)"))

p2
# ggsave(filename = "IGSF3.png", plot = p2, width = 8, height = 5, dpi = 300)
###  stat_compare_means(aes(group = type), label.x = 1.2, label.y = 1) + p.value


########################### TCGA melanoma & GTEx ###############################
########################### TCGA melanoma & GTEx ###############################
########################### TCGA melanoma & GTEx ###############################
Skin = filter(metadata, primary_site == "Skin")
z24 = intersect(rownames(Skin), colnames(tpm_with_gene))
Skin = tpm_with_gene[,z24]

table(metadata$primary_site)
table(metadata$body_site_detail)

metadata[1:5, 1:5]


setorder(metadata, primary_site)
metadata_skin <- metadata
metadata_skin <- metadata[metadata$primary_site == "Skin", ]
table(metadata_skin$body_site_detail)

Skin <- as.data.frame(t(Skin))
Skin[1:5, 1:5]
metadata_skin[1:5, 1:5]

combined_tissue <- merge(metadata_skin, Skin, by = "row.names")

rownames(combined_tissue) <- combined_tissue$Row.names
combined_tissue <- combined_tissue[, -1] 

combined_tissue[1:6, 1:6]
table(combined_tissue$body_site_detail)
skin_tissue <- subset(combined_tissue, body_site_detail != "Cells - Transformed fibroblasts")

skin_tissue <- t(skin_tissue[,-c(1:5)])    # 去除5列metadata

tpm_gtex <- 2^(skin_tissue) - 0.001
tpm_gtex[tpm_gtex < 0] <- 0
tpm_gtex <- as.data.frame(log2(tpm_gtex + 1))


############ 读取melanoma tumor数据 ##############
setwd("/Users/naturexv/Datasets/TCGA/TCGA_SKCM_TPM")
melanoma <- fread("TCGA-SKCM.star_tpm.tsv.gz.tsv", header = T, sep = "\t", check.names = F)
clinical <- fread("TCGA-SKCM.clinical.tsv.gz.tsv", header = T, sep = "\t", check.names = F)
genesymbol <- fread("gencode.v36.annotation.gtf.gene.probemap", header = T, sep = "\t", check.names = F)

head(genesymbol)
tpm[1:5, 1:5]

tpm_with_symbol <- merge(melanoma, genesymbol[, .(id, gene)], 
                         by.x = "Ensembl_ID", 
                         by.y = "id", 
                         all.x = TRUE)
tpm_with_symbol <- distinct(.data = tpm_with_symbol, gene, .keep_all = T)     # 根据gene列去除重复基因
tpm_with_symbol <- column_to_rownames(.data = tpm_with_symbol, var = "gene")   
tpm_with_symbol <- tpm_with_symbol[,-1]

colnames(tpm_with_symbol)    # 473样本
clinical$sample    # 476样本
colnames(clinical)
table(clinical$sample_type.samples)

colnames(tpm_with_symbol) %in% clinical$sample 

sample_recorded <- intersect(clinical$sample, colnames(tpm_with_symbol))

tumor_tpm <- tpm_with_symbol[,sample_recorded]     # 有survival记录的sample

############## merge
common_genes <- intersect(rownames(tpm_gtex), rownames(tumor_tpm))
length(common_genes)

gtex_subset <- tpm_gtex[common_genes, ]
tumor_subset <- tumor_tpm[common_genes, ]

combined_data <- cbind(gtex_subset, tumor_subset)   # 合并
 
group <- c(rep("Normal", ncol(gtex_subset)), 
           rep("Tumor", ncol(tumor_subset)))

metadata <- data.frame(
  sample_id = colnames(combined_data),
  group = group,
  row.names = colnames(combined_data)
)

library(limma)

dge_matrix <- combined_data # 此时应为 log2(TPM+1)

# 进行分位数标准化（非常重要，让两组数据的分布强制一致）
dge_matrix_norm <- normalizeBetweenArrays(dge_matrix)

# 建立线性模型
design <- model.matrix(~ group)
fit <- lmFit(dge_matrix_norm, design)
fit <- eBayes(fit)

# 4. 提取结果
res <- topTable(fit, coef = 2, n = Inf)

res[which(res$logFC >= 1 & res$adj.P.Val < 0.05),'sig'] <- 'up'
res[which(res$logFC <= -1 & res$adj.P.Val < 0.05),'sig'] <- 'down'
res[which(abs(res$logFC) <= 1 | res$adj.P.Val >= 0.05),'sig'] <- 'non-significant'

# write.table(x = res, file = "DEGs_melanoma.csv", sep = ",")

### boxplot
dge_matrix_norm <- as.data.frame(dge_matrix_norm)
selected_gene <- "IGSF3"
GENE = dge_matrix_norm[selected_gene,]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框

z$type <- c(rep("Normal", ncol(gtex_subset)), 
            rep("Tumor", ncol(tumor_subset)))
colnames(z) <- c("gene", "type")

p <- ggboxplot(z, 
               x = "type",       # gender
               y = "gene", 
               fill = "type",    # gender
               palette = color, 
               add = "jitter", 
               add.params = list(size = 0.5, alpha = 0.4)) +
  stat_compare_means(method = "t.test", 
                     label = "p.format", # 选项: "p.format" 显示数值, "p.signif" 显示星号
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label.x = 1.25, label.y = 11) + 
  ylim(0, 12) +
  labs(x = "", y = "IGSF3 Expression log2(TPM + 1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13),
        axis.text.y = element_text(size = 13), 
        axis.title.y = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
p

ggsave("IGSF3_T_N.png", width = 3, height = 5, dpi = 300)





