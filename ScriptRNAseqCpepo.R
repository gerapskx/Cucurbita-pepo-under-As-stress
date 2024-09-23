###zucchini on Arsenic####

library("tidyverse")
library(readxl)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

CpepoAsCSV <- read.csv("C_pepo_RNASeq_counts.csv", row.names = 1)
Cpepodata<- read.csv("ColDataCpepo.csv", row.names = 1)
head(CpepoAsCSV)
head(Cpepodata)

all(colnames(CpepoAsCSV) %in% rownames(Cpepodata))

all(colnames(CpepoAsCSV) == rownames(Cpepodata))

dds <-DESeqDataSetFromMatrix(countData = round(CpepoAsCSV),
                       colData = Cpepodata,
                       design = ~Treatment)
dds

keep <- rowSums(counts(dds)) >= 10

dds2 <- dds[keep,]
dds2

dds2$Treatment <- relevel(dds$Treatment, ref = "Control")
dds2$Treatment

dds2 <- DESeq(dds2)

res <- results(dds2)
res
summary(res)

res0.05 <- results(dds2, alpha = 0.05)
summary(res0.05)

resultsNames(dds2)

###########50 uM####################

res1 <- results(dds2, contrast = c("Treatment", "As50", "Control"))

res1Tableraw <- write.csv(res1, file="As50Resultsraw.csv")

res1Tableraw <- read.csv("As50Resultsraw.csv")

As50volcano <- ggplot(res1Tableraw, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 1 & -log10(pvalue) > 2, 
                       ifelse(log2FoldChange < 0, "grey", "lightcoral"), 
                       "lightblue")),
    size = 3
  ) +
  scale_color_manual(
    values = c("lightblue", "grey", "lightcoral"),
    labels = c("Down-regulated", "Not significant", "Up-regulated")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"  # Change the legend title here
  ) +
  theme_minimal()
genes_to_labelAs50 <- c("Cp4.1LG01g02220", "Cp4.1LG13g01840", "Cp4.1LG17g09350","Cp4.1LG04g08100","Cp4.1LG01g07360","Cp4.1LG08g07390", "Cp4.1LG03g06560", "Cp4.1LG02g16350", "Cp4.1LG20g03710", "Cp4.1LG03g10370")

As50volcano <- As50volcano +
  geom_label_repel(data = res1Tableraw[res1Tableraw$X %in% genes_to_labelAs50, ],
            aes(label = X), size = 2, nudge_x = 0.1)  # Adjust nudge_y if needed
As50volcano

summary(res1)

res1_filtered <- subset(res1, padj < 0.05)

summary(res1_filtered)

res1_filtered <- subset(res1_filtered, log2FoldChange > 1 | log2FoldChange < -1)

summary(res1_filtered)

res1Table <- write.csv(res1_filtered, file="As50Results.csv")


###########100 uM####################

res2 <- results(dds2, contrast = c("Treatment", "As100", "Control"))

res2Tableraw <- write.csv(res2, file="As100Resultsraw.csv")

res2Tableraw <- read.csv("As100Resultsraw.csv")

As100volcano <- ggplot(res2Tableraw, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 1 & -log10(pvalue) > 2, 
                       ifelse(log2FoldChange < 0, "grey", "lightcoral"), 
                       "lightblue")),
    size = 3
  ) +
  scale_color_manual(
    values = c("lightblue", "grey", "lightcoral"),
    labels = c("Down-regulated", "Not significant", "Up-regulated")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"  # Change the legend title here
  ) +
  theme_minimal()

genes_to_labelAs100 <- c("Cp4.1LG08g07390",
                         "Cp4.1LG02g16350",
                         "Cp4.1LG03g06560",
                         "Cp4.1LG20g03710",
                         "Cp4.1LG03g10120",
                         "Cp4.1LG01g07360",
                         "Cp4.1LG04g08100",
                         "Cp4.1LG17g09350",
                         "Cp4.1LG08g00890",
                         "Cp4.1LG01g02220")

As100volcano <- As100volcano +
  geom_label_repel(data = res2Tableraw[res2Tableraw$X %in% genes_to_labelAs100, ],
                   aes(label = X), size = 2, nudge_x = 0.1)  # Adjust nudge_y if needed
As100volcano


summary(res2)

res2_filtered <- subset(res2, padj < 0.05)

summary(res2_filtered)

res2_filtered <- subset(res2_filtered, log2FoldChange > 1 | log2FoldChange < -1)

summary(res2_filtered)

res2Table <- write.csv(res2_filtered, file="As100Results.csv")


###########200 uM####################

res3 <- results(dds2, contrast = c("Treatment", "As200", "Control"))

res3Tableraw <- write.csv(res3, file="As200Resultsraw.csv")

res3Tableraw <- read.csv("As200Resultsraw.csv")

As200volcano <- ggplot(res3Tableraw, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 1 & -log10(pvalue) > 2, 
                       ifelse(log2FoldChange < 0, "grey", "lightcoral"), 
                       "lightblue")),
    size = 3
  ) +
  scale_color_manual(
    values = c("lightblue", "grey", "lightcoral"),
    labels = c("Down-regulated", "Not significant", "Up-regulated")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"  # Change the legend title here
  ) +
  theme_minimal()

genes_to_labelAs200 <- c("Cp4.1LG15g04660",
                         "Cp4.1LG08g07390",
                         "Cp4.1LG05g10590",
                         "Cp4.1LG03g10370",
                         "Cp4.1LG03g14960",
                         "Cp4.1LG19g11220",
                         "Cp4.1LG15g04120",
                         "Cp4.1LG18g02900",
                         "Cp4.1LG17g01770",
                         "Cp4.1LG05g03760")

As200volcano <- As200volcano +
  geom_label_repel(data = res3Tableraw[res3Tableraw$X %in% genes_to_labelAs200, ],
                   aes(label = X), size = 2, nudge_x = 0.1)  # Adjust nudge_y if needed
As200volcano

summary(res3)

res3_filtered <- subset(res3, padj < 0.05)

summary(res3_filtered)

res3_filtered <- subset(res3_filtered, log2FoldChange > 1 | log2FoldChange < -1)

summary(res3_filtered)

res3Table <- write.csv(res3_filtered, file="As200Results.csv")

###########pca plot############

vsd <- vst(dds2)

pcaplot <- plotPCA(vsd, intgroup = "Treatment", ntop = 500, returnData = TRUE)

pcaplotI <- plotPCA(vsd, intgroup = "Treatment", ntop = 500, returnData = FALSE)
pcaplotI

pcaplot$PC1 <- (colnames())
pcaplotenhanced <- ggplot(pcaplot, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() +
  theme_minimal() +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_blank()) +
  labs(x = "PC1: 77% variance", y = "PC2: 17% variance") +
  ggforce::geom_mark_ellipse(aes(fill = Treatment), alpha = 0.3)
  
desired_order <- c("Control", "As50", "As100", "As200")
pcaplot$Treatment <- factor(pcaplot$Treatment, levels = desired_order)

pcaplotenhanced
