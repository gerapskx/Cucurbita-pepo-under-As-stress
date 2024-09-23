library(tidyverse)
library(DESeq2)
library(readxl)

CpepoAsCSV <- read.csv("CpepoMetabolome.csv", row.names = 1)
Cpepodata<- read.csv("ColDataCpepo.csv", row.names = 1)
head(CpepoAsCSV)
head(Cpepodata)

all(colnames(CpepoAsCSV) %in% rownames(Cpepodata))

all(colnames(CpepoAsCSV) == rownames(Cpepodata))

dds <-DESeqDataSetFromMatrix(countData = round(CpepoAsCSV),
                             colData = Cpepodata,
                             design = ~Group)
dds

keep <- rowSums(counts(dds)) >= 10

dds2 <- dds[keep,]
dds2

#pca plot

vsd <- varianceStabilizingTransformation(dds2)

pcaplot <- plotPCA(vsd, intgroup = "Group", ntop = 500, returnData = TRUE)
pcaplot
desired_order <- c("Control", "As50", "As100", "As200")
pcaplot$group <- factor(pcaplot$group, levels = desired_order)

pcaplot$PC1 <- (colnames())

pcaplotenhanced <- ggplot(pcaplot, aes(x = PC1, y = PC2, color = Group)) +
  geom_point() +
  theme_minimal() +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_blank()) +
  labs(x = "PC1: 64.8% variance", y = "PC2: 14.3% variance") +
  ggforce::geom_mark_ellipse(aes(fill = Group), alpha = 0.3)

pcaplotenhanced




