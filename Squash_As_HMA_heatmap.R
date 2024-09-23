library(ape)
library(ggplot2)
library(ComplexHeatmap)
library(readxl)

HMAtree <- ape::read.tree("newickHMAs")

dendogramHMA <- chronos(HMAtree)

dendogramHMA <- root(dendogramHMA, outgroup = "CpCuSOD", resolve.root = TRUE)

is.ultrametric(dendogramHMA)
is.rooted(dendogramHMA)

hc = as.hclust(dendogramHMA)

HMAproteins <- read_xlsx(path = "HMA-table.xlsx", sheet = 5)
HMAproteins <- as.data.frame(HMAproteins)

rownames(HMAproteins) <- HMAproteins[,1]
HMAproteins <- HMAproteins[,-1]
library(circlize)
Heatmap(matrix = HMAproteins, )
colors = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
circos.heatmap(HMAproteins, col=colors)
circos.clear()
circos.heatmap(HMAproteins, col = colors, cluster = dendogramHMA, dend.side = "inside", rownames.side = "outside")


