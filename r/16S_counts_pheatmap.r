library(pheatmap)
library(vegan)

args <- commandArgs(trailingOnly=T)
data <- read.table(args[1], sep="\t", header=T, row.names = 1, check.names=F)

#set minimum reads to be included, also convert to integer (from string)
MIN_READS <- strtoi(args[2])

#delete the taxonomy column
data <- data[, -ncol(data)]

#round the 16S normalized counts to whole numbers
data <- round(data, digits=0)

filtered <- data[, colSums(data) > MIN_READS]



rare <- rrarefy(t(filtered), MIN_READS)

#convert to rel abundance
relative <- rare / MIN_READS * 100
print(rowSums(relative))

#remove low-abundance samples
relative <- relative[order(rowSums(relative), decreasing=T), ]
relative <- relative[1:20, ]


#remove low-abundance OTUs
relative <- t(relative)
relative <- relative[order(rowSums(relative), decreasing=T), ]
final <- relative[1:10, ]

#rare <- rare[, rowSums(rare) > 500]

#log_rare <- log(final + 1)


palette <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)


#setEPS()
postscript("heatmap.ps")
pheatmap(as.matrix(final), cluster_rows=T, cluster_cols=F, col=palette, display_numbers=T, number_format="%.0f")
dev.off()
