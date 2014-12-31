#look into lattice: levelplot for making heatmaps



#
library(RColorBrewer);
library(gplots);
library(vegan);
#library(gtools);


MIN_READS <- 3000
#make palette
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100);



##read in the data file
myargs <- commandArgs(trailingOnly = TRUE);
myfile <- myargs[1];
print(myfile);

#this should be doable in the read statement but for some reason it isn't
#don't put sep="\t" and it works
data <- read.table(myfile, header = T, row.names = 1); 
#row.names(data) <- data[, 1]
#data <- data[, -1]


#read in the ordering file
order <- read.table(myargs[2], sep="\t", header = TRUE, row.names = 1);
order2 <- read.table(myargs[3], sep="\t", header = TRUE, row.names = 1);

#remove unmapped 
data <- data[-nrow(data), ]

#filter to number of reads
filtered <- data[, colSums(data) >= MIN_READS]

#rarefy data
rare <- rrarefy(t(filtered), MIN_READS)


#
##
###merge the data with the formatting frame
##
#


final <- merge(rare, order, all.y = T, by="row.names")
final2 <- merge(rare, order2, all.y = T, by="row.names")

#reset the headers
row.names(final) <- final[, 1]
row.names(final2) <- final2[, 1]
final <- final[, -1]
final2 <- final2[, -1]

#reorder by the order column
final <- final[order(final[, "order"]), ]
final2 <- final2[order(final2[, "order"]), ]

#replace NA with 0
final[is.na(final)] <- 0
final2[is.na(final2)] <- 0

#remove order column
final <- final[, -ncol(final)]
final2 <- final2[, -ncol(final2)]

#remove any isolates with no reads
final <- final[, colSums(final) > 0]
final2 <- final2[, colSums(final2) > 0]

#convert to log scale
final <- final + 1
final2 <- final2 + 1
final <- log(final)
final2 <- log(final2)

#transpose to allow for labeling at the top
final <- t(final)
final2 <- t(final2)


png("basic_heatmap1.png", 1000, 1000)
heatmap <- heatmap.2(final, Colv="False", dendrogram = "row", col = scaleyellowred, margins = c(11, 5), xlab = "Isolates", ylab = "Samples", trace = "none", density.info="none", keysize = ".5", lhei = c(2, 15)) 
dev.off()

png("basic_heatmap2.png", 1000, 1000)
heatmap <- heatmap.2(final2, Colv="False", dendrogram = "row", col = scaleyellowred, margins = c(11, 5), xlab = "Isolates", ylab = "Samples", trace = "none", density.info = "none", keysize = ".5", lhei = c(2, 15)) 
dev.off()

write.table(final, "table1.txt", sep="\t")
write.table(final2, "table2.txt", sep="\t")
#heatmap <- heatmap.2(as.matrix(switch.prop1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(11, 5), trace = "none", density.info = "none", xlab = "Isolates", ylab = "Samples", main = "Experiment A Heatmap", lhei = c(2, 8)) # this makes the colour-key legend a little thinner






#___________________________________________________________________________


#remove isolates with less than 1% abundance
#maxab <- apply(switch.prop, 2, max)
#n1 <- names(which(maxab < 0.01))
#switch.prop1 <- switch.prop[, -which(names(switch.prop) %in% n1)]

#
##cluster samples
#

# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
#switch.dist <- vegdist(switch.prop, method = "bray")
# Do average linkage hierarchical clustering. Other options are 'complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.
#row.clus <- hclust(switch.dist, "aver")

#
##cluster isolates
#
# you have to transpose the dataset to get the isolates as rows
#data.dist.g <- vegdist(t(data.prop.1), method = "bray")
#col.clus <- hclust(data.dist.g, "aver")

#transpose table
#switch <- t(data)
#switch <- data.frame(switch)
#colnames(switch) <- rownames(data)


#convert to proportions
#switch.prop <- data/rowSums(data)
#for (n in 1:ncol(data)) {
#    total <- sum(data[, n]);
#    if (total > 0) {
#        for (i in 1:nrow(data)) {
#            data[i, n] <- data[i, n] / total;
#        }
#    }
#}


