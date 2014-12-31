#!/usr/bin/env Rscript
library(RColorBrewer);
library(gplots);
library(vegan);
#library(gtools);

library(Heatplus)



MIN_READS <- 3000
#make palette
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100);



##read in the data file
myargs <- commandArgs(trailingOnly = TRUE);
myfile <- myargs[1];
print(myfile);

#this should be doable in the read statement but for some reason it isn't
data <- read.table(myfile, header = T, row.names = 1); 

#read in the ordering file
annotation <- read.table(myargs[2], sep="\t", header = TRUE, row.names = 1);
#order2 <- read.table(myargs[3], sep="\t", header = TRUE, row.names = 1);

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
#print(annotation)
order <- data.frame(order=annotation[, "order"])    #make a new dataframe to order
row.names(order) <- row.names(annotation)
annotation <- annotation[, -ncol(annotation)]       #remove the order column from annotation
print(order)
final <- merge(rare, order, all.y = T, by="row.names")

#final2 <- merge(rare, order2, all.y = T, by="row.names")

#reset the headers
row.names(final) <- final[, 1]
#row.names(final2) <- final2[, 1]
final <- final[, -1]
#final2 <- final2[, -1]

#print(names(final))

#reorder by the order column
final <- final[order(final[, "order"]), ]
#final2 <- final2[order(final2[, "order"]), ]

#replace NA with 0
final[is.na(final)] <- 0
#final2[is.na(final2)] <- 0

#remove order column
final <- final[, -ncol(final)]
#final2 <- final2[, -ncol(final2)]

#remove any isolates with no reads
final <- final[, colSums(final) > 0]
#final2 <- final2[, colSums(final2) > 0]

#convert to log scale
final <- final + 1
#final2 <- final2 + 1
final <- log(final)
#final2 <- log(final2)

#transpose to allow for labeling at the top
final <- t(final)
#final2 <- t(final2)



# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
#row_dist <- vegdist(final, method = "bray")
# Do average linkage hierarchical clustering. Other options are 'complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.
#row_clus <- hclust(row_dist, "aver")



# you have to transpose the dataset to get the samples as rows
#col_dist <- vegdist(t(final), method = "bray")
#col_clus <- hclust(col_dist, "aver")


png("heatmap.png", 1000, 1000)
heatmap <- annHeatmap2(final, 
dendrogram = list(Col = list(status="no"), Row = list(side = 4)),     #turn off column dendrograms
col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(59), breaks = 50,     #color options
legend = TRUE,
ann = list(Col = list(data = annotation)),
labels = list(Col = list(nrow = 12))
)
plot(heatmap)

# the annHeatmap2 function needs to be wrapped in the plot function in order to display the results
#plot(annHeatmap(final,
# the colours have to be fiddled with a little more than with the other functions. With 50 breaks in the heatmap densities, there needs to be 51 colours in the scale palette. Error messages from the plotting should help to determine how many colours are needed for a given number of breaks
#col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(51), breaks = 50,
# dendrogram controls how dendrograms are made and can be calculatated witin the function. The internal calculation can be overridden and externally calculated dendrograms can be displayed.
#dendrogram = list(Row = list(dendro = as.dendrogram(row_clus)), Col = list(dendro = as.dendrogram(col_clus))), legend = 3, # this puts the colour-scale legend on the plot. The number indicates the side on which to plot it (1 = bottom, 2 = left, 3 = top, 4 = right)
#labels = list(Col = list(nrow = 12)) # gives more space for the sample names
#))



# make a data frame of variables for annotation
# ann.dat <- data.frame(var1 = c(rep("cat1", 4), rep("cat2", 8)), var2 = rnorm(12,  mean = 50, sd = 20))

# plot(annHeatmap2(as.matrix(data.prop.1), col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(51), breaks = 50, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))), legend = 3, labels = list(Col = list(nrow = 12)), ann = list(Row = list(data = ann.dat))))
#heatmap <- heatmap.2(final, Colv="False", dendrogram = "row", col = scaleyellowred, margins = c(11, 5), xlab = "Isolates", ylab = "Samples", trace = "none", density.info="none", keysize = ".5", lhei = c(2, 15)) 
dev.off()

#png("basic_heatmap2.png", 1000, 1000)
#heatmap <- heatmap.2(final2, Colv="False", dendrogram = "row", col = scaleyellowred, margins = c(11, 5), xlab = "Isolates", ylab = "Samples", trace = "none", density.info = "none", keysize = ".5", lhei = c(2, 15)) 
#dev.off()

write.table(final, "table1.txt", sep="\t")
#write.table(final2, "table2.txt", sep="\t")
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


