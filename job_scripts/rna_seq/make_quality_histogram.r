
### Used to make histograms
#   USAGE: create_quality_histogram <quality_file> <name_file> [<output_file>]


args <- commandArgs(trailingOnly=T)

quality_file <- args[1]
name_conv <- args[2]

# set default name (and file path) the the base name if the quality file
title <- gsub("/?.*/", "", quality_file)
title <- gsub("_map_quality.txt", "", title)

print(title)
out <- paste0(title, ".png")
if (length(args) == 3) {   
    out <- args[3]
}

title <- gsub("/[^/]/", "", quality_file)
title <- gsub("_map_quality.txt", "", title)


qual <- read.table(quality_file, header=T, sep="\t", stringsAsFactors=F)
nc <- read.table(name_conv, header=T, row.names=1, sep="\t", stringsAsFactors=F, check.names=F)
qual[1:10, ]

row.names(nc) <- nc[, "contig"]
nc
qual$isolate <- "temp"
for (i in 1:nrow(qual)) {
    qual[i, "isolate"] <- nc[qual[i, "Genome"], "isolate"]
}

plotting <- qual[, names(qual) %in% c("quality", "isolate")]

names(plotting)
plotting[1:10, ]
library(ggplot2)
png(out)
p <- ggplot(plotting, aes(quality, fill=isolate)) + geom_histogram(alpha = 0.5, aes(y=..density..), position = 'identity') + ggtitle(title)
# ggsave(out, plot=p, width=7, height=7)
print(p)
dev.off()
